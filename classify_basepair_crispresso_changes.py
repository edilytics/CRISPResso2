import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import upsetplot
import zipfile

from collections import defaultdict
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPResso2Align

# Set up logging
logger = logging.getLogger(__name__)

def get_refpos_dict(ref_aln_seq):
    #given a reference alignment this returns a dictionary such that refpos_dict[ind] = returns the index of the alignment corresponding to the ind'th base in the reference
    refpos_dict = {}
    ref_pos = 0
    for ind, ref_base in enumerate(ref_aln_seq):
        if ref_base != '-':
            refpos_dict[ref_pos] = ind
            ref_pos += 1
    return refpos_dict

def get_refpos_values(ref_aln_seq, read_aln_seq):
    """
    Given a reference alignment this returns a dictionary such that refpos_dict[ind] is the value of the read at the position corresponding to the ind'th base in the reference
    Any additional bases in the read (gaps in the ref) are assigned to the first position of the ref (i.e. refpos_dict[0])
    For other additional bases in the ref (gaps in the read), the value is appended to the last position of the ref that had a non-gap base (to the left)

    For example:
    ref_seq =  '--A-TGC-'
    read_seq = 'GGAGTCGA'
    get_refpos_values(ref_seq, read_seq)
    {0: 'GGAG', 1: 'T', 2: 'C', 3: 'GA'}

    Args:
    - ref_aln_seq: str, reference alignment sequence
    - read_aln_seq: str, read alignment sequence

    Returns:
    - refpos_dict: dict, dictionary such that refpos_dict[ind] is the value of the read at the position corresponding to the ind'th base in the reference

    """
    refpos_dict = defaultdict(str)

    # First, if there are insertions in read, add those to the first position in ref
    if ref_aln_seq[0] == '-':
        aln_index = 0
        read_start_bases = ""
        while aln_index < len(ref_aln_seq) and ref_aln_seq[aln_index] == '-':
            read_start_bases += read_aln_seq[aln_index]
            aln_index += 1
        refpos_dict[0] = read_start_bases
        ref_aln_seq = ref_aln_seq[aln_index:]
        read_aln_seq = read_aln_seq[aln_index:]
        
    ref_pos = 0
    last_nongap_ref_pos = 0
    for ind in range(len(ref_aln_seq)):
        ref_base = ref_aln_seq[ind]
        read_base = read_aln_seq[ind]
        if ref_base == '-':
            refpos_dict[last_nongap_ref_pos] += read_base
        else:
            refpos_dict[ref_pos] += read_base
            last_nongap_ref_pos = ref_pos
            ref_pos += 1
    return refpos_dict

def process_folder(crispresso_output_folder, output_root, wt_ref_name=None, target_ref_seq=None, target_ref_skip_allele_count=0, consider_changes_outside_of_guide=False, consider_indels_outside_of_guide=False):
    """
    Args:
    - crispresso_output_folder: str, path to CRISPResso output folder
    - output_root: str, root name for the output files
    - wt_ref_name: str, name of the reference sequence for WT
    - target_ref_seq: str, target reference sequence. If this value is not provided, the target ref will be inferred from the most-frequent non-reference MODIFIED allele in the allele table.
    - target_ref_skip_allele_count: int, number of alleles in the allele table to skip before assigning target reference sequence. This is useful when the most-frequent allele is not the target allele.
    - consider_changes_outside_of_guide: bool, if True, consider changes outside of the guide region. If this flag is not set, only reads with changes in the quantification window are considered for analysis.
    - consider_indels_outside_of_guide: bool, if True, consider indels outside of the guide region. If this flag is not set, only reads with indels in the quantification window are considered "Indel" reads.
    
    """
    logger.debug('Processing ' + crispresso_output_folder)
    crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)

    # extablish reference and target sequences
    if args.wt_ref_name is None and len(crispresso2_info['results']['refs']) == 1:
        wt_ref_name = list(crispresso2_info['results']['refs'].keys())[0]

    if wt_ref_name not in crispresso2_info['results']['refs']:
        possible_refs = list(crispresso2_info['results']['refs'].keys())
        raise Exception('Reference sequence for reference "' + wt_ref_name + '" not found - was a reference sequence named "' + wt_ref_name + '" provided? (possible refs are: ' + ", ".join(possible_refs) + ')')
    ref_seq  = crispresso2_info['results']['refs'][wt_ref_name]['sequence']

    if target_ref_seq is not None:
        target_seq = target_ref_seq
    else:
        logger.info('Target reference sequence not provided. Inferring from allele table.')
        z = zipfile.ZipFile(os.path.join(crispresso_output_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
        zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
        df_alleles = pd.read_csv(zf,sep="\t")
        target_seq = ""
        seen_nonref_allele_count = 0
        for idx, allele in df_alleles.iterrows():
            if allele.Aligned_Sequence.replace("-","") != ref_seq and allele.Read_Status == 'MODIFIED':
                if seen_nonref_allele_count >= target_ref_skip_allele_count:
                    target_seq = allele.Aligned_Sequence.replace("-","")
                    break
                else:
                    logger.debug('Skipping allele ' + str(idx) + ' with sequence ' + allele.Aligned_Sequence)
                seen_nonref_allele_count += 1
        zf.close()
        if target_seq == "":
            raise Exception('Target reference sequence not found in allele table (all reads were equal to the reference sequence)')

    if consider_indels_outside_of_guide:
        if not crispresso2_info['running_info']['args'].write_detailed_allele_table:
            raise Exception('To use parameter --consider_indels_outside_of_guide, CRISPResso run must be run with the parameter --write_detailed_allele_table')
    
    # create reference/target read alignment
    aln_gap_incentive = crispresso2_info['results']['refs'][wt_ref_name]['gap_incentive']
    aln_gap_open_arg = crispresso2_info['running_info']['args'].needleman_wunsch_gap_open
    aln_gap_extend_arg = crispresso2_info['running_info']['args'].needleman_wunsch_gap_extend

    aln_matrix_loc = crispresso2_info['running_info']['args'].needleman_wunsch_aln_matrix_loc
    if aln_matrix_loc == 'EDNAFULL':
        aln_matrix = CRISPResso2Align.make_matrix()
    else:
        if not os.path.exists(aln_matrix_loc):
            raise Exception('Alignment matrix file not found at ' + aln_matrix_loc)
        aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    aln_target_seq, aln_ref_seq, aln_score = CRISPResso2Align.global_align(target_seq, ref_seq, matrix=aln_matrix, gap_incentive=aln_gap_incentive, gap_open=aln_gap_open_arg, gap_extend=aln_gap_extend_arg)

    logger.debug('Aligned target:    ' + aln_target_seq)
    logger.debug('Aligned reference: ' + aln_ref_seq)
           
    # get indices of reference sequence to include in analysis
    if consider_changes_outside_of_guide:
        ref_positions_to_include = [x for x in range(len(ref_seq))]
    else:
        ref_positions_to_include = crispresso2_info['results']['refs'][wt_ref_name]['include_idxs']

    # discover positions and bases that are different between reference and target
    ref_changes_dict = get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_changes_arr = []
    for idx in ref_positions_to_include:
        ref_base = ref_seq[idx]
        if ref_changes_dict[idx] != ref_base:
            bp_changes_arr.append((idx, ref_base, ref_changes_dict[idx]))

    logger.debug('Found ' + str(len(bp_changes_arr)) + ' base changes: ' + str(bp_changes_arr))

    # open allele table
    z = zipfile.ZipFile(os.path.join(crispresso_output_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf,sep="\t")


    # set up counters
    binary_allele_counts = defaultdict(int) # e.g. T,T,X,T > 100 where each item is a string of the base at each position in bp_changes_arr, and 'X' is nontarget
    category_allele_counts = defaultdict(int) # e.g. T,T,R,T > 100 where each item is a string of the base at each position in bp_changes_arr, and 'T' is Target, 'R' is Reference, 'D' is Deletion, 'I' is insertion, and 'N' is anything else
    precise_allele_counts = defaultdict(int) # e.g. A,A,C,AA > 100 where each item is a string of the base at each position in bp_changes_arr

    total_alleles = 0
    total_alleles_reads = 0
    total_alleles_on_ref = 0
    total_alleles_reads_on_ref = 0

    total_target_noindel_reads = 0
    total_target_indel_reads = 0
    total_reference_noindel_reads = 0
    total_reference_indel_reads = 0
    total_other_noindel_reads = 0
    total_other_indel_reads = 0

    target_base_counts = [0] * len(bp_changes_arr)
    reference_base_counts = [0] * len(bp_changes_arr)
    deletion_base_counts = [0] * len(bp_changes_arr)
    insertion_base_counts = [0] * len(bp_changes_arr)
    other_base_counts = [0] * len(bp_changes_arr)

    # iterate all alleles in input allele table
    for idx, allele in df_alleles.iterrows():
        total_alleles += 1
        total_alleles_reads += allele['#Reads']

        if allele.Reference_Name != wt_ref_name:
            continue
        total_alleles_on_ref += 1
        total_alleles_reads_on_ref += allele['#Reads']


        has_indel_guide = False
        if allele.n_deleted > 0:
            has_indel_guide = True
        if allele.n_inserted > 0:
            has_indel_guide = True

        has_indel = has_indel_guide
        if consider_indels_outside_of_guide:
            has_any_indel = False
            if allele.all_deletion_positions != '[]': # need to run with detailed allele output to get this value
                has_any_indel = True
            if allele.all_insertion_positions != '[]':
                has_any_indel = True

        ref_aln = allele.Reference_Sequence
        read_aln = allele.Aligned_Sequence
        ref_base_position_lookup = get_refpos_values(ref_aln, read_aln)

        binary_arr = []
        cat_arr = []
        val_arr = []
        for ind, (ref_ind, ref_base, mod_base) in enumerate(bp_changes_arr):
            base_at_pos = ref_base_position_lookup[ref_ind]
            this_binary = 'X'
            this_category = 'N'
            if base_at_pos == ref_base:
                this_category = 'R'
                reference_base_counts[ind] += allele['#Reads']
            elif base_at_pos == mod_base:
                this_category = 'T'
                this_binary = 'T'
                target_base_counts[ind] += allele['#Reads']
            elif base_at_pos == '-':
                this_category = 'D'
                deletion_base_counts[ind] += allele['#Reads']
            elif len(base_at_pos) != 1:
                this_category = 'I'
                insertion_base_counts[ind] += allele['#Reads']
            else:
                this_category = 'N'
                other_base_counts[ind] += allele['#Reads']
            binary_arr.append(this_binary)
            cat_arr.append(this_category)
            val_arr.append(base_at_pos)

        if cat_arr.count('R') == len(cat_arr):
            if not has_indel:
                total_reference_noindel_reads += allele['#Reads']
            else:
                total_reference_indel_reads += allele['#Reads']
        elif cat_arr.count('T') == len(cat_arr):
            if not has_indel:
                total_target_noindel_reads += allele['#Reads']
            else:
                total_target_indel_reads += allele['#Reads']
        else:
            if not has_indel:
                total_other_noindel_reads += allele['#Reads']
            else:
                total_other_indel_reads += allele['#Reads']
        
        binary_arr_str = "\t".join(binary_arr) + "\t" + str(has_indel)
        cat_arr_str = "\t".join(cat_arr) + "\t" + str(has_indel)
        val_arr_str = "\t".join(val_arr) + "\t" + str(has_indel)

        binary_allele_counts[binary_arr_str] += allele['#Reads']
        category_allele_counts[cat_arr_str] += allele['#Reads']
        precise_allele_counts[val_arr_str] += allele['#Reads']


    total_counts = [total_alleles_reads] * len(bp_changes_arr)

    with open(output_root + '.binary_allele_counts.txt','w') as fout:
        sorted_binary_allele_counts = sorted(binary_allele_counts.keys(), key=lambda x: binary_allele_counts[x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_changes_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_binary_allele_counts:
            fout.write(allele_str + '\t' + str(binary_allele_counts[allele_str]) + '\n')


    with open(output_root + '.category_allele_counts.txt','w') as fout:
        sorted_category_allele_counts = sorted(category_allele_counts.keys(), key=lambda x: category_allele_counts[x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_changes_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_category_allele_counts:
            fout.write(allele_str + '\t' + str(category_allele_counts[allele_str]) + '\n')
    
    with open(output_root + '.precise_allele_counts.txt','w') as fout:
        sorted_precise_allele_counts = sorted(precise_allele_counts.keys(), key=lambda x: precise_allele_counts[x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_changes_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_precise_allele_counts:
            fout.write(allele_str + '\t' + str(precise_allele_counts[allele_str]) + '\n')
            
    with open(output_root + '.arrays.txt','w') as fout:
        fout.write('Class\t' + "\t".join([str(x) for x in bp_changes_arr]) + '\n')
        fout.write('total_counts\t' + "\t".join([str(x) for x in total_counts]) + '\n')
        fout.write('reference_counts\t' + "\t".join([str(x) for x in reference_base_counts]) + '\n')
        fout.write('target_counts\t' + "\t".join([str(x) for x in target_base_counts]) + '\n')
        fout.write('deletion_counts\t' + "\t".join([str(x) for x in deletion_base_counts]) + '\n')
        fout.write('insertion_counts\t' + "\t".join([str(x) for x in insertion_base_counts]) + '\n')
        fout.write('other_counts\t' + "\t".join([str(x) for x in other_base_counts]) + '\n')

    with open(output_root + '.counts.txt','w') as fout:
        target_name = 'Target'
        fout.write("\t".join([wt_ref_name,wt_ref_name+"_indels",target_name,target_name+"_indels","other","other_indels"]) + '\n')
        fout.write("\t".join([str(x) for x in [total_reference_noindel_reads, total_reference_indel_reads, total_target_noindel_reads, total_target_indel_reads, total_other_noindel_reads, total_other_indel_reads]]) + '\n')


    logger.debug('Read ' + str(total_alleles) + ' alleles with ' + str(total_alleles_reads) + ' reads')
    logger.debug('Got ' + str(total_alleles_on_ref) + ' alleles on reference "' + wt_ref_name + '" with ' + str(total_alleles_reads_on_ref) + ' reads')

    if len(bp_changes_arr) > 0:
        plot_successful_counts_by_category(output_root, bp_changes_arr, reference_base_counts, target_base_counts, deletion_base_counts, insertion_base_counts, other_base_counts)
        plot_combination_upset(output_root, bp_changes_arr, binary_allele_counts)
    else:
        logger.error('No base changes found in ' + crispresso_output_folder + '. Not plotting results.')

    return({'wt': total_reference_noindel_reads, 'wt_indels': total_reference_indel_reads, 'target': total_target_noindel_reads, 'target_indels': total_target_indel_reads, 'other': total_other_noindel_reads, 'other_indels': total_other_indel_reads})


def plot_successful_counts_by_category(fig_root, bp_changes_arr, reference_base_counts, target_base_counts, deletion_base_counts, insertion_base_counts, other_base_counts):
    """
    Plot barplot showing edit status by category
    
    Args:
    - fig_root: str, root name for the figure file
    - bp_changes_arr: list of tuples, each tuple contains the index of the base change, the reference base, and the target edit base
    - reference_base_counts: list of int, counts of reference bases
    - target_base_counts: list of int, counts of target base changes
    - deletion_base_counts: list of int, counts of deletions
    - insertion_base_counts: list of int, counts of insertions
    - other_base_counts: list of int, counts of other bases
    """
    # Indices for the x-axis
    indices = np.arange(len(bp_changes_arr))

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 6)) 

    # Stacked bar plot
    ax.bar(indices, reference_base_counts, label='Reference')
    bottom_so_far = reference_base_counts

    ax.bar(indices, target_base_counts, label='Target', bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, target_base_counts)]

    ax.bar(indices, deletion_base_counts, label='Deletion', bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, deletion_base_counts)]

    ax.bar(indices, insertion_base_counts, label='Insertion', bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, insertion_base_counts)]

    ax.bar(indices, other_base_counts, label='Other', bottom=bottom_so_far)

    # Adding labels and title
    ax.set_xlabel('Base Index')
    ax.set_ylabel('Read Counts')
    ax.set_title('Base composition by category')

    ax.set_xticks(indices)
    ax.set_xticklabels([str(x[0]) for x in bp_changes_arr])
    ax.legend()

    # Display the plot
    plt.savefig(fig_root + '.base_composition_by_category.pdf')

def plot_combination_upset(fig_root, bp_changes_arr, binary_allele_counts):
    header_arr = []
    for ind, (ref_ind, ref_base, mod_base) in enumerate(bp_changes_arr):
        header_arr.append(str(ref_ind) + ':' + ref_base + '->' + mod_base)
    header_arr.append('has_indel')
    header_arr.append('cat_counts')

    df_by_combination_items = []
    for ref_comb in binary_allele_counts:
        arr_for_upset = ref_comb.split("\t")
        for i in range(len(arr_for_upset)-1): # don't check the last column (has indel)
            arr_for_upset[i] = arr_for_upset[i] == 'T'
        arr_for_upset[len(arr_for_upset)-1] = arr_for_upset[len(arr_for_upset)-1] == 'True'

        arr_for_upset.append(binary_allele_counts[ref_comb])
        df_by_combination_items.append(arr_for_upset)

    df_by_combination = pd.DataFrame(df_by_combination_items, columns=header_arr)
    df_by_combination.set_index(header_arr[:-1], inplace=True)
    ax_dict = upsetplot.UpSet(df_by_combination.cat_counts, show_counts=True, sort_categories_by='-input').plot()
    plt.title(fig_root + ' read counts with no indels')
    plt.savefig(fig_root + '.by_amplicon_combination.no_indels.pdf')
    plt.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Produce summary of editing frequencies for CRISPResso outputs')
    parser.add_argument('--crispresso_folder', type=str, help='CRISPResso output folder', required=True)
    parser.add_argument('--output_root', type=str, help='Root name for the output files', default=None)
    parser.add_argument('--wt_ref_name', type=str, help='Name of the reference sequence for WT', default=None)
    parser.add_argument('--target_ref_seq', type=str, help='Target reference sequence if it is not in the CRISPResso output. If this value is not provided, the target ref will be inferred from the most-frequent non-reference MODIFIED allele in the allele table)', default=None)
    parser.add_argument('--target_ref_skip_allele_count', type=int, help='Number of alleles in the allele table to skip before assigning target reference sequence. This is useful when the most-frequent allele is not the target allele.', default=0)
    parser.add_argument('--consider_changes_outside_of_guide', action='store_true', help='Consider changes outside of the guide region. If this flag is not set, only reads with changes in the quantification window are considered for analysis.')
    parser.add_argument('--consider_indels_outside_of_guide', action='store_true', help='Consider indels outside of the guide region. If this flag is not set, only reads with indels in the quantification window are considered "Indel" reads.')
    parser.add_argument('--debug', action='store_true', help='Print debug messages')
    args = parser.parse_args()

    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # console logger
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)



    if args.crispresso_folder.endswith('/'):
        args.crispresso_folder = args.crispresso_folder[:-1]

    crispresso_folders_to_analyze = []
    crispresso_folder_infos = {} # keep track of folders in this batch
    crispresso_folder_output_roots = {}

    crispresso_info_file = os.path.join(args.crispresso_folder, 'CRISPResso2_info.json')
    if os.path.exists(crispresso_info_file):
        try:
            run_data = CRISPRessoShared.load_crispresso_info(args.crispresso_folder)
            crispresso_folders_to_analyze.append(args.crispresso_folder)
            crispresso_folder_infos[args.crispresso_folder] = run_data
            if args.output_root is not None:
                crispresso_folder_output_roots[args.crispresso_folder] = args.output_root
            else:
                crispresso_folder_output_roots[args.crispresso_folder] = args.crispresso_folder + ".classify"
        except Exception as e:
            logger.error('Could not open CRISPResso2 info file in ' + args.crispresso_folder)

    batch_info_file = os.path.join(args.crispresso_folder, 'CRISPResso2Batch_info.json')
    if os.path.exists(batch_info_file):
        batch_data = CRISPRessoShared.load_crispresso_info(
            args.crispresso_folder, 'CRISPResso2Batch_info.json',print('set here 2')
        )
        if 'completed_batch_arr' in batch_data['results']:
            run_names = batch_data['results']['completed_batch_arr']
            for run_name in run_names:
                run_folder_loc = os.path.join(args.crispresso_folder, 'CRISPResso_on_%s'%run_name)
                try:
                    run_data = CRISPRessoShared.load_crispresso_info(run_folder_loc)
                    crispresso_folder_infos[run_folder_loc] = run_data
                    crispresso_folders_to_analyze.append(run_folder_loc)
                    if args.output_root is not None:
                        crispresso_folder_output_roots[run_folder_loc] = args.output_root + '.' + run_name
                    else:
                        crispresso_folder_output_roots[run_folder_loc] = args.crispresso_folder + ".classify." + run_name
                except Exception as e:
                    logger.error('Could not open CRISPResso2 info file in ' + run_folder_loc)
        else:
            logger.error('Could not process batch folder ' + args.crispresso_folder)

    if args.output_root is None:
        summary_output_file = args.crispresso_folder + ".classify.summary.txt"
    else:
        summary_output_file = args.output_root + ".summary.txt"
    
    with open(summary_output_file,'w') as fout:
        fout.write('Folder\tWT\tWT_indels\tTarget\Target_indels\tOther\tOther_indels\n')
        for folder in crispresso_folders_to_analyze:
            folder_category_counts = process_folder(folder, crispresso_folder_output_roots[folder], wt_ref_name=args.wt_ref_name, target_ref_seq=args.target_ref_seq, 
                                                    target_ref_skip_allele_count=args.target_ref_skip_allele_count, consider_changes_outside_of_guide=args.consider_changes_outside_of_guide,
                                                    consider_indels_outside_of_guide=args.consider_indels_outside_of_guide)
            fout.write(folder + '\t' + '\t'.join([str(folder_category_counts[x]) for x in ['wt','wt_indels','target','target_indels','other','other_indels']]) + '\n')
    
    logger.info('Finished analyzing ' + str(len(crispresso_folders_to_analyze)) + ' CRISPResso runs. Wrote summary to ' + summary_output_file)