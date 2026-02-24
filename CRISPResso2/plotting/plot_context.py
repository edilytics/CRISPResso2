"""PlotContext: read-only view into CRISPRessoCORE's computed data.

This is the ONLY part of the plugin architecture that lives in CRISPResso2.
CRISPRessoPro and plugins consume this to generate custom plots.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional


@dataclass
class PlotContext:
    """Read-only view into CORE's computed data, passed to CRISPRessoPro.

    Constructed once per run by CORE after the analysis and built-in plot
    phases complete. Wraps references to CORE's local variables -- no data
    is copied.

    Plugin data_funcs should treat all fields as read-only. Modifying the
    underlying data structures will corrupt the analysis results.
    """

    # === Run-level data (always available) ===

    args: object                    # argparse.Namespace -- all CLI arguments
    run_data: dict                  # crispresso2_info dict
    refs: dict                      # Per-amplicon reference data
    ref_names: list                 # Ordered list of amplicon/reference names

    # Alignment counts (ref_name -> int)
    counts_total: dict
    counts_modified: dict
    counts_unmodified: dict
    counts_discarded: dict
    counts_insertion: dict
    counts_deletion: dict
    counts_substitution: dict

    # Class counts (class_name -> count)
    class_counts: dict
    N_TOTAL: int

    # Allele data
    df_alleles: object              # pd.DataFrame

    # Per-position count vectors (ref_name -> numpy array)
    all_insertion_count_vectors: dict
    all_insertion_left_count_vectors: dict
    all_deletion_count_vectors: dict
    all_substitution_count_vectors: dict
    all_indelsub_count_vectors: dict
    all_substitution_base_vectors: dict
    all_base_count_vectors: dict

    # Quantification window vectors
    insertion_count_vectors: dict
    deletion_count_vectors: dict
    substitution_count_vectors: dict
    insertion_length_vectors: dict
    deletion_length_vectors: dict

    # Nucleotide frequency data (ref_name -> DataFrame)
    nucleotide_frequency_summary: dict
    nucleotide_percentage_summary: dict

    # Histograms (ref_name -> Counter)
    hists_frameshift: dict
    hists_inframe: dict

    # Frameshift / coding counts (ref_name -> int)
    counts_modified_frameshift: dict
    counts_modified_non_frameshift: dict
    counts_non_modified_non_frameshift: dict
    counts_splicing_sites_modified: dict

    # HDR / ref1-aligned vectors (populated when expected_hdr_amplicon_seq is set)
    ref1_all_insertion_count_vectors: dict = field(default_factory=dict)
    ref1_all_deletion_count_vectors: dict = field(default_factory=dict)
    ref1_all_substitution_count_vectors: dict = field(default_factory=dict)
    ref1_all_indelsub_count_vectors: dict = field(default_factory=dict)
    ref1_all_insertion_left_count_vectors: dict = field(default_factory=dict)
    ref1_all_base_count_vectors: dict = field(default_factory=dict)

    # Configuration
    custom_config: dict = field(default_factory=dict)

    # Utility
    _jp: Optional[Callable] = None  # Joins path with output directory
    save_png: bool = False
    output_directory: str = ""

    # === Scope fields (set by Pro during iteration) ===

    ref_name: Optional[str] = None
    sgRNA_ind: Optional[int] = None
