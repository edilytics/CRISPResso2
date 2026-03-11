"""PlotContext: view into CRISPRessoCORE's computed data for plugins.

This is the ONLY part of the plugin architecture that lives in CRISPResso2.
CRISPRessoPro and plugins consume this to generate custom plots.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import pandas as pd


@dataclass
class PlotContext:
    """View into CORE's computed data, passed to CRISPRessoPro and plugins.

    Constructed once per run by CORE after the analysis and built-in plot
    phases complete. Wraps references to CORE's local variables -- no data
    is copied.

    .. warning::

        All fields except ``ref_name``, ``sgRNA_ind``, and
        ``coding_seq_ind`` should be treated as read-only. Mutating the
        underlying dicts, arrays, or DataFrames will corrupt the analysis
        results and the output report. Python cannot enforce this at
        runtime -- it is a contract.
    """

    # === Run-level data (always available) ===

    args: argparse.Namespace        # All CLI arguments
    run_data: dict                  # crispresso2_info dict
    refs: dict                      # Per-amplicon reference data
    ref_names: list[str]            # Ordered list of amplicon/reference names

    # Alignment counts (ref_name -> int)
    counts_total: dict[str, int]
    counts_modified: dict[str, int]
    counts_unmodified: dict[str, int]
    counts_discarded: dict[str, int]
    counts_insertion: dict[str, int]
    counts_deletion: dict[str, int]
    counts_substitution: dict[str, int]

    # Class counts (class_name -> count)
    class_counts: dict[str, int]
    N_TOTAL: int

    # Allele data
    df_alleles: pd.DataFrame

    # Per-position count vectors (ref_name -> numpy array)
    all_insertion_count_vectors: dict[str, np.ndarray]
    all_insertion_left_count_vectors: dict[str, np.ndarray]
    all_deletion_count_vectors: dict[str, np.ndarray]
    all_substitution_count_vectors: dict[str, np.ndarray]
    all_indelsub_count_vectors: dict[str, np.ndarray]
    all_substitution_base_vectors: dict[str, np.ndarray]
    all_base_count_vectors: dict[str, np.ndarray]

    # Quantification window vectors
    insertion_count_vectors: dict[str, np.ndarray]
    deletion_count_vectors: dict[str, np.ndarray]
    substitution_count_vectors: dict[str, np.ndarray]
    insertion_length_vectors: dict[str, np.ndarray]
    deletion_length_vectors: dict[str, np.ndarray]

    # Histograms (ref_name -> Counter)
    hists_frameshift: dict
    hists_inframe: dict

    # Frameshift / coding counts (ref_name -> int)
    counts_modified_frameshift: dict[str, int]
    counts_modified_non_frameshift: dict[str, int]
    counts_non_modified_non_frameshift: dict[str, int]
    counts_splicing_sites_modified: dict[str, int]

    # === Fields with defaults (must follow all required fields) ===

    # Nucleotide frequency data (ref_name -> DataFrame)
    # Not populated by CRISPRessoCORE (always {}); reserved for Batch/Aggregate.
    nucleotide_frequency_summary: dict = field(default_factory=dict)
    nucleotide_percentage_summary: dict = field(default_factory=dict)

    # HDR / ref1-aligned vectors (populated when expected_hdr_amplicon_seq is set)
    ref1_all_insertion_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_deletion_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_substitution_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_indelsub_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_insertion_left_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_base_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)

    # Configuration
    custom_config: dict = field(default_factory=dict)

    # Utility
    _jp: Optional[Callable[[str], str]] = None  # Joins path with output directory
    save_png: bool = False
    output_directory: str = ""

    # === Additional data fields (populated by CORE when available) ===

    # Allele classification ordering (for plot_1b/1c)
    class_counts_order: list[str] = field(default_factory=list)

    # Allele homology data (for plot_1e)
    homology_scores: list[float] = field(default_factory=list)
    homology_counts: list[int] = field(default_factory=list)

    # Noncoding modification vectors (for plot_7; ref_name → numpy array)
    insertion_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)
    deletion_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)
    substitution_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)

    # Quantification-window substitution base vectors (for plot_10c)
    # Same shape as all_substitution_base_vectors but window-only
    substitution_base_vectors: dict[str, np.ndarray] = field(default_factory=dict)

    # Scaffold insertion sizes DataFrame (for plot_11c)
    df_scaffold_insertion_sizes: Optional[pd.DataFrame] = None

    # Alternate allele counts (for plots 10b/10c; populated when base_editor_output)
    # ref_name → {ref_nuc: {obs_nuc: count}}
    alt_nuc_counts: dict[str, dict] = field(default_factory=dict)       # quantification window
    alt_nuc_counts_all: dict[str, dict] = field(default_factory=dict)   # full amplicon

    # === Scope fields (set by Pro during iteration) ===

    ref_name: Optional[str] = None
    sgRNA_ind: Optional[int] = None
    coding_seq_ind: Optional[int] = None
