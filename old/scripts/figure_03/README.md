# Figure 3 scripts

Scripts for untreated malignant-cell plasticity, CNV sharing/burden visualization, and per-patient CANDLE plasticity summaries.

Recommended run order:
1. `01_analyze_untreated_subcompartment_monocle_trajectory.R`
2. `02_plot_cnv_plasticity_compact_panels.R`
3. `03_plot_per_patient_candle_summary.R`

The trajectory script uses blastema as the root state. The CNV script assumes that upstream CNV summary tables have already been generated.
