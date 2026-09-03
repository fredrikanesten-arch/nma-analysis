## Running the updated NMA workflow

Build the enriched long dataset (adds treatment/class/bridge flags):

```r
Rscript /home/runner/work/nma-analysis/nma-analysis/build_mmc5_ms_smd_bias_adj_dataset.R
```

Run the multi-analysis netmeta workflow (primary hybrid + sensitivities):

```r
Rscript /home/runner/work/nma-analysis/nma-analysis/netmeta_class_ms_control_arm_baseline_sd.R
```

Run a single analysis mode:

```r
Rscript /home/runner/work/nma-analysis/nma-analysis/netmeta_class_ms_control_arm_baseline_sd.R <input_csv> <mapping_csv> <output_dir> Placebo drop hybrid any_ad_any_psychotherapy
```
