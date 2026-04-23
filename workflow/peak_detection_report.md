---
title: "Peak Detection Across Project Datasets"
knit: !expr function(input_file, encoding) rmarkdown::render(input_file, encoding = encoding, output_dir = dirname(input_file), intermediates_dir = dirname(input_file), envir = new.env(parent = globalenv()))
output:
  html_document:
    toc: true
    toc_float: true
date: "2026-04-20"
---



# Goal

This report applies the promoted density-based peak detection workflow to each
project dataset under `data/`, using the gated DNA vectors exported under
`processed_data/`. For every dataset it detects one below-threshold reference
peak, up to two above-threshold peaks, and flags the primary outlier sample
using the same two-peak ratio rule that was developed in `dev/`.





# Dataset Overview


|dataset_id                       |dataset_name                     | n_samples| n_with_cen_peak| n_with_g1_peak| n_with_two_tumor_peaks| n_primary_outliers| median_ratio_above_below| median_count_ratio_above_below|
|:--------------------------------|:--------------------------------|---------:|---------------:|--------------:|----------------------:|------------------:|------------------------:|------------------------------:|
|anoxia-flowcytometry             |Anoxia_FlowCytometry             |        24|              24|             24|                     24|                  1|                  166.062|                         30.805|
|hypoxia-sum159                   |Hypoxia_SUM159                   |        20|              20|             20|                     20|                  1|                   35.807|                          2.695|
|polyploidization-ethanolfixation |Polyploidization_EthanolFixation |        21|              21|             21|                     21|                  1|                  417.609|                          5.462|

# Saved Outputs


|dataset_name                     |dataset_id                       |cache_source            |summary_csv                                                                                                                   |
|:--------------------------------|:--------------------------------|:-----------------------|:-----------------------------------------------------------------------------------------------------------------------------|
|Anoxia_FlowCytometry             |anoxia-flowcytometry             |existing_processed_data |C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/anoxia-flowcytometry/peak_detection_summary.csv             |
|Hypoxia_SUM159                   |hypoxia-sum159                   |existing_processed_data |C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/hypoxia-sum159/peak_detection_summary.csv                   |
|Polyploidization_EthanolFixation |polyploidization-ethanolfixation |existing_processed_data |C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/polyploidization-ethanolfixation/peak_detection_summary.csv |

# Per-Dataset Results

## Anoxia_FlowCytometry

Peak input source: `existing_processed_data`  
Cached summary: `C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/anoxia-flowcytometry/peak_detection_summary.csv`

Samples: 24. Two-peak calls: 24. Primary outliers: 1.


<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-1.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>
<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-2.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>

## Hypoxia_SUM159

Peak input source: `existing_processed_data`  
Cached summary: `C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/hypoxia-sum159/peak_detection_summary.csv`

Samples: 20. Two-peak calls: 20. Primary outliers: 1.


<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-3.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>
<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-4.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>

## Polyploidization_EthanolFixation

Peak input source: `existing_processed_data`  
Cached summary: `C:/Users/4473331/Documents/projects/023_ploidyFlow/processed_data/polyploidization-ethanolfixation/peak_detection_summary.csv`

Samples: 21. Two-peak calls: 21. Primary outliers: 1.


<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-5.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>
<div class="figure" style="text-align: center">
<img src="figure/dataset-sections-6.png" alt="plot of chunk dataset-sections" width="100%" />
<p class="caption">plot of chunk dataset-sections</p>
</div>
