---
title: "Flow Gating Workflow"
output:
  html_document:
    toc: true
    toc_float: true
date: "2026-04-20"
---



# Goal

This document runs the gating workflow independently for each dataset directory
under `data/`. Each dataset gets its own cache under `processed_data/` and its
own figure namespace under `figure/`, so additional datasets can be added
without colliding with existing outputs. The current inventory includes the
`Polyploidization_EthanolFixation` dataset alongside the anoxia and hypoxia
collections.


<script src="https://unpkg.com/panzoom@9.4.0/dist/panzoom.min.js"></script>
<style>
.plot-gallery {
  margin: 1rem 0 2rem 0;
}
.plot-gallery-controls {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  margin-bottom: 0.75rem;
}
.plot-gallery-controls button {
  padding: 0.3rem 0.7rem;
}
.plot-gallery-status {
  font-size: 0.95rem;
  color: #444;
}
.plot-slide {
  max-width: 100%;
}
.plot-caption {
  margin-top: 0.5rem;
  font-size: 0.95rem;
  color: #444;
  text-align: center;
}
</style>





# Dataset Inventory


|dataset_name                     |dataset_id                       | n_files|has_cache |has_all_reusable_plots |
|:--------------------------------|:--------------------------------|-------:|:---------|:----------------------|
|Anoxia_FlowCytometry             |anoxia-flowcytometry             |      24|TRUE      |TRUE                   |
|Hypoxia_SUM159                   |hypoxia-sum159                   |      20|TRUE      |TRUE                   |
|Polyploidization_EthanolFixation |polyploidization-ethanolfixation |      21|TRUE      |TRUE                   |

# Hypoxia Passaging Days


|sample_name                 |latest_matching_id           |latest_match_date | relative_day|
|:---------------------------|:----------------------------|:-----------------|------------:|
|Sample_SUM159_2N_C_A1_.fcs  |SUM-159_NLS_2N_C_A1_seedT1   |2024-11-18        |            0|
|Sample_SUM159_4N_C_A1_.fcs  |SUM-159_NLS_4N_C_A1_seedT1   |2024-11-18        |            0|
|Sample_SUM159_2N_O1_A6.fcs  |SUM-159_NLS_2N_O1_A6_seedT2  |2024-12-23        |           35|
|Sample_SUM159_2N_O2_A6.fcs  |SUM-159_NLS_2N_O2_A6_seedT2  |2024-12-23        |           35|
|Sample_SUM159_4N_O1_A6.fcs  |SUM-159_NLS_4N_O1_A6_seedT2  |2024-12-23        |           35|
|Sample_SUM159_4N_O2_A6.fcs  |SUM-159_NLS_4N_O2_A6_seedT2  |2024-12-23        |           35|
|Sample_SUM159_2N_C_A12.fcs  |SUM-159_NLS_2N_C_A12_seedT1  |2024-12-30        |           42|
|Sample_SUM159_4N_C_A12.fcs  |SUM-159_NLS_4N_C_A12_seedT1  |2024-12-30        |           42|
|Sample_SUM159_4N_O1_A12.fcs |SUM-159_NLS_4N_O1_A12_seedT2 |2025-06-11        |          205|
|Sample_SUM159_4N_O2_A12.fcs |SUM-159_NLS_4N_O2_A12_seedT2 |2025-06-11        |          205|
|Sample_SUM159_2N_O1_A12.fcs |SUM-159_NLS_2N_O1_A12_seedT3 |2025-06-13        |          207|
|Sample_SUM159_2N_O2_A12.fcs |SUM-159_NLS_2N_O2_A12_seedT3 |2025-06-13        |          207|
|Sample_SUM159_2N_O1_A18.fcs |SUM-159_NLS_2N_O1_A18_seedT2 |2025-08-20        |          275|
|Sample_SUM159_2N_O2_A18.fcs |SUM-159_NLS_2N_O2_A18_seedT2 |2025-08-20        |          275|
|Sample_SUM159_4N_O2_A19.fcs |SUM-159_NLS_4N_O2_A19_seedT2 |2025-09-10        |          296|
|Sample_SUM159_4N_O1_A19.fcs |SUM-159_NLS_4N_O1_A19_seedT3 |2025-09-12        |          298|
|Sample_SUM159_2N_O1_A23.fcs |SUM-159_NLS_2N_O1_A23_seedT2 |2025-09-29        |          315|
|Sample_SUM159_2N_O2_A23.fcs |SUM-159_NLS_2N_O2_A23_seedT2 |2025-09-29        |          315|
|Sample_SUM159_4N_O1_A22.fcs |SUM-159_NLS_4N_O1_A22_seedT2 |2025-10-06        |          322|
|Sample_SUM159_4N_O2_A22.fcs |SUM-159_NLS_4N_O2_A22_seedT2 |2025-10-06        |          322|

# Helper Functions



# Build or Load Dataset Results



# Cross-Dataset Gate Counts


|dataset_id                       |population                           | n_samples| min_count| median_count| max_count|
|:--------------------------------|:------------------------------------|---------:|---------:|------------:|---------:|
|anoxia-flowcytometry             |/Margin_Clean                        |        24|     12967|      33475.0|     51306|
|anoxia-flowcytometry             |/Margin_Clean/PeacoQC_Clean          |        24|     12967|      33475.0|     51306|
|anoxia-flowcytometry             |/Margin_Clean/PeacoQC_Clean/Singlets |        24|     12085|      31513.5|     50537|
|hypoxia-sum159                   |/Margin_Clean                        |        20|     24705|      30138.5|     60660|
|hypoxia-sum159                   |/Margin_Clean/PeacoQC_Clean          |        20|     21750|      28324.0|     46750|
|hypoxia-sum159                   |/Margin_Clean/PeacoQC_Clean/Singlets |        20|     20946|      26961.5|     42275|
|polyploidization-ethanolfixation |/Margin_Clean                        |        21|      3469|      13156.0|     27687|
|polyploidization-ethanolfixation |/Margin_Clean/PeacoQC_Clean          |        21|      2250|      12876.0|     27687|
|polyploidization-ethanolfixation |/Margin_Clean/PeacoQC_Clean/Singlets |        21|      2223|      12686.0|     26919|

# Cross-Dataset Metadata



# Figure Galleries


## FSC vs SSC

<div class="plot-gallery" id="gallery-fsc-ssc"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-fsc-ssc')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-fsc-ssc')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/fsc-ssc-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/fsc-ssc-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/fsc-ssc-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## FSC-A vs Time

<div class="plot-gallery" id="gallery-fsc-time"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-fsc-time')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-fsc-time')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/fsc-time-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/fsc-time-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/fsc-time-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## DNA Scatter

<div class="plot-gallery" id="gallery-dna-scatter"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-dna-scatter')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-dna-scatter')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/dna-scatter-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/dna-scatter-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/dna-scatter-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## DNA Density

<div class="plot-gallery" id="gallery-dna-density"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-dna-density')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-dna-density')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/dna-density-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/dna-density-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/dna-density-plots-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Other Fluorescence Area Histograms

<div class="plot-gallery" id="gallery-fluorescence-area-histograms"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-fluorescence-area-histograms')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-fluorescence-area-histograms')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/fluorescence-area-histograms-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/fluorescence-area-histograms-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/fluorescence-area-histograms-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Other Fluorescence vs DNA-A

<div class="plot-gallery" id="gallery-fluorescence-vs-dna"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-fluorescence-vs-dna')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-fluorescence-vs-dna')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/fluorescence-vs-dna-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/fluorescence-vs-dna-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/fluorescence-vs-dna-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## DNA Threshold Fraction Summary

<div class="plot-gallery" id="gallery-dna-fraction-summary"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-dna-fraction-summary')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-dna-fraction-summary')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/dna-fraction-summary-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/dna-fraction-summary-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/dna-fraction-summary-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Raw DNA-A Peak Locations

<div class="plot-gallery" id="gallery-raw-dna-peak-locations"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-raw-dna-peak-locations')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-raw-dna-peak-locations')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/raw-dna-peak-locations-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/raw-dna-peak-locations-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/raw-dna-peak-locations-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Raw DNA-A Peak Location Covariation

<div class="plot-gallery" id="gallery-raw-dna-peak-location-covariation"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-raw-dna-peak-location-covariation')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-raw-dna-peak-location-covariation')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/raw-dna-peak-location-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/raw-dna-peak-location-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/raw-dna-peak-location-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Raw DNA-A Peak Width Covariation

<div class="plot-gallery" id="gallery-raw-dna-peak-width-covariation"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-raw-dna-peak-width-covariation')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-raw-dna-peak-width-covariation')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/raw-dna-peak-width-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/raw-dna-peak-width-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/raw-dna-peak-width-covariation-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>

## Raw DNA-A Tail Fraction vs Above-Threshold Peak

<div class="plot-gallery" id="gallery-raw-dna-peak-vs-fraction"><div class="plot-gallery-controls"><button type="button" onclick="plotGalleryPrev('gallery-raw-dna-peak-vs-fraction')">Prev</button><button type="button" onclick="plotGalleryNext('gallery-raw-dna-peak-vs-fraction')">Next</button><span class="plot-gallery-status"></span></div><div class="plot-slide" data-index="0" style="display:block;"><img src="../figure/anoxia-flowcytometry/raw-dna-peak-vs-fraction-above-threshold-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Anoxia_FlowCytometry (anoxia-flowcytometry)</div></div><div class="plot-slide" data-index="1" style="display:none;"><img src="../figure/hypoxia-sum159/raw-dna-peak-vs-fraction-above-threshold-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Hypoxia_SUM159 (hypoxia-sum159)</div></div><div class="plot-slide" data-index="2" style="display:none;"><img src="../figure/polyploidization-ethanolfixation/raw-dna-peak-vs-fraction-above-threshold-1.png" class="zoomable-plot gallery-img" style="width:100%;" /><div class="plot-caption">Polyploidization_EthanolFixation (polyploidization-ethanolfixation)</div></div></div>


<script>
window.plotGalleryState = window.plotGalleryState || {};

function initPlotGallery(id) {
  var gallery = document.getElementById(id);
  if (!gallery) return;

  var slides = gallery.querySelectorAll(".plot-slide");
  var statusEl = gallery.querySelector(".plot-gallery-status");
  if (!slides.length) return;

  window.plotGalleryState[id] = {
    i: 0,
    n: slides.length,
    panzoomInit: {}
  };

  updatePlotGallery(id);
}

function updatePlotGallery(id) {
  var gallery = document.getElementById(id);
  var state = window.plotGalleryState[id];
  if (!gallery || !state) return;

  var slides = gallery.querySelectorAll(".plot-slide");
  var statusEl = gallery.querySelector(".plot-gallery-status");

  Array.prototype.forEach.call(slides, function(slide, idx) {
    slide.style.display = idx === state.i ? "block" : "none";
  });

  if (statusEl) {
    statusEl.textContent = (state.i + 1) + " / " + state.n;
  }

  var activeSlide = slides[state.i];
  if (!activeSlide) return;

  var img = activeSlide.querySelector("img.zoomable-plot");
  if (!img || typeof panzoom === "undefined") return;

  if (!state.panzoomInit[state.i]) {
    panzoom(img, { maxZoom: 10, minZoom: 1 });
    state.panzoomInit[state.i] = true;
  }
}

function plotGalleryNext(id) {
  var state = window.plotGalleryState[id];
  if (!state) return;
  state.i = (state.i + 1) % state.n;
  updatePlotGallery(id);
}

function plotGalleryPrev(id) {
  var state = window.plotGalleryState[id];
  if (!state) return;
  state.i = (state.i - 1 + state.n) % state.n;
  updatePlotGallery(id);
}

document.addEventListener("DOMContentLoaded", function() {
  document.querySelectorAll(".plot-gallery").forEach(function(el) {
    initPlotGallery(el.id);
  });
});
</script>
