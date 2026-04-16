# Provisional Current DNA Content Flow Cytometry Prep Protocol

## Status

This document is a working reconstruction of the current sample-preparation
workflow used for the fixed-cell DNA-content experiments in this repository. It
is not yet a confirmed wet-lab SOP.

The sequence below is consistent with the project's current modeling
assumptions, result interpretation, and archived checkpoint notes:

- target cells are formaldehyde-fixed
- commercially pre-fixed CEN are added after target-cell fixation
- permeabilization and DNA staining are performed on the mixed target/CEN
  suspension
- DNA content is acquired from the violet-excited DNA channel

Operational parameters that are not confirmed from the repo alone are marked as
`Needs wet-lab confirmation`.

## Objective

Prepare cultured cancer cell suspensions for downstream DNA-content flow
cytometry analysis using a post-fix CEN spike-in and FxCycle Violet staining.

## Reagents and Materials

- Phosphate-buffered saline (PBS), cold
- Formaldehyde solution
- Permeabilization buffer
  `Needs wet-lab confirmation`: exact buffer identity and concentration
- Commercial pre-fixed chicken erythrocyte nuclei (CEN) standard
- FxCycle Violet stain
- Trypan Blue solution
- Hemocytometer or automated cell counter
- Flow cytometry tubes

## Phase I: Sample Harvest and Wash

1. Harvest target cells to obtain a uniform single-cell suspension.
2. Stain a small aliquot with Trypan Blue and count cells to estimate
   concentration and viability.
3. Adjust the suspension with cold PBS to a target density of about
   `1 x 10^6 cells/mL`.
   `Needs wet-lab confirmation`: target density and acceptable operating range.
4. Centrifuge to pellet the cells and aspirate the supernatant.
5. Resuspend the pellet in cold PBS to wash away residual culture media and
   serum proteins.
6. Centrifuge again and aspirate the PBS wash.

## Phase II: Fixation

1. Resuspend the washed cell pellet in formaldehyde solution.
2. Incubate for about `15 minutes at room temperature in the dark`.
   `Needs wet-lab confirmation`: exact formaldehyde concentration, fixation
   duration, temperature, and whether dark incubation is required.
3. Centrifuge the fixed cell suspension.
4. Aspirate the fixative and wash the pellet with PBS to halt the reaction.

## Phase III: Internal Reference Spiking

1. Prepare a standardized aliquot of the commercially pre-fixed CEN standard.
2. Spike the CEN aliquot into the suspension of formaldehyde-fixed target
   cells.
3. Mix gently to ensure even distribution.

Note:
- This post-fix spike-in arrangement is consistent with the current repo
  assumptions and with the observed batch-structure discussed in the checkpoint
  analysis.
- It also means the reference and target cells do not share the same fixation
  history, which is a known interpretive caveat in this project.

## Phase IV: Permeabilization and Staining

1. Add the designated permeabilization buffer to the mixed target/CEN
   suspension and incubate according to the buffer requirements.
   `Needs wet-lab confirmation`: exact buffer composition, concentration, and
   incubation duration.
2. Add FxCycle Violet stain.
   Provisional working assumption: `1 uL per 1 x 10^6 cells`.
   `Needs wet-lab confirmation`: whether the lab is using the manufacturer
   baseline ratio or an internally optimized titration.
3. Incubate the samples in the dark at room temperature for about `30 minutes`
   to allow dye binding.
   `Needs wet-lab confirmation`: exact staining time and temperature.

## Phase V: Data Acquisition

1. Acquire the samples on the flow cytometer using the `405 nm` violet laser.
2. Record the violet DNA area channel used by the analysis pipeline.

Note:
- In the current repo workflow, the default DNA channels are
  `450-Violet C-A` and `450-Violet C-H`.
- FSC/SSC and time channels are also used for QC and gating in
  [workflow/data_prep.Rmd](../workflow/data_prep.Rmd).

## Project-Specific Caveats

- The current analysis treats separately fixed controls as potentially
  problematic. In both the anoxia and hypoxia-style datasets, the `2N` and `4N`
  controls appear systematically different from the main sample set.
- The leading hypothesis recorded in the repo is that fixation history matters:
  the controls were formaldehyde-fixed separately, whereas permeabilization and
  staining were performed contemporaneously with the rest of the samples.
- The current preferred model, `observed_cen_scaled`, was motivated in part by
  covariance between the observed CEN peak and the tumor G0/G1 peak. That makes
  protocol consistency around fixation, CEN handling, permeabilization, and
  staining especially important.

## Open Items To Confirm

- Formaldehyde concentration
- Exact fixation time and temperature
- Whether fixation is performed in the dark
- Permeabilization buffer identity, concentration, and incubation time
- FxCycle Violet volume per test or per cell count
- Post-fix wash details
- Centrifugation speeds and durations
- Whether the target cell input density is actually standardized at
  `1 x 10^6 cells/mL`
