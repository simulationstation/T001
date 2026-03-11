# Strict Action Plan

## Objective
Replace the current source-list-first pilot with a stricter difference-imaging workflow that can survive the specific objections raised in review.

## Required Changes
1. Register every image pair empirically from stars before subtraction.
2. Detect on the residual image itself, not only on a pre-image source list.
3. Carry explicit mask, edge, chip-gap, and non-detection provenance through every measurement.
4. Run forced photometry on registered difference images for every surviving candidate across all available same-filter epochs.
5. Build a blinded benchmark from targeted supernova observations already present in the archive and evaluate recall outside the detector logic.

## Acceptance Criteria
1. The upgraded detector must write pair-level registration diagnostics: shift, residual, and number of matched anchors.
2. The upgraded detector must output residual detections for both fading and brightening events.
3. Surviving fade candidates must have registered difference light curves, including low-S/N and non-detection epochs.
4. The blinded benchmark truth table must live in a separate output tree and must not be consumed by the detector stage.
5. The benchmark report must quantify at least:
   - truth groups
   - benchmark pairs
   - recovered truth groups
   - recovery fraction
   - localization error

## Execution
1. Implement a new `difference_upgrade` pipeline stage with empirical registration and residual-first detection.
2. Run the upgraded fade search on the existing archive pair queue.
3. Run difference follow-up on the upgraded surviving fade candidates.
4. Build the hidden supernova benchmark and run the same detector in blind `both-sign` mode on the benchmark pair set.
5. Write a final report summarizing:
   - upgraded fade candidates
   - whether `NGC 5861` survives
   - whether `MESSIER 101` survives
   - benchmark recovery performance

