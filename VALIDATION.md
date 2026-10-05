# Verification of COQ Python 1.0.0

Verification date: **5 October 2026**.

## Completed checks

| Check | Result |
|---|---|
| Numerical and file-format suite | **26 passed**, no failures |
| GUI workflow smoke test | **18 passed**, no errors |
| Python compilation | Application, numerical engine, and tests compiled successfully |
| Gaussian edge/slit MTF | MTF50 and MTF10 within the 2.5% analytical test tolerance for both horizontal and vertical orientations |
| White-noise NPS | Mean NNPS and integrated NPS agree with expected noise/variance within the 2.5% statistical test tolerance; tested with no window and Hann window |
| Response fitting | All four models recover synthetic dose values and achieve R² above 0.999999 |
| Original-JAR FFT comparison | Relative L2 difference **2.0529 × 10⁻⁷** against `2*abs(FFT)^2`; maximum absolute difference **1.4305 × 10⁻⁶** for the tested 32 × 32 float fragment |
| Original-JAR AGD g/c comparison | Original-index mode agrees at all eight PMMA selections at HVL 0.36 mm Al; maximum absolute coefficient difference **5.5511 × 10⁻¹⁷** |
| Original-JAR mean/sample-SD helpers | Identical values for the tested integer sequence: mean **30.5**, sample SD **17.46424919657298** |
| GUI layout | Rendered and visually inspected at 1480 × 940; tabs, ROI selection, image display and result plot checked |

## Numerical and file-format coverage

The suite exercises:

- Response fitting and forward/inverse consistency for linear, log(K), log(1+cK), and square-root models.
- Edge and slit MTF with a known Gaussian line-spread function in both orientations; Nyquist limits; rejection of an unslanted feature.
- White-noise NPS normalization, variance integration, window-power correction, quadratic detrending and invalid fragment/spacing handling.
- DQE from known curves, common frequency support, positive kerma and positive NNPS requirements.
- Bad-pixel coordinates, partial tiles and full tested-pixel counts.
- Uniformity's global formula and maximum absolute local deviation.
- Five-ROI dark mean/sample-SD statistics.
- Pooled-background CNR and overlap rejection.
- Additive/multiplicative lag equations and image order.
- AGD polynomial mapping, manual factors, guarded HVL range and inverse-square geometry.
- Tube-table interpolation and unavailable-zero rejection.
- RAW big-endian values, header offsets, interframe gaps and truncated-file rejection.
- DICOM preference for detector spacing, preserved raw MONOCHROME1 signal, optional rescale and dGy-to-mGy conversion.
- TIFF/NumPy input and CSV/JSON export/import.

## GUI workflow coverage

The smoke test runs actual Tkinter event processing and background calculation callbacks using synthetic images. It verifies:

1. Edge MTF.
2. Slit MTF.
3. NNPS.
4. DQE.
5. Importing exported MTF/NNPS and recalculating DQE.
6. JSON export.
7. CSV and PNG plot export.
8. Uniformity.
9. Bad-pixel analysis.
10. Dark analysis.
11. Additive lag.
12. Multiplicative lag.
13. CNR.
14. Batch CNR.
15. Response fitting and automatic linearization.
16. AGD with entered kerma.
17. Tube-output lookup and AGD.
18. Batch AGD from the supplied template.

It also verifies that clearing the image list preserves a previously computed MTF for use with a later flat-field acquisition. Dialogs are replaced with temporary file paths so export/import checks do not require manual interaction.

## Original-JAR comparisons

The inspected user-supplied JAR was run locally with ImageJ classes for selected numerical methods. Its SHA-256 is recorded in `coq_legacy_factors.json`. The direct comparisons establish the FFT scaling convention, original g/c index mapping and sample-statistics convention. They do **not** establish agreement of the complete GUI, MTF tracing, uniformity grid, or clinical DQE workflow.

AGD polynomial comparison at HVL 0.36 mm Al:

| PMMA label (mm) | Original Java g | Original Java c |
|---:|---:|---:|
| 20 | 0.420280 | 0.89235200 |
| 30 | 0.420280 | 0.94665792 |
| 40 | 0.297688 | 0.94665792 |
| 45 | 0.210328 | 1.09120880 |
| 50 | 0.182056 | 1.07960960 |
| 60 | 0.157144 | 1.15991040 |
| 70 | 0.123424 | 1.24662208 |
| 80 | 0.099672 | 1.25667040 |

The default named-thickness mode deliberately selects a different polynomial at some thicknesses. For example, at PMMA 45 mm and HVL 0.36 mm Al, default g=0.182056 and c=1.07960960, whereas original-index mode gives g=0.210328 and c=1.09120880. Neither comparison validates these legacy approximations as a substitute for verified contemporary conversion tables.

## Environment

Tests were run under Linux with Python 3.12 and a virtual X display for Tkinter. Runtime package versions used were:

| Package | Version |
|---|---|
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| Matplotlib | 3.10.8 |
| pydicom | 3.0.2 |
| tifffile | 2026.9.20 |
| Pillow | 12.3.0 |
| pytest, tests only | 9.1.1 |
| JPype1, Java comparisons only | 1.7.1 |

The application does not require pytest, JPype, Java, or a virtual X server for normal desktop use. The Windows batch launchers were inspected but not executed on Windows.

## Remaining validation

No clinical detector-image dataset was supplied. The release therefore has not been validated against actual system response, manufacturer's preprocessing, beam-specific q, measured AGD coefficients, or complete ImageJ COQ output on clinical QC acquisitions. See README.md for deliberate implementation differences, exposure requirements, imported-NNPS unit conversion and supported metadata limits.

For deployment, compare source-image ROI values first, then response inversion, MTF profiles/curves, NPS normalization and extraction, and finally DQE/AGD. Keep acquisition settings and pixel spacing identical and document the tolerances appropriate to the intended protocol. The supplied `Data_AGD.csv` remains an example calibration table, not a measured characterization of your equipment.
