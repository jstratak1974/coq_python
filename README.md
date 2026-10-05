# COQ Python: detector characterization and quality assurance

A standalone Python application with a Tkinter graphical interface, based on the supplied **COQ ImageJ plugin**, its user manual, and inspection of its compiled calculation classes.

Python was selected because the application can run without ImageJ, Java, or a MATLAB licence, and provides DICOM support, interactive ROIs, scientific plots, and portable numerical code. The GUI and calculation engine are separate so that the same functions can also be used in notebooks or automated analyses.

**This is an independent reimplementation, not a byte-for-byte translation or an officially validated replacement for COQ.** Several calculation details differ deliberately and are documented below. It has been tested with synthetic data and selected original-JAR calculations. Clinical detector images and a complete end-to-end comparison with ImageJ COQ were not supplied. Establish agreement with your own QC datasets before routine use; this release does not claim IEC compliance or medical-device certification.

## Files

| File | Purpose |
|---|---|
| `coq_python.py` | Tkinter GUI and application entry point |
| `coq_core.py` | Numerical functions and image readers; usable without launching the GUI |
| `requirements.txt` | Python dependencies |
| `Data_AGD.csv` | Your supplied two-section HVL and tube-output example table |
| `coq_legacy_factors.json` | g/c polynomial coefficients and s factors recovered from the JAR |
| `agd_batch_template.csv` | Example batch input for mammography AGD |
| `install_windows.bat` | Install dependencies with the Windows Python launcher |
| `run_windows.bat` | Launch the application on Windows |
| `tests/test_core.py` | Analytical, regression, and image-format tests |
| `tests/gui_smoke_test.py` | GUI workflow checks on a system with a display |
| `VALIDATION.md` | Actual verification results and remaining limitations |

Keep the Python files, JSON coefficients, and CSV table in the same directory.

## Installation and launch

### Windows

1. Install Python 3.10 or later from [python.org](https://www.python.org/downloads/). Include the Tcl/Tk component and Python launcher.
2. Extract the complete ZIP into a folder.
3. Double-click `install_windows.bat`, or run:

   ```powershell
   py -m pip install -r requirements.txt
   ```

4. Double-click `run_windows.bat`, or run:

   ```powershell
   py coq_python.py
   ```

If your installation exposes `python` instead of `py`, use that command in the same folder. The application runs locally and does not upload image data.

### Linux

Ensure that Tkinter and virtual-environment support are installed for your Python. For Ubuntu:

```bash
sudo apt install python3-tk python3-venv
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python coq_python.py
```

### macOS

Use a Python distribution that includes Tk support, create a virtual environment, install the requirements, and launch `python coq_python.py`.

### Open files at startup

```bash
python coq_python.py flat_01.dcm flat_02.dcm
```

The script itself requires a desktop/display and a screen of at least 1280 × 800 pixels. Importing `coq_core` for numerical work does not require a display. Optional tests require `pytest`; Java comparisons do not form part of the runtime dependencies.

## Implemented functions

| COQ workflow | Python implementation |
|---|---|
| Response curve | Linear, logarithmic with/without offset inside the logarithm, square-root fits; manual coefficients; calibration JSON import/export; inverse response linearization |
| Edge MTF | Vertical/horizontal or automatic orientation; line tracing; oversampled ESF; LSF; presampled MTF; MTF50/MTF10 |
| Slit MTF | Slanted slit tracing; oversampled LSF; optional finite-slit-width correction |
| NNPS | Multiple flat fields; overlapping fragments; horizontal, vertical, directional-average and radial spectra; 2D NPS/NNPS; mean/plane/quadratic detrending |
| DQE | Stored or imported MTF/NNPS; COQ beam constants or custom q; explicit kerma; common-frequency interpolation |
| Bad pixels | Local 10 mm tiles by default; 20% deviation threshold; coordinates and tested-pixel counts |
| Uniformity | Signal/SNR maps; local and global nonuniformity; legacy comparison limits |
| Dark analysis | Five ROIs: centre and four corners; mean and sample standard deviation |
| Lag | Additive and multiplicative lag with documented image order |
| CNR | Two-ROI or five-ROI method; pooled background noise; active-image or batch calculation |
| Mammography AGD | Entered incident kerma or calibrated tube output; g/c/s factors; original-index comparison mode; batch CSV input |
| Image tools | DICOM/TIFF/NumPy/RAW input; multiframe images; selectable ROIs; display window; linearized-image export |
| Results | CSV and JSON export with settings; PNG/PDF/SVG plots |

The modality selector follows the documented COQ workflow: CR excludes dark analysis and lag, CR uniformity reports global values, and angiography offers additive lag. AGD/CNR tools remain accessible for expert use; choose mammography when interpreting these measurements.

## Image handling and geometry

- Supported inputs: grayscale DICOM (`.dcm`, `.ima`, or files selected through All files), TIFF, `.npy`, grayscale PNG/BMP/JPEG, and RAW.
- TIFF, NumPy and DICOM stacks are loaded as separate frames. Image order is preserved and can be changed with **Move up / Move down**.
- RGB images are rejected instead of being silently converted. Rendered PNG/JPEG images may have lost detector signal information; use original QC images for MTF/NPS/DQE.
- No automatic resizing or resampling is performed. NNPS and lag reject mismatched image geometry.
- DICOM **Imager Pixel Spacing** is preferred, then **Nominal Scanned Pixel Spacing**, then **Pixel Spacing**. The latter may describe patient-plane calibration rather than detector sampling; verify the displayed values.
- `dy` is row spacing and `dx` is column spacing, both in **mm**. Missing spacing is left blank. Enter measured values for non-DICOM images.
- MTF and physical QC dimensions support anisotropic spacing. This release restricts the 1D NNPS/DQE extraction to square pixels.
- Stored pixel values are used by default. **Apply DICOM rescale slope/intercept when opening** explicitly enables that transformation. No VOI LUT, presentation LUT, gamma adjustment, or display inversion is applied to numerical data.
- MONOCHROME1 affects the displayed grayscale only. Windowing and the image toolbar do not change calculation values.
- Compressed DICOM requires a suitable optional pydicom decoder. If decoding fails, export uncompressed source images or install the decoder corresponding to the transfer syntax. See [pydicom's pixel-data documentation](https://pydicom.github.io/pydicom/stable/guides/user/image_data_handlers.html).

### RAW input

Choose **Import RAW** and supply rows, columns, number of frames, header offset, gap between frames, numeric dtype, and byte order. `<` means little endian; `>` means big endian. Spacing must be entered afterward. Incomplete frames are rejected.

## ROIs

Drag with the left mouse button on the image, or enter `X`, `Y`, `W`, `H` and click **Set**. The numerical ROI is always `(x, y, width, height)`, with **zero-based pixel coordinates** and exclusive upper array boundaries. Click **Full** to analyse the complete active image. **Magic ROI** places a user-sized ROI centrally.

The current selection is the Analysis ROI used for response fitting, MTF, NNPS and QC. **Store** saves it under a chosen name; **Recall** restores the saved selection. CNR and lag require the following named ROIs:

| ROI name | CNR | Lag |
|---|---|---|
| Signal / Inside | Aluminium/detail signal | Inside the test-phantom region |
| Background 1 / Outside | Background | Outside the test-phantom region |
| Background 2–4 | Additional backgrounds for the five-ROI method | Unused |

Named ROIs are overlaid on the image. Out-of-bounds ROIs cause an error rather than being cropped silently. Place CNR signal and background ROIs in distinct, nonoverlapping regions.

## Response fitting and linearization

1. Load at least four flood images acquired at distinct measured detector air-kerma levels.
2. Select the same central ROI in every image, ideally approximately 20 × 20 mm as in the manual.
3. Enter one air-kerma value per loaded image, in microGy, separated by commas.
4. Select the model and click **Fit response from all images**.
5. Review coefficients, fit curve and R². A successful fit enables **Use response linearization**.

The response models are:

| Model | Signal S as a function of kerma K | Inverse used for linearization |
|---|---|---|
| Linear | `S = a + b*K` | `K = (S-a)/b` |
| Log(1+cK) | `S = a + b*ln(1+c*K)` | `K = expm1((S-a)/b)/c` |
| Log(K) | `S = a + b*ln(K)` | `K = exp((S-a)/b)` |
| Square root | `S = a + b*sqrt(K)` | `K = ((S-a)/b)^2`, on the valid response branch |
| Identity | `S = K` | Values unchanged |

Natural logarithms are used. The logarithmic offset model requires positive `c`; increasing or decreasing response slopes are allowed. Non-finite output, zero slopes and ambiguous square-root inversion are rejected. R² describes fit agreement, not measurement uncertainty. Stored calibration ranges are recorded but do not prove validity outside the fitted exposure range.

Manual coefficients can be applied and saved/loaded as JSON. Changing calibration clears stored MTF/NNPS measurements so that DQE cannot accidentally combine obsolete calibration states. **Save linearized** exports numerical float TIFF or NumPy values using the current calibration, irrespective of the display checkbox; it does not create a DICOM object or alter the source image.

## MTF: edge or slit

1. Open an unprocessed edge/slit QC image and verify physical detector spacing.
2. Draw a rectangle around a single straight feature with a slight slant. A 1–5° angle is preferred; the application rejects angles outside 0.5–15° and weak/poorly traced features.
3. Select **edge** or **slit**. Orientation can be automatic, vertical, or horizontal.
4. Set oversampling (default 10) and profile length (default 10 mm).
5. Click **Calculate MTF**.

The line is fitted from subpixel row traces. Pixel samples are binned by **physical normal distance** to the fitted line; empty bins are interpolated only if at most 10% of bins are missing. The profile is limited to the length present in every selected row. ESF differentiation uses a first difference. A Hamming window is applied to the LSF, and the Fourier amplitude is normalized to MTF(0)=1. The optional edge correction divides by the first-difference transfer response `sinc(f*delta)`.

For slit data, the tails define a constant background and the oversampled profile is used as the LSF. Enter the measured slit width in mm to divide by its `sinc(f*width)` response; zero disables this correction. Corrections near the first zero are rejected.

Frequency is in **lp/mm**, numerically equivalent to cycles/mm. MTF is reported through Nyquist `1/(2*dx)` for a vertical feature or `1/(2*dy)` for a horizontal feature. MTF50/MTF10 use the first downward crossing. If a crossing is absent in the measured range, the result is `null`, not an extrapolated value.

Inspect **Curve**, **ESF / slit profile**, and **LSF** plot views. The output also reports angle, trace residual, used rows, and empty-bin fraction.

## NNPS and DQE

### NNPS

Use one or more original flat-field images at the **same exposure and detector configuration**. Select a sufficiently large, uniform Analysis ROI, verify spacing and linearization, then click **Calculate NNPS**.

Defaults: 256 × 256 pixel fragments, 50% overlap, subtraction of the mean, no FFT window, one central line omitted, and seven lines averaged on each side of each axis. Plane/quadratic detrending and a Hann window are optional, and change the spectrum relative to the legacy plugin.

For a fragment with `N = rows*columns` pixels, pixel area `A = dy*dx`, mean signal `mu`, residual image `r`, window `w`, and window power `U = mean(w^2)`:

```text
NPS2D  = A/(N*U) * |FFT2(r*w)|^2
NNPS2D = NPS2D / mu^2
```

NPS has units **signal²·mm²**, and NNPS has units **mm²**. Each fragment is normalized by its own mean before averaging. The optional COQ-style rejection excludes fragments whose range is greater than or equal to the absolute mean. Accepted and rejected counts are recorded; a zero mean or rejection of all fragments causes an error.

Horizontal/vertical spectra average bands around their axes and both frequency signs. Radial averaging uses annular bins and omits the configured axis lines. `omit=1` excludes the central line; `lines=7` uses offsets 1–7 on either side. DC and Nyquist are excluded from the extracted 1D curve. The complete 2D arrays are included in JSON exports; CSV contains the 1D curves.

### DQE workflow with separate images

1. Calculate MTF from the edge/slit acquisition.
2. **Clear** the image list and open the flat-field acquisition. Clearing/opening images preserves stored MTF/NNPS curves.
3. Apply the same calibration state, verify spacing and exposure, select the flat-field Analysis ROI, and calculate NNPS.
4. Select the NNPS direction, appropriate beam or **Custom**, detector air kerma in microGy, and `q/kerma` in **1/mm²/microGy**.
5. Click **Calculate DQE from stored MTF + NNPS**.

```text
q = air_kerma_microGy * q_per_microGy
DQE(f) = MTF(f)^2 / (q * NNPS(f))
```

The RQA and mammography beam constants are copied from the supplied COQ JAR and are labelled legacy values. The mammography presets are for the explicitly described beam conditions, including added 2 mm Al; they are not universal q values for every clinical spectrum. Enter a verified custom value for other beams.

Linear interpolation is the default. **COQ quartic** reproduces the original choice of a fourth-degree polynomial fit for MTF interpolation. Both options are limited to overlapping frequency support; no curve is extrapolated. DQE values above 1 are preserved and highlighted so that inconsistent inputs remain visible.

### Importing curves

Simple CSV/TSV curves should contain frequency in lp/mm and MTF or NNPS in the second column. Python NNPS exports contain separate horizontal, vertical, average, and radial columns. Original COQ four-column tables can also be read, but their numerical **units must first be reconciled** as described below. Imported files do not provide trustworthy acquisition/calibration equivalence; verify these conditions yourself.

## Uniformity, defective pixels and dark noise

### Bad pixels

The defaults match the compiled plugin: 10 mm × 10 mm tiles, a 15 mm outer margin, and pixels satisfying:

```text
abs(pixel - tile_mean) > 0.20 * abs(tile_mean)
```

Partial edge tiles are included. Coordinates refer to the original image, including the Analysis ROI offset. Zero-mean tiles are skipped and counted explicitly. This is a local deviation test; it does not reconstruct a manufacturer's detector defect map.

### Uniformity

The defaults are 30 mm × 30 mm ROIs, overlapped by 50%. Each ROI yields its mean, sample SD and SNR. For a signal/SNR map `z`:

```text
Global NU = 200 * (max(z)-min(z)) / (max(z)+min(z))
Local NU  = max[100 * abs(z_center - mean(eight_neighbours))
                    / abs(mean(eight_neighbours))]
```

Only interior grid ROIs contribute to local nonuniformity. With fewer than three ROI rows/columns, local results are unavailable. Legacy comparison limits are LSNU 4%, GSNU 8%, LSNRNU 8%, and GSNRNU 20%. These labels are retained for comparison, not presented as current universal acceptance criteria. The evaluated fraction of the selected area and maps are reported.

### Dark analysis

The compiled JAR uses **five 55 mm square ROIs with a 30 mm margin**, despite the manual's dark-analysis paragraph repeating the defective-pixel instructions. This release follows that five-ROI operation, with configurable size/margin, reporting mean and sample SD in the centre and four corners. Smaller demonstration images require reduced ROI size/margin. No unsupported dark-noise pass/fail limit is invented.

## Lag and CNR

### Lag

Images must have identical geometry. Select the ROIs on the test image and preserve acquisition timing and processing state.

| Mode | Required image-list order |
|---|---|
| Additive | Test phantom image, then nonirradiated image acquired afterward |
| Multiplicative | Test phantom image, uniform image acquired before the test, uniform image acquired after the test |

With `I` denoting inside-phantom mean and `O` outside-phantom mean:

```text
Additive lag (%)       = 100 * abs(I_dark - O_dark) / O_test
Multiplicative lag (%) = 100 * abs((I_before-O_before) - (I_after-O_after))
                              / ((O_before+O_after)/2)
```

The legacy COQ comparison limit is `<0.5%`. The test image defines the phantom footprint in multiplicative lag; its means do not enter the final formula.

### CNR

Store one detail ROI and either one or four background ROIs. The background pixels are pooled, including any variation between the background ROIs:

```text
CNR = abs(mean(background) - mean(detail)) / sample_SD(background)
```

Click **Calculate CNR on active image**, or **Calculate CNR on all images** after positioning corresponding ROIs. When phantom positions differ between acquisitions, analyse each image with appropriately repositioned ROIs rather than using one batch geometry.

## AGD / mean glandular dose

AGD is calculated as:

```text
AGD_mGy = incident_air_kerma_mGy * g * c * s
```

Use **incident air kerma without backscatter** at the phantom/breast entrance. A generic entrance-dose/ESAK value cannot automatically be substituted if it includes backscatter. The GUI deliberately displays DICOM dose fields for verification without silently treating them as incident kerma.

### Entered kerma

Select **Entered incident kerma**, enter the measured kerma in mGy, HVL and PMMA thickness, then choose a factor mode. A protocol-specific AGD acceptance limit can optionally be entered; this release does not apply the original incomplete hardcoded AGD limit array.

### Tube output

Select **Tube output** and provide output in microGy/mAs, its reference distance, mAs, and focus-to-entrance distance:

```text
K_mGy = mAs * output_microGy_per_mAs / 1000
                * (reference_distance_mm / entrance_distance_mm)^2
```

The entrance distance must be explicitly checked. `SID - object_thickness` assumes the entrance surface is at that geometry and does not account for an unknown detector/support offset.

**Tube lookup** reads the supplied table at 500 mm; **Load tube table** accepts the same two-section format. kVp interpolation is linear inside the header range (25–34 kV); extrapolation is rejected. Zero entries mean unavailable values and are rejected. The original third column is retained as a descriptor; kVp is read from the numerical column headers. Replace the example table with measured data for the actual equipment before using output-based kerma.

### Factor modes

| Mode | Behaviour |
|---|---|
| Thickness matched COQ | Uses the polynomial whose embedded name matches the entered PMMA thickness |
| Original COQ index | Reproduces the supplied JAR's actual shifted g/c index mapping for comparison |
| Manual g/c/s | Uses explicitly entered, independently verified conversion factors |

Legacy polynomial modes accept PMMA 20, 30, 40, 45, 50, 60, 70 and 80 mm. A conservative application guard restricts their HVL to 0.25–0.65 mm Al; this guard is **not** evidence of a validated polynomial range. Unsupported spectra/thicknesses require manual factors. `s` values are 1.000, 1.017, 1.061, 1.044, 1.042 and 1.050 for Mo/Mo, Mo/Rh, Rh/Rh, Rh/Al, W/Rh and W/Al respectively. This release does not substitute new Dance tables or extend these coefficients to DBT, W/Ag, or patient-specific glandularity.

**DICOM acquisition fields** copies available kVp/mAs and a proposed entrance distance; it does not infer that patient Body Part Thickness is PMMA thickness. DICOM Organ Dose (0040,0316) is in dGy and is displayed in mGy using ×100. Entrance Dose in mGy (0040,8302) and its derivation are displayed for review. Multiframe/DBT dose values can represent aggregate acquisition quantities, so they are not automatically used as per-frame incident kerma.

### Batch CSV

Use `agd_batch_template.csv` with these columns:

```text
label,pmma_mm,hvl_mm_al,anode_filter,incident_kerma_mGy,mas,
output_uGy_mAs,reference_distance_mm,entrance_distance_mm,g,c,s,
legacy_index,acceptance_limit_mGy
```

The actual header is a single line. Supply either incident kerma or output/mAs/distance. Enter all three `g,c,s` for manual factors; leave all three blank to use the legacy coefficients. `legacy_index=true` requests the original index mapping. Invalid rows stop the batch and report their CSV line number.

## Documented differences from the supplied plugin

| Topic | Supplied JAR | Python release |
|---|---|---|
| NNPS units | Accumulates squared Hartley-amplitude values after signal normalization; the inspected path has no pixel-area/fragment-size normalization | Applies explicit two-sided physical NPS normalization |
| MTF implementation | COQ-specific trace, binning and profile routines | Physical normal-distance binning, robust row trace, Hamming LSF window and optional derivative correction |
| DQE interpolation | Fourth-degree polynomial fit | Linear interpolation by default; original quartic choice available |
| Uniformity grid | Legacy centred/margin grid can omit additional border regions | Complete fitting windows within the selected region, with evaluated coverage reported |
| Local nonuniformity | Signed running maximum in the inspected routine can mishandle negative extrema | Maximum absolute deviation |
| AGD index mapping | PMMA labels and g/c coefficient indices are inconsistent | Named-thickness mapping by default; original mapping available explicitly |
| AGD limits | Six-entry array used with an eight-thickness selection | User-entered protocol limit; no automatic historical AGD limits |
| Dark analysis | Five-ROI statistics in code; inconsistent manual paragraph | Five-ROI statistics with configurable dimensions |
| Invalid input | Several legacy paths can yield empty spectra or invalid array access | Explicit geometry, calibration, range and normalization errors |
| Image stack geometry | ImageJ can resize differing images according to preferences | No automatic resizing; geometry mismatch rejected |

For **mean detrending and no FFT window**, the supplied JAR's raw NNPS spectrum uses `2*abs(FFT)^2`. Convert it to the physical NNPS convention used here with:

```text
NNPS_physical_mm2 = NNPS_COQ_raw * (dy*dx) / (2*fragment_size^2)
```

This convention was checked directly by running the JAR FFT and comparing it with NumPy. Fragment selection, frequency-band averaging, precision and radial binning can still cause differences. **Do not import an unconverted legacy NNPS spectrum into a physical DQE calculation.**

The supplied JAR's actual AGD mapping is:

| PMMA label (mm) | g polynomial selected | c polynomial selected |
|---:|---:|---:|
| 20 | 20 | 20 |
| 30 | 20 | 30 |
| 40 | 30 | 30 |
| 45 | 40 | 40 |
| 50 | 45 | 45 |
| 60 | 50 | 50 |
| 70 | 60 | 60 |
| 80 | 70 | 70 |

These names identify embedded coefficient arrays, not a newly validated dosimetry table. The original-index mode is provided to support comparison, not to endorse that mapping.

## Using the calculation engine directly

```python
import coq_core as coq

edge = coq.load_images("edge.dcm")[0]
mtf = coq.mtf(edge.pixels, edge.spacing, roi=(200, 300, 256, 512))

flats = [coq.load_images(p)[0] for p in ["flat1.dcm", "flat2.dcm"]]
nnps = coq.nps([im.pixels for im in flats], flats[0].spacing,
               roi=(100, 100, 1024, 1024), fragment=256)

dqe = coq.dqe(mtf["frequency"], mtf["mtf"],
              nnps["frequency"], nnps["average"],
              kerma_uGy=2.5, q_per_uGy=coq.BEAMS["RQA 5"])
coq.export_result("dqe.csv", "DQE", dqe)
```

Adapt the ROIs to actual image bounds and apply the appropriate response inverse to both datasets before characterization. `coq_core` accepts NumPy arrays and does not infer acquisition validity.

## Verification and known boundaries

Run:

```bash
python -m pip install pytest
python -m pytest -q tests/test_core.py
python tests/gui_smoke_test.py
```

The numerical tests include Gaussian edge/slit MTF in both orientations, white-noise NPS normalization and variance integration, all response inverses, DQE overlap, g/c/s and inverse-square calculations, ROI statistics, lag image order, RAW offset/gap/endian handling, DICOM spacing/rescale/dose units, and export/import.

Synthetic demonstrations are built into the GUI. Their expected behaviour is useful for checking installation, but synthetic agreement does not validate a clinical acquisition protocol.

Not included: CDMAM/contrast-detail scoring, DBT angular dose factors, complete enhanced-DICOM acquisition metadata extraction, patient-specific dose models, automatic manufacturer's defect correction, uncertainty propagation, and an exhaustive reproduction of every ImageJ import codec, preference or macro expression. Calibration supports the four documented response models rather than arbitrary executable ImageJ macros.

## Attribution and references

Original COQ developers: **Bruno Donini, Stefano Rivetti, Nico Lanconelli and Marco Bertolini**.

- Donini B, Rivetti S, Lanconelli N, Bertolini M. *Free software for performing physical analysis of systems for digital radiography and mammography*. Medical Physics. 2014;41:051903. [doi:10.1118/1.4870955](https://doi.org/10.1118/1.4870955).
- [COQ original project](https://www.medphys.it/down_dqe.htm).
- User-supplied `COQ_user_manual.pdf` and `COQ_.jar`; numerical classes inspected locally to resolve ambiguities in the manual.
- [DICOM DX Detector Module: detector spacing](https://dicom.nema.org/Medical/Dicom/current/output/chtml/part03/sect_C.8.11.4.html).
- [DICOM Breast Projection Image Modules: dose quantities and derivation](https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_c.8.31.html).

The package contains new Python code and recovered numerical coefficients; it does not redistribute the original Java implementation or third-party bundled Java libraries. Preserve the original project attribution when sharing derived work. The SHA-256 of the supplied JAR is recorded in `coq_legacy_factors.json` so comparison results can be tied to that particular plugin build.
