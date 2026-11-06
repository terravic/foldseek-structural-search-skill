# Dashboard Design Specification

This document defines the visual, structural, and technical specifications for standalone interactive HTML5 analytical dashboards.

---

## 1. Visual Hierarchy and Typography Rules

| Design Dimension | Specification | Implementation Standard |
|---|---|---|
| **Typography** | Distinct typographic hierarchy across prose, headers, and tabular values. | Pair a sans-serif UI font stack (`DM Sans`, `system-ui`, `sans-serif`), a serif display heading stack (`Instrument Serif`, `Georgia`, `serif`), and a tabular monospaced font stack (`DejaVu Sans Mono`, `ui-monospace`, `monospace`). |
| **Surfaces and Depth** | Flat, solid surface tokens with clear structural boundaries. | Use solid background tokens and `1px solid` borders (`6px` to `8px` radius) for separation. |
| **Container Layout** | Single-level bento-style grid organization. | Partition panes using flat grid cells and thin `1px` dividers rather than nested card-in-card containers. |
| **Status Indicators** | High-signal state reporting. | Highlight warning, error, or critical states explicitly; render nominal items using standard monospaced or muted tokens. |
| **Background Surfaces** | Uncluttered viewport backgrounds for coordinate views and workspaces. | Reserve gridlines strictly for quantitative chart axes and coordinate plots. |
| **Control Surfaces** | Compact text-labeled action controls and section headers. | Use clear text hierarchy for section titles, select menus, action buttons, and theme controls. |

---

## 2. Runtime Environment and Technical Architecture

* **Self-Contained Execution:** Dashboards execute inside standard browsers and sandboxed webviews using standard HTML5, CSS custom properties, and client-side JavaScript.
* **Dynamic State Updates:** Interface components support incremental state updates, local file loading, and live data refresh without resetting scroll position or active camera state unless a new coordinate model is loaded.
* **Standalone Output (`dashboard.html`):** Build workflows emit a self-contained `dashboard.html` file in the target output directory or project root that can be opened directly from disk or served over HTTP with optional Server-Sent Events (SSE).
* **Empirical Data Rendering:** Views render computed structural coordinates, alignment metrics, provenance citations, and summary statistics from the active dataset.

---

## 3. Light, Dark, and System Theme Specifications

The dashboard natively supports light mode, high-contrast dark mode, and an automatic system mode that tracks the operating system appearance preference.

* **Three-State Segmented Theme Selector:** The theme switcher in the top navigation bar is a compact segmented control consisting of three monospaced text buttons placed side by side:
  * **Light Mode (`#theme-btn-light`):** Text label `Light` (`title="Light Theme"`, `aria-label="Light Theme"`, `data-theme-mode="light"`).
  * **Dark Mode (`#theme-btn-dark`):** Text label `Dark` (`title="Dark Theme"`, `aria-label="Dark Theme"`, `data-theme-mode="dark"`).
  * **System Mode (`#theme-btn-system`):** Text label `System` (`title="System Theme"`, `aria-label="System Theme"`, `data-theme-mode="system"`).
  * The active mode button applies the `active` class and `aria-pressed="true"`, and persists the selected mode (`light`, `dark`, or `system`) in `localStorage`.
* **Cross-Platform System Theme Detection:**
  * **Browser / Webview Layer:** When `system` mode is selected, resolve the active theme via `window.matchMedia('(prefers-color-scheme: dark)')` and register a `'change'` listener so operating system appearance updates propagate immediately.
  * **Host OS Layer (Backend Fallback):** When compiled by a local CLI script, inspect native desktop environment appearance settings (`GTK_THEME`, `COLORFGBG`, `gsettings`, or `kdeglobals`) and embed `system_theme` in the JSON payload as a fallback for restricted webviews.
* **Color Schemes:**
  * **Light Mode:**
    * Background: `#FAFAF9`
    * Cards & Surfaces: `#FFFFFF`
    * Borders: `#E7E5E4` (`1px solid`, `6px` to `8px` border-radius)
    * Text Primary: `#1C1917`
    * Text Secondary: `#787774`
  * **Dark Mode:**
    * Background: `#0D0E11`
    * Cards & Surfaces: `#16171D`
    * Borders: `#292A30` (`1px solid`, `6px` to `8px` border-radius)
    * Text Primary: `#ECE8E2`
    * Text Secondary: `#9A9895`

---

## 4. Dashboard Layout and Component Rules

* **Section Headings:** Keep component titles concise and descriptive of the underlying analytical domain.
* **Tables:** Render numerical and tabular records in a monospaced font (`DejaVu Sans Mono`, `ui-monospace`, or `monospace`) with `font-variant-numeric: tabular-nums` and `1px` horizontal row dividers.
* **3D and Quantitative Visualizations:**
  * Provide interactive 3D molecular coordinate rendering via WebGL paired with an automatic software 3D projection fallback when hardware WebGL is unavailable.
  * Use clear component titles such as `Structural Data Diagram` and `Homology Landscape`.
  * Keep coordinate projection and bond calculations deterministic for responsive rotation, zoom, and residue inspection.

---

## 5. Data Visualization Color Palette

Container surfaces remain neutral while quantitative charts, alignment tracks, and 3D molecular views apply a five-role semantic color palette legible against both light and dark surfaces:

### Semantic Color Tokens

| Color Role | Hex Code | Visual Application |
|---|---|---|
| **Primary Teal** | `#0D9488` | Primary data series, baseline coverage depth, alpha-helices, high-identity homologs. |
| **Secondary Amber** | `#D97706` | Secondary data series, coil/loop regions, partial or local motif matches. |
| **Comparison Indigo** | `#4F46E5` | Comparison series, beta-sheets, remote twilight-zone structural homologs. |
| **Validated Emerald** | `#059669` | High-probability fold validation (`Prob >= 0.90`) and upper confidence bounds. |
| **Alert Rose** | `#E11D48` | Low-confidence hits (`Prob < 0.50`), active residue selection anchors, and error states. |

### Visualization Guidelines

* **Continuous Scales:** Transition across `#0D9488`, `#4F46E5`, `#D97706`, and `#E11D48` (or `#059669` for confidence metrics) so residue and significance gradients remain distinct in both themes.
* **Molecular Representations:** Support secondary-structure ribbons, covalent sidechain sticks, spacefill spheres, and translucent molecular surfaces colored by sequence spectrum, secondary structure, B-factor / pLDDT, alignment depth, or chain.
* **Chart Backgrounds and Gridlines:** Keep SVG plot backgrounds transparent so they inherit the active light or dark surface token, using `#E7E5E4` in light mode and `#292A30` in dark mode for axis gridlines.

---

## 6. Verification Checklist

1. Does the build emit a self-contained `dashboard.html` file that operates directly from disk and over HTTP?
2. Does the theme control provide `#theme-btn-light`, `#theme-btn-dark`, and `#theme-btn-system` with real-time OS preference tracking?
3. Do Light Mode and Dark Mode apply the exact background, surface, border, and text hex tokens defined in Section 3?
4. Do all 2D SVG charts, alignment tracks, and 3D molecular views apply the semantic palette from Section 5 and update cleanly across theme changes?
5. Does the dashboard display empirical structure coordinates, alignment statistics, on-demand `.pdb` / `.cif` inspection tables, and provenance citations?