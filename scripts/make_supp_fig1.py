#!/usr/bin/env python3
"""
Supplementary Figure S1 - study design flowchart.

Reproduces the flowchart that was previously drawn by hand, so that it can be
regenerated whenever the underlying numbers change. Every count is read from
output/flowchart_numbers.csv (written by scripts/flowchart_numbers.R) rather than
typed in, which is what went wrong before: the hand-drawn version still showed
the pre-correction enrichment counts.

Two things make this look like the PowerPoint original rather than a plot:

  1. One data unit = one point. figsize is set to (x-range/72, y-range/72), so a
     box 22 units tall is 22 pt tall and a fontsize of 6 is 6 units. Box geometry
     and type sizes are then directly comparable, which they are not by default.
  2. Boxes are sized FROM their text. measure() renders each label, reads its
     extent, and the box is that extent plus a fixed padding - so no box is
     larger than its contents need, and corners use a generous rounding_size to
     match PowerPoint's rounded rectangle.

Colours sampled from the original figure so the reproduction matches:
    dataset boxes  #fff2cd    complex trait GWAS  #ffff00
    outcome GWAS   #f4cdcd    database/green      #d2e3cc
    grey panels    #ededed    filter boxes        #ffffff

Usage:  python3 scripts/make_supp_fig1.py
Output: figures/supp_fig1_revised.png  (600 dpi)
        figures/supp_fig1_revised.pdf  (vector)
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Nimbus Sans (a Helvetica clone) rather than matplotlib's DejaVu Sans: S2-S4 and the
# captions are Helvetica/Arial, and DejaVu made S1 the odd one out. fonttype 42 embeds
# the font as TrueType, which journals require. The >= and <= glyphs are used directly
# instead of mathtext ($\geq$), which would pull in a Computer Modern face as well.
matplotlib.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Nimbus Sans", "Helvetica", "Liberation Sans", "Arial"],
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

ROOT = os.environ.get("MRCOLOC_ROOT", "/home/mohd/mrcoloc")
CSV  = os.path.join(ROOT, "output", "flowchart_numbers.csv")

# ---- numbers, read from the pipeline -----------------------------------------
def load_numbers():
    n = {}
    with open(CSV) as fh:
        for row in csv.DictReader(fh):
            n[row["category"]] = row["value"]
    return n

N = load_numbers()
def i(k):  return int(float(N[k]))
def f(k):  return "{:,}".format(i(k))

# ---- palette (sampled from the original) -------------------------------------
C_DATA   = "#fff2cd"
C_YELLOW = "#ffff00"
C_PINK   = "#f4cdcd"
C_GREEN  = "#d2e3cc"
C_GREY   = "#ededed"
C_WHITE  = "#ffffff"
EDGE     = "#7f7f7f"
LINE     = "#595959"

# ---- canvas: 1 data unit == 1 pt ---------------------------------------------
X0, X1 = 104, 852
Y0, Y1 = 14, 548
fig, ax = plt.subplots(figsize=((X1 - X0) / 72.0, (Y1 - Y0) / 72.0))
ax.set_xlim(X0, X1); ax.set_ylim(Y1, Y0); ax.axis("off")   # y inverted: 0 at top
fig.canvas.draw()                                          # renderer for measure()
REND = fig.canvas.get_renderer()

PADX, PADY = 7, 5      # space between text and box edge, in pt
ROUND      = 4         # corner radius; PowerPoint's rounded rectangle is ~this

# All type sizes in one place, in pt (== data units). Raised from the first pass,
# where the two grey panels in particular were set too small to read once the page
# is scaled down to sit above its caption. Fixed-width boxes are checked for
# overflow at the end of the script, so these can be edited without guessing.
FS = {
    "stage":    8.0,   # the four stage labels down the left
    "dataset":  6.5,   # the eight pGWAS boxes (the tightest fit in the figure)
    "ctg":      7.5,   # complex trait GWAS
    "outcome":  7.0,   # the three outcome GWAS sources
    "mrtests":  7.5,   # 47.2 million MR tests
    "notes":    6.5,   # the two right-hand method/filter notes
    "green":    7.5,   # unique target-trait pairs
    "pharma":   7.5,   # match with Pharmaprojects
    "filter":   7.0,   # the centre column of filters
    "panel_h":  8.5,   # grey panel headings
    "panel_t":  7.0,   # grey panel body text
}
OVERFLOW = []          # (label, needed width, available width), filled by box()

def measure(text, fs, bold=False, lsp=1.3):
    """width, height of a rendered label, in data units (== pt)."""
    t = ax.text(0, 0, text, fontsize=fs, linespacing=lsp,
                fontweight="bold" if bold else "normal")
    bb = t.get_window_extent(renderer=REND)
    t.remove()
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1))
    return abs(x1 - x0), abs(y1 - y0)

def box(text, y, cx=None, x=None, w=None, h=None, fs=6, bold=False, fc=C_WHITE,
        align="center", lsp=1.3, padx=None, pady=None):
    """Draw a box sized to its text. Pass w/h to force a dimension, which is how
    rows of boxes with unequal line counts are kept the same size. (x, y, w, h)."""
    padx = PADX if padx is None else padx
    pady = PADY if pady is None else pady
    tw, th = measure(text, fs, bold, lsp)
    if w is not None and tw + 2 * padx > w:
        OVERFLOW.append((text.split("\n")[0], tw + 2 * padx, w))
    bw = tw + 2 * padx if w is None else w
    bh = th + 2 * pady if h is None else h
    bx = cx - bw / 2 if cx is not None else x
    ax.add_patch(FancyBboxPatch((bx, y), bw, bh,
                                boxstyle="round,pad=0,rounding_size=%g" % ROUND,
                                linewidth=0.8, edgecolor=EDGE, facecolor=fc, zorder=2))
    if align == "center":
        ax.text(bx + bw / 2, y + bh / 2, text, ha="center", va="center",
                fontsize=fs, fontweight="bold" if bold else "normal",
                linespacing=lsp, zorder=3)
    else:
        ax.text(bx + padx, y + bh / 2, text, ha="left", va="center",
                fontsize=fs, fontweight="bold" if bold else "normal",
                linespacing=lsp, zorder=3)
    return bx, y, bw, bh

def panel(x, y, w, h):
    ax.add_patch(Rectangle((x, y), w, h, linewidth=0.8, edgecolor=EDGE,
                           facecolor=C_GREY, zorder=1))

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=7, linewidth=0.8, color=LINE, zorder=4))

def elbow(x1, y1, x2, y2):
    """down, across, then into the target with an arrowhead"""
    mid = y1 + (y2 - y1) * 0.55
    ax.plot([x1, x1], [y1, mid], color=LINE, lw=0.8, zorder=1)
    ax.plot([x1, x2], [mid, mid], color=LINE, lw=0.8, zorder=1)
    arrow(x2, mid, x2, y2)

SPINE = 533            # the vertical axis of the flow

# ---- stage labels (left column) ----------------------------------------------
box("pGWAS",             y=28,  x=118, fs=FS["stage"], bold=True, fc=C_DATA)
box("OUTCOME GWAS",      y=120, x=118, fs=FS["stage"], bold=True, fc=C_PINK)
box("MR-COLOC ANALYSES", y=222, x=118, fs=FS["stage"], bold=True, fc=C_WHITE)
box("DATABASE\nCONSTRUCTION &\nENRICHMENT ANALYSES",
                         y=300, x=118, fs=FS["stage"], bold=True, fc=C_GREEN, lsp=1.4)

# ---- row 1: the eight proteogenomic datasets --------------------------------
datasets = [("SCALLOP_2020", "n = 90 proteins"), ("PIETZNER_2020", "n = 186"),
            ("HILLARY_2019", "n = 92"),          ("SUN_2018", "n = 3,622"),
            ("SUHRE_2017", "n = 1,124"),         ("OLLI_2017", "n = 41"),
            ("FOLKERSEN_2017", "n = 83"),        ("UKBPPP_2023", "n = 2,923")]
# Width is derived from the widest label rather than fixed, because the row has to
# tile eight equal boxes across the frame and FOLKERSEN_2017 is the binding
# constraint - at a fixed 77 it did not fit even at the smaller original size. A
# tighter horizontal pad is used here than elsewhere; this row is the only place
# where labels come close to their box edges.
DS_GAP, DS_PADX = 3, 3
ds_labels = ["%s\n(%s)" % (name, nn) for name, nn in datasets]
dw = max(measure(t, FS["dataset"])[0] for t in ds_labels) + 2 * DS_PADX
row_w = 8 * dw + 7 * DS_GAP
# Right-align the row to the frame; xlim clips, so a row that overruns X1 would
# lose the last box. 176 is where the pGWAS stage label ends.
DS_X0 = max(176, X1 - 6 - row_w)
if DS_X0 + row_w > X1:
    print("  !! dataset row is %.1f wide and cannot fit the frame" % row_w)
x, centres = DS_X0, []
for lab in ds_labels:
    _, _, _, dh = box(lab, y=28, x=x, w=dw, fs=FS["dataset"], fc=C_DATA, padx=DS_PADX)
    centres.append(x + dw / 2); x += dw + DS_GAP
DROW_BOT = 28 + dh

# ---- row 2: complex trait GWAS ----------------------------------------------
CTG_Y = 78
for cx in centres:
    elbow(cx, DROW_BOT, SPINE, CTG_Y)
box("Complex trait GWAS\n(n = %s)" % f("gwas_total"),
    y=CTG_Y, cx=SPINE, fs=FS["ctg"], bold=True, fc=C_YELLOW)

# ---- row 3: outcome GWAS sources -------------------------------------------
MR_Y = 218
outs = [("GWAS catalog",           f("gwas_catalog"),      418),
        ("pan-UKB +\nNEALE/SAIGE", f("gwas_panukb_neale"), SPINE),
        ("Finngen (r5 + r12)",     f("gwas_finngen"),      648)]
labels = ["%s\n(n = %s)" % (label, nn) for label, nn, _ in outs]
OH = max(measure(t, FS["outcome"])[1] for t in labels) + 2 * PADY   # one height for the row,
for lab, (_, _, ox) in zip(labels, outs):               # pan-UKB wraps to 3 lines
    box(lab, y=120, cx=ox, w=110, h=OH, fs=FS["outcome"], fc=C_PINK)
    arrow(ox, 120 + OH, SPINE, MR_Y)

# ---- row 4: MR tests + method notes ----------------------------------------
_, _, _, mrh = box("%.1f million MR tests" % (i("mr_total_tests") / 1e6),
                   y=MR_Y, cx=SPINE, fs=FS["mrtests"])
box("Mendelian randomization\n"
    "  -  Genome-wide (p \u2264 5e-8)\n"
    "  -  r2 < 0.05\n"
    "  -  HEIDI-outlier-flag = T\n"
    "Genetic colocalization\n"
    "  -  cis & trans H4 \u2265 0.8",
    y=196, x=684, w=152, fs=FS["notes"], align="left", lsp=1.45)

GREEN_Y = 300
arrow(SPINE, MR_Y + mrh, SPINE, GREEN_Y)

# ---- row 5: target-trait pairs + additional filtering ----------------------
_, _, _, gh = box("%s unique target-trait pairs\n(%s colocalizing)"
                  % (f("ttpairs_total"), f("ttpairs_with_any_coloc")),
                  y=GREEN_Y, cx=SPINE, fs=FS["green"], fc=C_GREEN)
box("Additional filtering\n"
    "  -  MR p-value < 0.05/47M tests\n"
    "  -  Remove medically irrelevant\n     traits (n = %s)" % f("traits_excluded"),
    y=300, x=684, w=152, fs=FS["notes"], align="left", lsp=1.45)

PHARMA_Y = 388
arrow(SPINE, GREEN_Y + gh, SPINE, PHARMA_Y)

# ---- row 6: Pharmaprojects --------------------------------------------------
_, _, _, ph = box("Match with Pharmaprojects\n(Minikel et al)",
                  y=PHARMA_Y, cx=SPINE, fs=FS["pharma"], bold=True)

# ---- centre column of filters, built first so the panels can match its height
FILT_Y, FILT_GAP = 434, 5
fy = FILT_Y
for label in ["MeSH similarity \u2265 0.8", "L2G share \u2265 0.5",
              "Restricted to genes with\nmeasured proteins", "Phase 1 baseline"]:
    _, _, _, fh = box(label, y=fy, cx=SPINE, w=152, fs=FS["filter"])
    fy += fh + FILT_GAP
FILT_BOT = fy - FILT_GAP
arrow(SPINE, PHARMA_Y + ph, SPINE, FILT_Y)

# ---- row 7: the two analysis panels ----------------------------------------
PAN_Y, PAN_PAD = 418, 8
left_txt  = ("-    %s drug targets overlap\n-    %s target-indication pairs"
             % (f("path1_drug_targets"), f("path1_ti_matched")))
right_txt = ("-    %s pQTL-supported T-I pairs (S[G])\n"
             "        -    %s pQTL-supported T-I pairs\n"
             "             from Phase 1 > Launch\n"
             "             (A[G])\n"
             "-    %s unsupported (S![G])\n"
             "        -    %s unsupported T-I pairs\n"
             "             from Phase > Launch\n"
             "             (A![G])"
             % (f("enrichment_supported_phase1"), f("enrichment_supported_launched"),
                f("enrichment_unsupported_phase1"), f("enrichment_unsupported_launched")))

# Panel height is whichever is taller: the filter column it brackets, or its own
# contents. Previously it tracked the filter column alone, which is why raising the
# body type risked overrunning the grey box.
HEAD_H  = measure("pQTL enrichment analysis", FS["panel_h"], bold=True)[1]
body_h  = max(measure(left_txt,  FS["panel_t"], lsp=2.00)[1],
              measure(right_txt, FS["panel_t"], lsp=1.45)[1])
PAN_H   = max(FILT_BOT - PAN_Y + 6, PAN_PAD + HEAD_H + 10 + body_h + PAN_PAD)
panel(232, PAN_Y, 220, PAN_H)          # left: Supplementary Table 16
panel(612, PAN_Y, 226, PAN_H)          # right: pQTL enrichment analysis

for px, heading in [(242, "Supplementary Table 16"), (622, "pQTL enrichment analysis")]:
    ax.text(px, PAN_Y + PAN_PAD, heading, fontsize=FS["panel_h"],
            fontweight="bold", va="top", zorder=3)
BODY_Y = PAN_Y + PAN_PAD + HEAD_H + 10

ax.text(258, BODY_Y + 6, left_txt,  fontsize=FS["panel_t"], va="top", zorder=3,
        linespacing=2.00)
ax.text(628, BODY_Y, right_txt, fontsize=FS["panel_t"], va="top", zorder=3,
        linespacing=1.45)

if PAN_Y + PAN_H > Y1:
    print("  !! grey panels overrun the canvas: %.1f > %.1f" % (PAN_Y + PAN_H, Y1))

# ---- save --------------------------------------------------------------------
for out, kw in [("figures/supp_fig1_revised.png", dict(dpi=600)),
                ("figures/supp_fig1_revised.pdf", {})]:
    fig.savefig(os.path.join(ROOT, out), bbox_inches="tight", pad_inches=0.04,
                facecolor="white", **kw)
    print("  -> %s" % out)

if OVERFLOW:
    print("\n  !! text wider than its fixed box (reduce FS or widen the box):")
    for lab, need, have in OVERFLOW:
        print("     %-34s needs %.1f, has %.1f" % (lab, need, have))
else:
    print("\n  all fixed-width boxes fit their text")

print("\n  numbers used (from output/flowchart_numbers.csv):")
for k in ["gwas_total","mr_total_tests","ttpairs_total","ttpairs_with_any_coloc",
          "traits_excluded","path1_drug_targets","path1_ti_matched",
          "enrichment_supported_phase1","enrichment_supported_launched",
          "enrichment_unsupported_phase1","enrichment_unsupported_launched"]:
    print("    %-34s %s" % (k, N[k]))
