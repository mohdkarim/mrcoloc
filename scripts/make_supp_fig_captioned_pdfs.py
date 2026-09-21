#!/usr/bin/env python3
"""
Build one captioned PDF per supplementary figure, matching what was submitted.

The submitted captioned pages (figures/submitted_supplementary_figures_with_captions/)
were assembled by hand in PowerPoint, one at a time, and are inconsistent as a
result: left margins range from 17 to 63 pt and caption type from 9 to 11 pt across
the four. This script keeps what was deliberate - the 720 x 540 pt page, the
Helvetica/Arial caption with a bold "Supplementary Figure Sn." lead-in, figure above
caption - and standardises the rest.

Layout: the caption is typeset into a box first and measured, then the figure is
given whatever vertical space is left and scaled to fit it (keepaspectratio), so
captions of different lengths all produce a full page with no overflow.

CAPTION TEXT is carried over verbatim from the submitted pages. S2's caption is the
only one that asserts numbers - the 6 + 3 excluded therapeutic areas - and that was
re-derived from the corrected data before reuse: still 6 with no pQTL-supported
pairs at Phase I (endocrine, immune, infection, ophthalmology, other, psychiatry)
and 3 with supported pairs but no launches (neurology, oncology, signs/symptoms),
8 included. The point estimates inside S2 did move (cardiovascular 3.08 -> 2.31,
metabolic 4.46 -> 3.72), but the caption does not quote them.

Usage:  python3 scripts/make_supp_fig_captioned_pdfs.py
Output: figures/supplementary_figures_with_captions/figS{1..4}_*.pdf
"""
import glob, os, re, shutil, subprocess, sys

ROOT    = os.environ.get("MRCOLOC_ROOT", "/home/mohd/mrcoloc")
FIGDIR  = os.path.join(ROOT, "figures")
OUTDIR  = os.path.join(FIGDIR, "supplementary_figures_with_captions")
TEXBIN  = "/home/mohd/.TinyTeX/bin/x86_64-linux"
GS      = "/usr/bin/gs"

# page and type, uniform across all four
PAPER_W, PAPER_H = "720bp", "540bp"   # bp = 1/72 in, the PDF point; "pt" is 1/72.27
MARGIN           = "28bp"
CAP_SIZE, CAP_LEAD = 11.0, 13.2
GAP              = "14pt"      # between figure and caption

# Ceiling on a forest plot's on-page body text. Panels are fitted to the page, which
# is what makes the type as large as the space allows - but fitting alone made it
# uneven: S4 is drawn with base_size 14 in R against 11 for S2 and S3, so on the same
# page its text reached ~18 pt while S2 sat near 10 pt. Capping at 14 pt pulls the
# outlier back without shrinking it so far that a two-row plot floats on a mostly
# empty page. S1 has no single base size (a flowchart of differing label sizes), so
# it is fitted with no cap.
#
# The underlying inconsistency is in mrcoloc_paper_2025_supp_figures.R, where S4 asks
# for base_size 14 and the others 11. Worth aligning there if that script is revised.
MAX_BODY_PT = 14.0

# ---- the four figures. source = the regenerated panel, not the submitted page ---
FIGURES = [
 ("figS1_flowchart", "supp_fig1_revised.pdf", "S1", None,   # flowchart: a range of label sizes, so fit to page
  "Study design flowchart illustrating the two analytical paths: (1) systematic "
  "identification of pQTL-supported target-trait pairs through Mendelian randomization "
  "of eight proteogenomic datasets against \\textasciitilde{}8,000 complex trait GWAS, and (2) "
  "enrichment analysis comparing clinical success rates of pQTL-supported versus "
  "unsupported target-indication pairs using the Minikel et al. framework. Numbers at "
  "each step indicate the count of datasets, traits, associations, or T-I pairs "
  "retained after applying the indicated filters."),

 ("figS2_enrichment_by_ta", "figS2_enrichment_by_ta.pdf", "S2", 11.0,   # forest_plot(base_size = 11)
  "Phase I-Launch relative success (RS) of pQTL-supported versus unsupported "
  "target-indication (T-I) pairs by therapeutic area. RS was computed within each "
  "therapeutic area using the Katz log method. Therapeutic areas were included if they "
  "had at least one pQTL-supported T-I pair entering Phase I and at least one "
  "pQTL-supported pair that reached launch (RS \\textgreater{} 0). Nine therapeutic areas were "
  "excluded: six with no pQTL-supported T-I pairs at Phase I (endocrine, immune, "
  "infection, ophthalmology, other, psychiatry) and three with pQTL-supported pairs but "
  "no launches, yielding RS = 0 (neurology, oncology, signs/symptoms). The background "
  "was restricted to T-I pairs where the target was measured on at least one proteomic "
  "platform (Olink or SomaScan). Supported = A[G]/S[G] = launched/Phase I among "
  "pQTL-supported pairs; Unsupported = A[!G]/S[!G] = launched/Phase I among "
  "pQTL-unsupported pairs. Error bars: 95\\% CI (Katz log method)."),

 ("figS3_enrichment_by_phase", "figS3_enrichment_by_phase.pdf", "S3", 11.0,   # forest_plot(base_size = 11)
  "Relative success (RS) of pQTL-supported versus unsupported target-indication (T-I) "
  "pairs by clinical phase transition. RS was computed separately for each individual "
  "phase transition (Preclinical \\textgreater{} I, I \\textgreater{} II, II \\textgreater{} III, III \\textgreater{} Launch) and for "
  "cumulative transitions (Preclinical \\textgreater{} Launch, I \\textgreater{} Launch). The denominator for "
  "each transition includes only T-I pairs that entered the starting phase. The "
  "background was restricted to T-I pairs where the target was measured on at least one "
  "proteomic platform (Olink or SomaScan). Supported = A[G]/S[G] = succeeded/entered "
  "among pQTL-supported pairs; Unsupported = A[!G]/S[!G] = succeeded/entered among "
  "pQTL-unsupported pairs. Error bars: 95\\% CI (Katz log method)."),

 ("figS4_background_comparison", "figS4_background_comparison.pdf", "S4", 14.0,   # forest_plot(base_size = 14) - the odd one out
  "Sensitivity analysis comparing Phase I-Launch relative success (RS) of "
  "pQTL-supported T-I pairs under two background definitions: (1) measured proteins "
  "background, restricted to T-I pairs where the target was assayed on at least one "
  "proteomic platform used in the pQTL studies (Olink or SomaScan), matching the main "
  "analysis; and (2) full Minikel background, using all T-I pairs from Minikel et al. "
  "without restriction to measured proteins. pQTL-supported T-I pairs are identical in "
  "both analyses; only the unsupported comparator group changes. A[G]/S[G] = "
  "launched/Phase I among pQTL-supported pairs; A[!G]/S[!G] = launched/Phase I among "
  "pQTL-unsupported pairs. Error bars: 95\\% CI (Katz log method)."),
]

TEMPLATE = r"""\documentclass[10pt]{article}
\usepackage[paperwidth=%(pw)s,paperheight=%(ph)s,margin=%(margin)s]{geometry}
\usepackage[T1]{fontenc}
\usepackage[scaled]{helvet}          %% Helvetica: metrically Arial, as submitted
\renewcommand{\familydefault}{\sfdefault}
\usepackage{graphicx}
\pagestyle{empty}
\setlength{\parindent}{0pt}
\graphicspath{{%(figdir)s/}}

\newsavebox{\capbox}
\newlength{\figh}

\begin{document}
%% Typeset the caption first, then give the figure the space that is left. Doing it
%% this way means a long caption shrinks the figure instead of overflowing the page.
\savebox{\capbox}{\parbox[b]{\textwidth}{%%
  \fontsize{%(capsize).1f}{%(caplead).1f}\selectfont\raggedright
  \textbf{Supplementary Figure %(label)s.} %(caption)s}}
\setlength{\figh}{\dimexpr\textheight-\ht\capbox-\dp\capbox-%(gap)s-4bp\relax}
\typeout{MEASURED figh=\the\figh\space textwidth=\the\textwidth}

\vbox to \textheight{\boxmaxdepth=0pt \offinterlineskip
  \vfil
  \hbox to \textwidth{\hfil
    \includegraphics[trim=%(trim)s,clip,scale=%(scale).5f]{%(src)s}\hfil}
  \vfil
  \vskip %(gap)s
  \usebox{\capbox}}
\end{document}
"""


def ink_trim(path):
    """`trim=` values that crop a panel to its ink, via Ghostscript's bbox device.

    All four panels carry 8-24 pt of white margin of their own. Left in, that margin
    is scaled along with the figure and eats into the space available on the page, so
    the plot text ends up smaller than it needs to be. Cropping here rather than
    regenerating the panels keeps the sources untouched.
    """
    info = subprocess.run(["pdfinfo", path], capture_output=True, text=True).stdout
    W, H = [float(v) for v in re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", info).groups()]
    p = subprocess.run([GS, "-dQUIET", "-dNOPAUSE", "-dBATCH", "-sDEVICE=bbox", path],
                       capture_output=True, text=True)
    m = re.search(r"%%HiResBoundingBox:\s+([\d.-]+) ([\d.-]+) ([\d.-]+) ([\d.-]+)",
                  p.stderr + p.stdout)
    if not m:
        return "0 0 0 0", W, H
    x0, y0, x1, y1 = [float(v) for v in m.groups()]
    pad = 2.0                                     # don't shave the ink itself
    l, b = max(0.0, x0 - pad), max(0.0, y0 - pad)
    r, t = max(0.0, W - x1 - pad), max(0.0, H - y1 - pad)
    return "%.1fbp %.1fbp %.1fbp %.1fbp" % (l, b, r, t), W - l - r, H - b - t


def embed_fonts(path):
    """Ghostscript pass to embed every font.

    The R-generated panels (figS2-S4) reference Helvetica WITHOUT embedding it - the
    pdf() device treats it as a base-14 font and assumes the viewer has it. Journals
    reject that, and it is why the submitted panels were at the mercy of whatever
    Helvetica the typesetter had. Ghostscript substitutes Nimbus Sans, which is
    metrically identical, and embeds it. Vector content is preserved.
    """
    tmp = path + ".embed"
    p = subprocess.run([GS, "-dNOPAUSE", "-dBATCH", "-dQUIET", "-sDEVICE=pdfwrite",
                        "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
                        "-dSubsetFonts=true", "-dAutoRotatePages=/None",
                        "-sOutputFile=" + tmp, path], capture_output=True, text=True)
    if p.returncode != 0 or not os.path.exists(tmp):
        print("     !! font embedding failed: %s" % p.stderr.strip()[:200])
        return False
    os.replace(tmp, path)
    q = subprocess.run(["pdffonts", path], capture_output=True, text=True).stdout
    unembedded = [l.split()[0] for l in q.splitlines()[2:]
                  if l.strip() and l.split()[-4:-3] == ["no"]]
    if unembedded:
        print("     !! still unembedded: %s" % ", ".join(unembedded))
        return False
    return True


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    env = dict(os.environ, PATH=TEXBIN + os.pathsep + os.environ.get("PATH", ""))
    tmp = os.path.join(OUTDIR, "_build")
    os.makedirs(tmp, exist_ok=True)
    rc = 0

    for name, src, label, base_pt, caption in FIGURES:
        src_path = os.path.join(FIGDIR, src)
        if not os.path.exists(src_path):
            print("  !! missing source: %s" % src); rc = 1; continue

        trim, trim_w, trim_h = ink_trim(src_path)
        texf = os.path.join(tmp, "page_" + name + ".tex")

        def build(scale):
            with open(texf, "w") as fh:
                fh.write(TEMPLATE % dict(pw=PAPER_W, ph=PAPER_H, margin=MARGIN,
                                         figdir=FIGDIR, capsize=CAP_SIZE,
                                         caplead=CAP_LEAD, gap=GAP, trim=trim,
                                         scale=scale, label=label, caption=caption,
                                         src=src))
            return subprocess.run([os.path.join(TEXBIN, "pdflatex"),
                                   "-interaction=nonstopmode", "-halt-on-error",
                                   "-output-directory", tmp, texf],
                                  capture_output=True, text=True, env=env)

        # Pass 1 measures the caption, which is what determines how much height is
        # left; the figure scale cannot be worked out before that is known.
        p = build(1.0)
        m = re.search(r"MEASURED figh=([\d.]+)pt textwidth=([\d.]+)pt", p.stdout)
        if not m:
            print("  !! could not measure the caption box for %s" % name); rc = 1; continue
        figh, tw = float(m.group(1)) / 1.00375, float(m.group(2)) / 1.00375
        fit = min(tw / trim_w, figh / trim_h)
        scale = fit if base_pt is None else min(fit, MAX_BODY_PT / base_pt)

        p = build(scale)                                  # pass 2, at the real scale
        built = os.path.join(tmp, "page_" + name + ".pdf")
        if p.returncode != 0 or not os.path.exists(built):
            print("  !! pdflatex failed for %s" % name)
            print("\n".join(p.stdout.splitlines()[-25:])); rc = 1; continue

        final = os.path.join(OUTDIR, name + ".pdf")
        shutil.move(built, final)
        ok = embed_fonts(final)
        bad = [l for l in p.stdout.splitlines() if "Overfull" in l or "Underfull" in l]
        if base_pt is None:
            eff, why = "fitted to page", ""
        else:
            eff = "body text %.1f pt" % (scale * base_pt)
            why = "  [page-limited]" if scale < MAX_BODY_PT / base_pt - 1e-9 else "  [capped]"
        print("  -> %-30s %-15s at %3.0f%% of native%s%s%s" % (
              name + ".pdf", "fonts embedded" if ok else "FONT PROBLEM",
              100 * scale, "  " + eff, why,
              "  (%d box warning(s))" % len(bad) if bad else ""))

    shutil.rmtree(tmp, ignore_errors=True)
    print("\n  output: figures/supplementary_figures_with_captions/")
    return rc

sys.exit(main())
