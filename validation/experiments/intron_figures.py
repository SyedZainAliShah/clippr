"""Draw the intron figures from `figure_data.json`. Nothing is computed here.

Every number comes from the JSON the generator wrote, so a figure cannot drift from the
measurement behind it. Run `intron_figure_data.py` first.

Output is SVG, in the wiki's own palette, sized for a 5xl content column.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "work/introns/figure_data.json")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "work/introns/figures")
OUT.mkdir(parents=True, exist_ok=True)
F = json.loads(DATA.read_text(encoding="utf-8"))

# The wiki's daisyUI theme, from src/containers/App/App.css. Figures that use these read as
# part of the page rather than as something pasted into it.
GREEN, GREEN_L, GREEN_D = "#65b93d", "#cfebc1", "#428725"
PINK, PINK_L, PINK_D = "#ff66c4", "#e0b0d0", "#a8337e"
CYAN, CYAN_L, CYAN_D = "#5cd6ff", "#d0f3ff", "#007aa3"
INK, MUTED, FAINT = "#111111", "#666666", "#dddddd"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",       # keep text as text, so it stays crisp and selectable
    "figure.facecolor": "white",
})


def save(fig, name):
    path = OUT / (name + ".svg")
    fig.savefig(path, format="svg", bbox_inches="tight", transparent=False)
    plt.close(fig)
    print("  " + str(path))


# ------------------------------------------------------------------ 1. the two introns
def fig_gc_band():
    """Why the weaker intron is the one we can use: the stronger one leaves the band."""
    I = F["introns"]["introns"]
    lo, hi = F["introns"]["band"]
    win = F["introns"]["window_nt"]
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.6), sharey=True,
                             gridspec_kw={"hspace": 0.45})

    for ax, (name, colour, dark) in zip(axes, (("rbcS2i1", GREEN, GREEN_D),
                                               ("LHCBM1i2", PINK, PINK_D))):
        v = I[name]
        trace = v["gc_trace"]
        x = list(range(len(trace)))

        # The band is the rule. Draw it first so the trace sits inside or breaks out of it.
        ax.axhspan(lo, hi, color=GREEN_L, alpha=0.45, lw=0, zorder=0)
        ax.axhline(hi, color=MUTED, lw=1, ls="--", zorder=1)
        ax.plot(x, trace, color=colour, lw=2, zorder=3, solid_capstyle="round")

        # Mark every window that breaks the ceiling, not just the worst one.
        over = [(i, g) for i, g in enumerate(trace) if g > hi or g < lo]
        if over:
            ax.scatter([i for i, _ in over], [g for _, g in over], s=26, color=PINK_D,
                       zorder=4, linewidths=0)

        worst, at = v["worst_window_gc"], v["worst_window_at"]
        ax.scatter([at], [worst], s=70, facecolor="white", edgecolor=dark,
                   linewidths=2, zorder=5)

        verdict = ("inside the band at every window"
                   if v["passes"] else
                   str(v["windows_breaching"]) + " of " + str(v["windows"])
                   + " windows breach the ceiling")
        ax.set_title(name + "  ·  " + str(v["length_nt"]) + " nt  ·  "
                     + str(v["fold_enhancement"]) + "× enhancement  ·  "
                     + verdict,
                     loc="left", fontsize=11.5, color=dark, pad=9)

        note = ("worst window " + format(worst, ".3f")
                + ("  —  " + format(v["margin_to_ceiling"], ".3f") + " under the ceiling"
                   if v["passes"] else
                   "  —  " + format(-v["margin_to_ceiling"], ".3f") + " over the ceiling"))
        # Park the label in the empty lower half so it never crosses the trace, which it
        # did when it was offset from the point itself.
        ax.annotate(note, xy=(at, worst), xytext=(at + len(trace) * 0.08, 0.435),
                    fontsize=9.5, color=dark,
                    arrowprops={"arrowstyle": "-", "color": dark, "lw": 1,
                                "connectionstyle": "angle,angleA=0,angleB=90,rad=4"})

        ax.set_xlim(-3, max(len(trace), 100) + 3)
        ax.set_ylim(0.38, 0.74)
        ax.set_xlabel("start of " + str(win) + "-nt window (nt)", fontsize=9.5)
        ax.set_ylabel("GC", fontsize=9.5)

    axes[0].text(0.995, 0.06, "allowed band " + format(lo, ".2f") + "–"
                 + format(hi, ".2f"), transform=axes[0].transAxes, ha="right",
                 fontsize=9, color=GREEN_D)
    fig.suptitle("An intron has to pass the same GC rule as everything else we order",
                 x=0.012, ha="left", fontsize=13.5, y=1.015)
    save(fig, "01_intron_gc_band")


# -------------------------------------------------------- 2. where the fidelity goes
def fig_waterfall():
    """The budget is spent before the introns arrive."""
    D = F["decomposition"]
    lvl = D["levels"]["level0"]
    floor = F["architectures"]["floor"]
    ceiling = lvl["acceptor_alone"]
    final = lvl["with_introns"]
    n = len(D["chosen"])

    backbone_cost = 1.0 - ceiling
    intron_cost = ceiling - final
    headroom = ceiling - floor

    # Two panels, because the whole decision happens inside a 0.07-wide strip that a full
    # 0-1 axis flattens into nothing. The left panel keeps the honest scale; the right one
    # magnifies the only region anybody is arguing about.
    fig, (ax, az) = plt.subplots(1, 2, figsize=(11, 4.8),
                                 gridspec_kw={"width_ratios": [1.55, 1], "wspace": 0.3})

    bars = [
        ("A perfect\nreaction", 0.0, 1.0, FAINT, INK),
        ("Empty\nbackbone", 1.0 - backbone_cost, backbone_cost, PINK, PINK_D),
        (str(n) + " intron\njunctions", final, intron_cost, CYAN, CYAN_D),
        ("What is\nleft", 0.0, final, GREEN, GREEN_D),
    ]
    for i, (_, bottom, height, colour, dark) in enumerate(bars):
        ax.bar(i, height, bottom=bottom, width=0.56, color=colour,
               edgecolor=dark, linewidth=1.2, zorder=3)
        if i in (1, 2):
            ax.annotate("−" + format(height, ".3f"),
                        xy=(i, bottom + height / 2), ha="center", va="center",
                        fontsize=10.5, color=dark, fontweight="bold", zorder=5)
        else:
            ax.annotate(format(height, ".3f"), xy=(i, bottom + height + 0.018),
                        ha="center", fontsize=11.5, color=dark, fontweight="bold", zorder=5)

    for i, y in ((0, 1.0), (1, ceiling), (2, final)):
        ax.plot([i + 0.28, i + 0.72], [y, y], color=MUTED, lw=1, ls=":", zorder=2)

    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([b[0] for b in bars], fontsize=9.5)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("predicted assembly fidelity", fontsize=10)
    ax.set_title("Where the fidelity goes", loc="left", fontsize=12, pad=8)

    # ---- the zoom: ceiling, floor and result, with the two quantities that decide it
    az.add_patch(Rectangle((0.08, floor), 0.84, ceiling - floor, color=GREEN_L,
                           alpha=0.5, lw=0, zorder=0))
    for y, colour, label in ((ceiling, PINK_D, "backbone ceiling " + format(ceiling, ".4f")),
                             (floor, MUTED, "our floor " + format(floor, ".2f")),
                             (final, GREEN_D, "with " + str(n) + " introns "
                              + format(final, ".4f"))):
        az.axhline(y, color=colour, lw=1.6,
                   ls="--" if y == floor else "-", zorder=3)
        az.annotate(label, xy=(0.085, y + (-0.0045 if y == floor else 0.0022)),
                    fontsize=9.5, color=colour,
                    va="top" if y == floor else "bottom")

    # Headroom against bill. These two arrows are the argument.
    az.annotate("", xy=(0.30, ceiling), xytext=(0.30, floor),
                arrowprops={"arrowstyle": "<->", "color": GREEN_D, "lw": 1.4})
    az.annotate("room above\nthe floor\n" + format(headroom, ".4f"),
                xy=(0.335, (ceiling + floor) / 2), fontsize=9.5, color=GREEN_D, va="center")
    az.annotate("", xy=(0.72, ceiling), xytext=(0.72, final),
                arrowprops={"arrowstyle": "<->", "color": CYAN_D, "lw": 1.4})
    az.annotate("what " + str(n) + " introns\ncost\n" + format(intron_cost, ".4f"),
                xy=(0.755, (ceiling + final) / 2), fontsize=9.5, color=CYAN_D, va="center")

    az.set_xlim(0, 1)
    az.set_ylim(final - 0.018, ceiling + 0.018)
    az.set_xticks([])
    az.set_ylabel("predicted assembly fidelity", fontsize=10)
    az.set_title("The cost is larger than the room", loc="left", fontsize=12, pad=8)

    fig.suptitle("At level 0 the backbone costs more than the introns do  ·  "
                 + "9S design, " + F["matrix"] + ", vector " + D["vector"].split(",")[0],
                 x=0.012, ha="left", fontsize=13, y=1.02)
    save(fig, "02_fidelity_waterfall")


# ------------------------------------------------------------------- 3. the cliff
def fig_cliff():
    """What introns cost each architecture, and why the big ones were never free to choose."""
    A = F["architectures"]
    floor, cap = A["floor"], A["oligo_cap_nt"]
    archs = list(A["architectures"])
    fig, ax = plt.subplots(figsize=(9.2, 5))

    ax.axhspan(0, floor, color=PINK_L, alpha=0.18, lw=0, zorder=0)
    ax.axhline(floor, color=PINK_D, lw=1.5, ls="--", zorder=4)
    # Right-aligned and above the line: the left side is where 9S crowds the floor, which
    # is exactly the region this label kept landing in.
    ax.annotate("our floor " + format(floor, ".2f"), xy=(len(archs) - 0.38, floor - 0.02),
                ha="right", va="top", fontsize=9.5, color=PINK_D)

    for i, arch in enumerate(archs):
        v = A["architectures"][arch]
        base, got, k = v["baseline_fidelity"], v["chosen_fidelity"], len(v["chosen"])
        good = got >= floor
        colour = GREEN_D if good else PINK_D

        # A dumbbell: where the reaction starts without introns, where it lands with them.
        ax.plot([i, i], [got, base], color=MUTED, lw=2, zorder=2, solid_capstyle="round")
        ax.scatter([i], [base], s=130, facecolor="white", edgecolor=MUTED,
                   linewidths=2, zorder=5)
        ax.scatter([i], [got], s=170, color=colour, zorder=5, linewidths=0)

        ax.annotate(format(base, ".3f"), xy=(i + 0.12, base), fontsize=9.5,
                    color=MUTED, va="center")
        ax.annotate(format(got, ".4f"), xy=(i + 0.12, got), fontsize=11,
                    color=colour, va="center", fontweight="bold")

        # Whether the count was chosen or imposed is the point of the figure, so it sits
        # with the result rather than under the axis where the tick labels live.
        forced = len(v["fidelity_by_count"]) == 1 or min(
            int(n) for n in v["fidelity_by_count"]) == k
        # The rightmost label would run off the canvas, so it reads back toward the centre.
        last = i == len(archs) - 1
        ax.annotate(str(k) + " introns" + (" (the minimum)" if forced else " (best of "
                    + str(len(v["fidelity_by_count"])) + ")"),
                    xy=(i - 0.12 if last else i + 0.12, got - 0.05), fontsize=9,
                    ha="right" if last else "left",
                    color=PINK_D if forced else MUTED, va="center")

    ax.set_xticks(range(len(archs)))
    ax.set_xticklabels([a + "\n" + str(A["architectures"][a]["cds_nt"]) + " nt, "
                        + str(A["architectures"][a]["protein_aa"]) + " aa" for a in archs],
                       fontsize=10.5)
    ax.set_xlim(-0.5, len(archs) - 0.3)
    ax.set_ylim(0, 0.92)
    ax.set_ylabel("predicted assembly fidelity", fontsize=10)
    ax.scatter([], [], s=130, facecolor="white", edgecolor=MUTED, linewidths=2,
               label="without introns")
    ax.scatter([], [], s=150, color=PINK_D, linewidths=0, label="with introns")
    ax.legend(frameon=False, fontsize=9.5, loc="lower left", bbox_to_anchor=(0.0, 0.02))
    ax.set_title("Every architecture lands under our floor, and the intron count is not a "
                 "choice\n"
                 + "both ends scored on " + F["matrix"] + "; exons capped at 320 nt so the "
                 "oligo stays under the " + str(cap) + " nt order limit",
                 loc="left", fontsize=12.5, pad=12)
    fig.subplots_adjust(bottom=0.2)
    save(fig, "03_architecture_cliff")


# ---------------------------------------------------------------- 4. the host gate
def fig_hosts():
    """An intron is not a neutral addition in a host that cannot splice it."""
    H = F["hosts"]
    names = list(H["hosts"])
    pretty = {
        "c_reinhardtii_nuclear": "C. reinhardtii  nuclear",
        "c_reinhardtii_chloroplast": "C. reinhardtii  chloroplast",
        "n_tabacum_chloroplast": "N. tabacum  chloroplast",
        "e_coli": "E. coli",
        "s_cerevisiae": "S. cerevisiae",
        "a_thaliana_nuclear": "A. thaliana  nuclear",
    }
    cols = ["Has a spliceosome", "Evidence rbcS2i1 works here"]
    xs = [0.0, 1.45]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))

    for r, host in enumerate(names):
        v = H["hosts"][host]
        for c, ok in enumerate((v["spliceosome"], v["evidence"])):
            ax.add_patch(Rectangle((xs[c], -r - 0.42), 1.2, 0.84,
                                   facecolor=GREEN_L if ok else PINK_L,
                                   edgecolor="white", lw=2, zorder=2))
            ax.annotate("yes" if ok else "no", xy=(xs[c] + 0.6, -r), ha="center", va="center",
                        fontsize=10.5, color=GREEN_D if ok else PINK_D,
                        fontweight="bold", zorder=3)
        # Every one of these is a species name, E. coli included.
        ax.annotate(pretty[host], xy=(-0.12, -r), ha="right", va="center", fontsize=10.5,
                    style="italic")

    for c, label in enumerate(cols):
        ax.annotate(label, xy=(xs[c] + 0.6, 0.72), ha="center", fontsize=10, color=MUTED)

    ax.set_xlim(-2.5, 2.75)
    ax.set_ylim(-len(names) + 0.1, 1.5)
    ax.axis("off")
    ax.set_title("An intron is only defined for three of our six hosts, and only evidenced "
                 "for one", loc="left", fontsize=12.5, x=-0.0, y=1.06)
    ax.annotate("In a host with no spliceosome the intron is not ignored. rbcS2i1 is "
                + str(H["intron_nt"]) + " nt, so it shifts the frame by "
                + str(H["frame_shift"]) + ",\nand it carries a stop in all three frames. "
                "That host makes no protein, not less of it.",
                xy=(-2.45, -len(names) + 0.55), fontsize=9.5, color=PINK_D, va="top")
    save(fig, "04_host_gate")


if __name__ == "__main__":
    print("figures from " + str(DATA) + ":")
    fig_gc_band()
    fig_waterfall()
    fig_cliff()
    fig_hosts()
