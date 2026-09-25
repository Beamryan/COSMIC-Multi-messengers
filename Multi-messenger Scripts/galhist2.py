import sys, pathlib
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
BINS   = np.arange(10, 45, 5)          # 5‑mag bins: [5‑10), [10‑15), …, [50‑55)
LIMITS = [16.0, 20.7, 24.5, 27.0, 31.0]  # survey limits to annotate
SECTIONS = ["ThinDisk", "ThickDisk", "Bulge"]     # gxComponent names in your CSV

# ----------------------------------------------------------------------
def safe_read(path):
    path = pathlib.Path(path)
    if not path.exists() or path.suffix.lower() != ".csv":
        sys.exit("Input file must be an existing .csv")
    return pd.read_csv(path)

# ----------------------------------------------------------------------
def one_histogram(df, label):
    """Make one histogram for a given gxComponent."""
    mags = df["mag2"].to_numpy()
    counts, edges = np.histogram(mags, bins=BINS)

    fig, ax = plt.subplots(figsize=(7,4))
    ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge",
           color="bisque", alpha=0.8, edgecolor="slategray")
    ax.set_xlabel("Apparent magnitude")
    ax.set_ylabel("Number of stars")
    ax.set_title(f"{label.capitalize()}: Magnitude distribution")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.grid(ls="--", alpha=0.4)

    # annotate survey limits
    for lim in LIMITS:
        n_visible = (mags < lim).sum()
        ax.axvline(lim, ls=":", lw=1, color="red")
        ax.text(lim-0.2, 200,
                f"N<{lim:g} = {n_visible:,}",
                color="black", ha="right", va="center", rotation=90, fontsize=14)

    fig.tight_layout()
    outfile = f"hist_{label}.png"
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    print(f"   saved  {outfile}")

    # print counts table to console
    print(f"\n{label.upper()} 5‑mag bin counts")
    for l, r, c in zip(edges[:-1], edges[1:], counts):
        print(f" {l:>2.0f}–{r:<2.0f}: {c:>7,}")
    print()

# ----------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python galaxyHistogram.py visible_stars.csv")

    df = safe_read(sys.argv[1])
    missing = set(SECTIONS) - set(df["gxComponent"].unique())
    if missing:
        sys.exit(f"gxComponent entries missing from file: {', '.join(missing)}")

    for sec in SECTIONS:
        sub = df[df["gxComponent"] == sec]
        print(f"\n=== {sec.capitalize()} component === (N = {len(sub):,})")
        one_histogram(sub, sec)

if __name__ == "__main__":
    main()