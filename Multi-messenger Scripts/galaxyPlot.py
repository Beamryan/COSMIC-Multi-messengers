import pathlib, sys, configparser, pandas as pd, matplotlib.pyplot as plt

# Desired telescope limits to add plot lines
MAGNITUDE_LIMITS = [
        (20.7, 'Gaia 10% Precision ≈ 21',':', 'black'),
        (24.5, 'Vera Rubin Single Exposure ≈ 24.5',':', 'black'),
        (27.9, 'HST F606W ≈ 28 ',':', 'black'),
        (31.0, 'JWST NIRCam ≈ 31',':', 'black')]

# Checks if data and parameter files exist and are of the right type
def checkFiles(dataFile, paramFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if not paramFile.exists():
        sys.exit(f"Error! Parameter file not found: {paramFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

    if paramFile.suffix.lower() != ".ini":
        sys.exit("Error! Parameter file must end with .ini")

# Reads in axis parameters to be plotted from .ini file, stores as tuple type [x,y]
def read_params(paramFile):
    cfg = configparser.ConfigParser()
    cfg.read(paramFile)

    if "axes" not in cfg:
        sys.exit("Parameter file needs an [axes] section")

    axesSection = cfg["axes"]

    try:
        x_list = [s.strip() for s in axesSection["x"].split(",")]
        y_list = [s.strip() for s in axesSection["y"].split(",")]
        wdtype_list = [s.strip() for s in axesSection["wdtype"].split(",")]
    except KeyError as e:
        sys.exit(f"Section [axes] needs x, y, and wdtype entries. Missing: {e}")

    if not (len(x_list) == len(y_list) == len(wdtype_list)):
        sys.exit("x, y, and wdtype must have the same number of entries")

    return list(zip(x_list, y_list, wdtype_list))

# Checks for columns required to make requested plots 
def ensure_axes(df, axes):
    required = {col for pair in axes for col in pair}
    missing = list(required - set(df.columns)) # Adds missing columns to list if cant be found in dataframe
    if missing:
        sys.exit(f"Missing column(s): {', '.join(missing)}") 

# Makes the plots using the dataframe and corresponding axis titles.
def scatter(df, xAxis, yAxis, wdTypeCol, out_prefix):

    wd_types = pd.unique(df[wdTypeCol])
    wd_types = [wd for wd in wd_types if pd.notna(wd)]

    cmap = plt.colormaps.get_cmap("tab10")
    wd2color = {wd: cmap(i) for i, wd in enumerate(sorted(wd_types))}

    for gx_name, sub in df.groupby("gxComponent", sort=False):

        fig, ax = plt.subplots()

        for wd in wd_types:
            mask = (sub[wdTypeCol] == wd)

            if mask.any():
                ax.scatter(
                    sub.loc[mask, xAxis],
                    sub.loc[mask, yAxis],
                    s=9,
                    alpha=0.7,
                    color=wd2color[wd],
                    label=wd
                )

        ax.set_xlabel(xAxis, fontsize=12)
        ax.set_ylabel(yAxis, fontsize=12)
        ax.set_title(f"{gx_name} — {yAxis} vs {xAxis}", fontsize=14)
        ax.grid(ls="--", lw=0.4, alpha=0.5)
        ax.legend(title=wdTypeCol, markerscale=1.4, fontsize=12, loc="lower left")

        ax.set_ylim(10, 50)
        ax.set_xlim(0, 20)

        if "mag" in yAxis.lower():
            for y, lab, style, col in MAGNITUDE_LIMITS:
                ax.axhline(y, ls=style, lw=1, color=col)
                ax.text(
                    ax.get_xlim()[1],
                    y + 0.1,
                    lab,
                    va="bottom",
                    ha="right",
                    color=col,
                    fontsize=12
                )

            ax.invert_yaxis()

        out_png = f"{out_prefix}_{gx_name}_{yAxis}_{wdTypeCol}.png"
        fig.tight_layout()
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

        print(f"Saved {out_png}")

# Main Program for plotting graphs of given star parameters
def main():

    # Checks for the desired script inputs
    if len(sys.argv) != 3:
        sys.exit("Wrong inputs! Desired: python galaxyPlot.py input.csv Plots.ini")

    # Assigns file variables from inputs
    dataFile = pathlib.Path(sys.argv[1])
    paramFile = pathlib.Path(sys.argv[2])

    # Checks that files exist and are of appropriate type
    checkFiles(dataFile, paramFile)

    # Reads parameter file for axes
    print(f"Reading from {paramFile} ...")
    axes = read_params(paramFile)
    print("Complete!")
    
    # Reads data file to form dataframe
    print(f"Reading from {dataFile} ...")
    df = pd.read_csv(dataFile)
    print("Complete!")

    # Ensures columns exist to make desired plots with given axes
    ensure_axes(df, axes)

    # For each axes create the scatter plot
    print("Creating plots!")
    for xAxis, yAxis, wdTypeCol in axes:
        scatter(df, xAxis, yAxis, wdTypeCol, f"{xAxis}_{yAxis}")  
    print("All Done!")

if __name__ == "__main__":
    main()