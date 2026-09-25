import pathlib, sys, pandas as pd, matplotlib.pyplot as plt, numpy as np

# Checks if data file exists and is of the right type
def checkFiles(dataFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

def scatter(dataFile):

    mags1 = dataFile["mag1"].to_numpy() # apparent magnitude of primary
    mags2 = dataFile["mag2"].to_numpy() # apparent magnitude of secondary

    bins = np.arange(5, 40, 1)

    fig, ax = plt.subplots(figsize=(7,4))
    ax.hist(mags1, bins=bins, alpha=0.7, label="Primary", color="steelblue")
    ax.hist(mags2, bins=bins, alpha=0.5, label="Secondary", color="goldenrod")
    ax.set_xlabel("Apparent magnitude")
    ax.set_ylabel("Number of binaries")
    ax.set_title("COSMIC binaries: magnitude distribution")

    limit = 15.0
    n1 = (mags1 < limit).sum()
    n2 = (mags2 < limit).sum()
    msg = f"{n1:,} primaries  and  {n2:,} secondaries brighter than {limit}"

    ax.text(0.98, 0.98, msg,       
        transform=ax.transAxes,   
        ha="right", va="top",
        color="black", fontsize=9)

    ax.invert_xaxis()            # optional: bright-to-faint from left to right
    ax.set_yscale('log')
    ax.grid(ls="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig("AppMagHistogram.png", dpi=300)
    plt.close(fig) 
    print(f"Saved AppMagHistogram.png!")

# Main Program for plotting graphs of given star parameters
def main():

    # Checks for the desired script inputs
    if len(sys.argv) != 2:
        sys.exit("Wrong inputs! Desired: python galaxyHistogram.py input.csv")

    # Assigns file variables from inputs
    dataFile = pathlib.Path(sys.argv[1])

    # Checks that files exist and are of appropriate type
    checkFiles(dataFile)
    
    # Reads data file to form dataframe
    print(f"Reading from {dataFile} ...")
    df = pd.read_csv(dataFile)
    print("Complete!")

    # For each axes create the scatter plot
    print("Creating plots!")
    scatter(df)  
    print("All Done!")

if __name__ == "__main__":
    main()