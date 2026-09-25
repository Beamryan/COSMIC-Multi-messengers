import sys, pathlib, numpy as np, pandas as pd
import matplotlib.pyplot as plt

# Checks if data file exists and is of the right type
def checkFile(dataFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

def main():

    # Checks for correct number of script arguments, needs data file 
    if len(sys.argv) != 2:
        sys.exit("Wrong inputs! Desired: python 3dPlot.py data.csv")

    # Finds and assigns file variables from script argument paths
    dataFile = pathlib.Path(sys.argv[1])

    # Checks that data file exists and has correct type
    checkFile(dataFile)

    # Reads in dataFrame from data file
    print(f"Reading from {dataFile} ...")
    df = pd.read_csv(dataFile)
    print("Complete!")

    x = df["xHx(kpc)"].to_numpy()
    y = df["yHx(kpc)"].to_numpy()
    z = df["zHx(kpc)"].to_numpy()
    A = df["extinction"].to_numpy()

    # Normalise A_V → colour map (Turbo)
    norm = (A - A.min()) / (np.ptp(A) + 1e-9)
    cmap = plt.cm.turbo
    col  = cmap(norm)

    # Build the figure
    fig = plt.figure()
    ax  = fig.add_subplot(projection='3d')

    xs = np.zeros_like(x)
    ys = np.zeros_like(y)
    zs = np.zeros_like(z)

    # Line quiver from Sun (0,0,0) to star (x,y,z)
    ax.quiver(xs, ys, zs, x, y, z, arrow_length_ratio = 0.01, colors=col, linewidths=1, alpha=0.4, zorder=0)

    # Labels, limits, view angle
    ax.set_xlabel("X [kpc]"); ax.set_ylabel("Y [kpc]"); ax.set_zlabel("Z [kpc]")
    ax.set_title("White-dwarf binaries – extinction-coloured rays")

    ax.set_xlim(-15, 35)   # show ±20 kpc in X
    ax.set_ylim(-30, 20)   # same for Y
    ax.set_zlim(-6, 8) 

    # colour-bar
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array(A)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.7, pad=0.05)
    cbar.set_label("$A_V$ [mag]")
    fig.tight_layout()

    ax.view_init(elev=20, azim=60)
    fig.savefig(f"elevatedPlot_{dataFile.stem}.png", dpi=300)

    ax.view_init(elev=90, azim=90)
    ax.axis('off')
    ax.grid(False)
    fig.savefig(f"topdownPlot_{dataFile.stem}.png", dpi=300, transparent=True)

    ax.set_yticks([])
    ax.view_init(elev=0, azim=90)
    fig.savefig(f"sidewayPlot_{dataFile.stem}.png", dpi=300)

    plt.close(fig) 


if __name__ == "__main__":
    main()