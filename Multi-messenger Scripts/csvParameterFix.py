import sys, pathlib, numpy as np, pandas as pd

# Checks if data file exists and is of the right type
def checkFiles(dataFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")


def csvConvert(dataFile):

    df = pd.read_csv(dataFile)
    df['D(pc)'] = df['dist(kpc)'] * 1000

    df['Porb(s)'] = df['porb_final(day)'] * 86400

    df['xGx(pc)'] = df['xGx(kpc)'] * 1000
    df['yGx(pc)'] = df['yGx(kpc)'] * 1000
    df['zGx(pc)'] = df['zGx(kpc)'] * 1000

    df.to_csv(dataFile, index=False)


# Main Program
def main():

    # Checks for correct number of script arguments, needs data file
    if len(sys.argv) != 2:
        sys.exit("Wrong inputs! Desired: python csvParameterFix.py data.csv")

    # Finds and assigns file variables from script argument paths
    dataFile = pathlib.Path(sys.argv[1])

    # Checks that data file exist has correct type
    checkFiles(dataFile)

    #
    print(f"Working on {dataFile} ...")
    csvConvert(dataFile)
    print("Complete!")

if __name__ == "__main__":
    main()
