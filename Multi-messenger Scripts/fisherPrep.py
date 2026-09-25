import pandas as pd
import pathlib, sys

# Checks if data and parameter files exist and are of the right type
def checkFiles(dataFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

# Checks to assure that required columns exist inside data file
def columnCheck(dataFile):

    # Required columns for later calculations
    requiredCol = ["xGx(kpc)", "yGx(kpc)", "zGx(kpc)", "dist(kpc)", "porb_final(day)"]

    # Iterates through data file columns to see if each from required list are not found
    missingCol = [c for c in requiredCol if c not in dataFile.columns]
    if missingCol:
        sys.exit(f"Error! Missing columns in {dataFile}: {', '.join(missingCol)}")

# Prepares csv with appropriate column headers for fisher code 
def appendColumns(df, csv):

    df['D(pc)'] = df['dist(kpc)'] * 1000

    df['Porb(s)'] = df['porb_final(day)'] * 86400
    df['fgw(Hz)'] = 2 * (1 / df['Porb(s)'])

    df['xHx(kpc)'] = df['xGx(kpc)'] + 7.86
    df['yHx(kpc)'] = df['yGx(kpc)'] 
    df['zHx(kpc)'] = df['zGx(kpc)'] - 0.027

    df.to_csv(csv, index=False)

# Main Program for plotting graphs of given star parameters
def main():

    # Checks for the desired script inputs
    if len(sys.argv) != 2:
        sys.exit("Wrong inputs! Desired: python fisherPrep.py input.csv")

    # Assigns file variables from inputs
    csv = pathlib.Path(sys.argv[1])
    dataFile = pd.read_csv(csv)

    # Checks that files exist and are of appropriate type
    checkFiles(csv)

    # Checks file has correct column headers
    columnCheck(dataFile)

    # Prepares csv with appropriate column headers for fisher code
    appendColumns(dataFile, csv)

if __name__ == "__main__":
    main()