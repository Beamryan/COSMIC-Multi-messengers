import sys, pathlib, numpy as np, pandas as pd, math
from astropy.io import fits 
from astropy import units as u
from scipy.interpolate import RegularGridInterpolator

# Checks if data and density map files exist and are of the right type
def checkFiles(dataFile, densityFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if not densityFile.exists():
        sys.exit(f"Error! Density file not found: {densityFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

    if densityFile.suffix.lower() != ".fits":
        sys.exit("Error! Density file must end with .fits")

# Checks to assure that required columns exist inside data file
def columnCheck(bcm):

    # Required columns for later calculations
    requiredCol = ["mass_1(Msun)", "mass_2(Msun)", "lum_1(Lsun)", "lum_2(Lsun)", 
                   "D(pc)", "Porb(s)", "xGx(pc)", "yGx(pc)", "zGx(pc)"]

    # Iterates through data file columns to see if each from required list are not found
    missingCol = [c for c in requiredCol if c not in bcm.columns]
    if missingCol:
        sys.exit(f"Error! Missing columns in {bcm}: {', '.join(missingCol)}")

# Stores x,y,z position values of compact objects for later calculations
def calibrateInputFile(bcm):

    # Suns position in Galactocentric frame
    SUN_X_GC = 7.86 * u.kpc             
    SUN_Z_GC = 0.027 * u.kpc         

    # Extracts each position vector column and converts to numpy array in parsec units
    x_gc = bcm["xGx(pc)"].to_numpy() / 1000 * u.kpc
    y_gc = bcm["yGx(pc)"].to_numpy() / 1000 * u.kpc
    z_gc = bcm["zGx(pc)"].to_numpy() / 1000 * u.kpc

    # Shifts position of stars from Galactocentric to HelioCentric frame
    x_hc = x_gc + SUN_X_GC
    y_hc = y_gc
    z_hc = z_gc - SUN_Z_GC

    bcm["xHx(kpc)"] = x_hc.value
    bcm["yHx(kpc)"] = y_hc.value
    bcm["zHx(kpc)"] = z_hc.value

    return [x_hc,y_hc,z_hc]

# Calibrates density file with given HI and H2 number density to mass density in g*cm^-3 
# returns mass density galactic grid in kpc
#
# File layout:
# No.    Name      Ver    Type      Cards   Dimensions   Format
#  0  PRIMARY       1 PrimaryHDU       4   ()      
#  1  MEAN OF HI NUMBER DENSITY   1 ImageHDU   25   (1251, 1251, 125)  float32 (x,y,z)   
#  2  MEAN OF H2 NUMBER DENSITY   1 ImageHDU   25   (1251, 1251, 125)  float32 (x,y,z) 
#  3  STANDARD DEVIATION OF HI NUMBER DENSITY   1 ImageHDU   25   (1251, 1251, 125)  float32 (x,y,z)    
#  4  STANDARD DEVIATION OF H2 NUMBER DENSITY   1 ImageHDU   25   (1251, 1251, 125)  float32 (x,y,z)  
def calibrateDensityFile(hdul, hdu_index=1, include_h2=True):

    # Finds number density of HI
    n_hi   = hdul[hdu_index].data.astype("f4") # float32 type
    header = hdul[hdu_index].header

    # Finds number density of H2 if desired and adds them to total
    if include_h2:
        n_h2 = hdul[hdu_index + 1].data.astype("f4") # float32 type
        n_tot = n_hi + 2.0 * n_h2
    else:
        n_tot = n_hi

    protonMass  = 1.6735575e-24 # Proton mass in grams
    massDensity  = n_tot * 1.4 * protonMass # Number density converted to mass density  (include He via 1.4 factor)

    num_x, num_y, num_z = massDensity.shape # Number of values on each axis
    pixelIndices = [np.arange(num_x), np.arange(num_y), np.arange(num_z)] # Array of pixel indices 

    # Converts pixel indices into absolute physical galactic coordinates in kpc
    # CRPIX#: References the pixel number inside given cube axis
    # CDELT#: References the size in kpc for each pixel
    # CRVAL#: Physical coordinate at reference pixel
    # Offests from reference pixel using CRPIX, then converts that to distance using CDELT,
    # then shifts absolute position using CRVAL
    x_grid = (pixelIndices[0] - (header["CRPIX3"] - 1)) * header["CDELT3"] + header["CRVAL3"]
    y_grid = (pixelIndices[1] - (header["CRPIX2"] - 1)) * header["CDELT2"] + header["CRVAL2"]
    z_grid = (pixelIndices[2] - (header["CRPIX1"] - 1)) * header["CDELT1"] + header["CRVAL1"]

    # Checks to make sure grid is oriented appropriately for data input
    assert massDensity.shape == (len(x_grid), len(y_grid), len(z_grid))
    print("x_grid:", x_grid.min(), x_grid.max())
    print("y_grid:", y_grid.min(), y_grid.max())
    print("z_grid:", z_grid.min(), z_grid.max())

    massDensityGrid = RegularGridInterpolator((x_grid, y_grid, z_grid), massDensity, bounds_error=False, fill_value=0.0)
    return massDensityGrid

# Calculates white dwarf type using mass
def calculateWdType(bcm):

    mass1 = bcm["mass_1(Msun)"] # defines mass variables from dataframe
    mass2 = bcm["mass_2(Msun)"]

    # creates dataframe column value for white dwarf types depending on mass
    bcm["wdType1"] = np.where((mass1 > 0.1) & (mass1 < 0.45),  "He",
                    np.where((mass1 >= 0.45) & (mass1 < 1.05), "C/O",
                    np.where((mass1 >= 1.05) & (mass1 < 1.55), "O/Ne", "unknown")))
    
    bcm["wdType2"] = np.where((mass2 > 0.1) & (mass2 < 0.45),  "He",
                    np.where((mass2 >= 0.45) & (mass2 < 1.05), "C/O",
                    np.where((mass2 >= 1.05) & (mass2 < 1.55), "O/Ne", "unknown")))

# Calculates gravitational wave amplitude
def calculateAmpGW(bcm):

    c = 299792458 # Speed of light (m * s^-1)
    gravCon = 6.6743 * 10**(-11) # Gravitational constant (m^3 * kg^-1 * s^-2)
    massSun = 1.9891 * 10**30 # Mass of the sun (kg)
    pc_to_m = (u.pc.to(u.m)) # Parsec to meters ratio

    mass1_mSun = bcm["mass_1(Msun)"] * massSun # Converts masses given in solar mass units to kilograms
    mass2_mSun = bcm["mass_2(Msun)"] * massSun

    orbitalPeriod = bcm["Porb(s)"].to_numpy()
    freq = 2 * (1/orbitalPeriod) # Frequency of orbiting binaries is TWICE the normal frequency (two bodies)
    bcm["fgw(Hz)"] = freq

    Dist = bcm["D(pc)"] * pc_to_m # Converts distance given in parsec units to meters

    # Calculates the chirpMass for grav. wave amp. later
    chirpMass = (((mass1_mSun * mass2_mSun)**(3/5))/((mass1_mSun + mass2_mSun)**(1/5)))
    
    # Calculates the gravitational wave amplitude
    bcm["gw_Amp"] = (4 * ((gravCon * chirpMass)**(5/3)) * ((math.pi * freq)**(2/3))) / ((c**4) * Dist)


# pos_xyz : (x,y,z) array in kpc  (Galactocentric)
# interpolatedMassDensity: Mass Density Grid
# opacity: Mass absorption coefficient 
# n_steps: Number of sample points along path
# returns optical depth  (dimensionless)
def calcOpticalDepth(pos, massDensityGrid, opacity=180, n_steps=500):

    kpc_to_cm = (u.kpc.to(u.cm)) # KiloParsec to cm ratio
    d_kpc = np.linalg.norm(pos) # Straight line distance in pc
    unitVector = pos / d_kpc # Unit vector toward star

    spaced_kpc = np.linspace(0.0, d_kpc, n_steps) # Linearly spaced sample points across given distance
    spaced_cm  = spaced_kpc * kpc_to_cm 
    pos_xyz_kpc = (unitVector[:, None] * spaced_kpc[None, :]).T # Array [n,3] of equally spaced way points in 3D (kpc) from sun to star

    # Integral of k * rho * ds here ds = x
    # Finds mass density at each point along path from observer to star
    opticalDepth = np.trapezoid(opacity * massDensityGrid(pos_xyz_kpc), x=spaced_cm)  
    return opticalDepth

# Calculates the apparent magnitude of star using calculated luminosity, distance, and dust extinction
def calcMagnitude(lum, dist_pc, extinction):
    
    Msun = 4.83 # Absolute magnitude of sun
    return Msun - 2.5*np.log10(lum) + 5.0*np.log10(dist_pc/10) + extinction # Apparent magnitude calculation

# Saves stars that are still visible after dust extinction calculation to csv file
def saveVisibleStars(bcm, limit):

    # Decides on visible stars
    visible = [(bcm["mag1"] < limit), (bcm["mag2"] < limit)]

    # Creates data frame of desired visible star parameters
    visible = pd.DataFrame({
    "#bin_num":  bcm.loc[(visible[0] | visible[1]), "#bin_num"],
    "gxComponent":  bcm.loc[(visible[0] | visible[1]), "gxComponent"],
    "tphys":  bcm.loc[(visible[0] | visible[1]), "tphys"],
    "xHx(kpc)":  bcm.loc[(visible[0] | visible[1]), "xHx(kpc)"],
    "yHx(kpc)":  bcm.loc[(visible[0] | visible[1]), "yHx(kpc)"],
    "zHx(kpc)":  bcm.loc[(visible[0] | visible[1]), "zHx(kpc)"],
    "D(pc)":  bcm.loc[(visible[0] | visible[1]), "D(pc)"],
    "dist(kpc)":  bcm.loc[(visible[0] | visible[1]), "dist(kpc)"],
    "mass_1(Msun)":  bcm.loc[visible[0], "mass_1(Msun)"],
    "mass_2(Msun)":  bcm.loc[visible[1], "mass_2(Msun)"],
    "lum_1(Lsun)":  bcm.loc[visible[0], "lum_1(Lsun)"],
    "lum_2(Lsun)":  bcm.loc[visible[1], "lum_2(Lsun)"],
    "mag1":  bcm.loc[visible[0], "mag1"],
    "mag2":  bcm.loc[visible[1], "mag2"],
    "wdType1":  bcm.loc[visible[0], "wdType1"],
    "wdType2":  bcm.loc[visible[1], "wdType2"],
    #"SNR":  bcm.loc[(visible[0] | visible[1]), "SNR"],
    "gw_Amp":  bcm.loc[(visible[0] | visible[1]), "gw_Amp"],
    "extinction":  bcm.loc[(visible[0] | visible[1]), "extinction"]})

    visible.to_csv(f"visible_stars{limit}.csv", index=False) # Saves data file

# Main Program for calculating apparent magnitude including dust extinction adjustment
def main():

    # Checks for correct number of script arguments, needs data file and density file
    if len(sys.argv) != 3:
        sys.exit("Wrong inputs! Desired: python magCalc.py data.csv densityFit.fits")

    # Finds and assigns file variables from script argument paths
    dataFile = pathlib.Path(sys.argv[1])
    densityFile = pathlib.Path(sys.argv[2])

    # Checks that data and density file exist and have correct type
    checkFiles(dataFile, densityFile)

    # Reads in dataFrame from data file
    print(f"Reading from {dataFile} ...")
    bcm = pd.read_csv(dataFile)
    print("Complete!")

    # Opens density file and assigns to hdul parameter 
    print(f"Reading from {densityFile} ...")
    hdul = fits.open(densityFile)
    print("Complete!")

    # Checks that required columns for future math exist
    columnCheck(bcm)

    # Stores x,y,z position values of compact objects for later calculations
    print(f"Calibrating coordinates from {dataFile} ...")
    [x,y,z] = calibrateInputFile(bcm)
    print("Complete!")

    # Creates mass density grid in kpc from provided .fits density file
    print(f"Calibrating coordinates and density values from {densityFile} ...")
    massDensityGrid = calibrateDensityFile(hdul)
    print("Complete!")

    # Calculating white dwarf types
    print(f"Calculating white dwarf types ...")
    calculateWdType(bcm)  
    print("Complete!")

    print("Calculating gravitational wave amplitude ...")
    calculateAmpGW(bcm) # Calculates gravitational wave amplitude  
    print("Complete!")

    # Calculates optical depth of straight line path from center of galaxy to star
    print("Calculating extinction using star coordinates and density grid ...")
    opticalDepth = np.array([calcOpticalDepth(pos, massDensityGrid) for pos in np.column_stack((x.value, y.value, z.value))])
    extinction = 1.086 * opticalDepth # Extinction linearly related to optical depth
    print("Complete!")
    
    bcm["extinction"] = extinction # Adds extinction column to datafile table

    # Recalculates or adds apparent magnitude column to datafile table using extinction correction
    print("Calculating apparent magnitude using extinction values ...")
    bcm["mag1"] = calcMagnitude(bcm["lum_1(Lsun)"].to_numpy(), bcm["D(pc)"].to_numpy(), extinction)   
    bcm["mag2"] = calcMagnitude(bcm["lum_2(Lsun)"].to_numpy(), bcm["D(pc)"].to_numpy(), extinction)
    print("Complete!")    

    bcm.to_csv(dataFile, index=False) # Saves modified data file

    saveVisibleStars(bcm, 50) # Saves visible stars
    saveVisibleStars(bcm, 40) # Saves visible stars
    saveVisibleStars(bcm, 20) # Saves visible stars

    print(f"Done! Wrote {dataFile} with magnitudes and saved visible star file!")

if __name__ == "__main__":
    main()
