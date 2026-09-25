import pandas as pd
import numpy as np
import pathlib, sys
import sympy as sy
from astropy import units as u

MSUN_TO_KG = 1.989 * (10**30)
KPC_TO_M = 3.0857 * (10**19)
MILHZ_TO_HZ = 1e-3

GravConstSI = 6.6743 * (10**-11)
cSI = 299792458 
GravConstMcKpcmHz = GravConstSI * ((1/KPC_TO_M)**3) * (MSUN_TO_KG) * ((1/MILHZ_TO_HZ)**2)
cKpcmHz = cSI * (1/KPC_TO_M) * (1/MILHZ_TO_HZ)

# Checks if data and parameter files exist and are of the right type
def checkFiles(dataFile):

    if not dataFile.exists():
        sys.exit(f"Error! Input file not found: {dataFile}")

    if dataFile.suffix.lower() != ".csv":
        sys.exit("Error! Input file must end with .csv")

# Checks to assure that required columns exist inside data file
def columnCheck(dataFile):

    # Required columns for later calculations
    requiredCol = ["xHx(kpc)", "yHx(kpc)", "zHx(kpc)", "mass_1(Msun)", "mass_2(Msun)", "D(pc)", "fgw(Hz)"]

    # Iterates through data file columns to see if each from required list are not found
    missingCol = [c for c in requiredCol if c not in dataFile.columns]
    if missingCol:
        sys.exit(f"Error! Missing columns in {dataFile}: {', '.join(missingCol)}")

# Calculates chirp mass using primary and secondary star
def calcChirpMass(dataFile):

    mass1 = dataFile["mass_1(Msun)"]
    mass2 = dataFile["mass_2(Msun)"]

    # Calculates chirp masses in solar units
    chirpMasses = ((mass1 * mass2) ** (3/5)) / ((mass1 + mass2) ** (1/5))

    return chirpMasses

# Calculates gravitational wave amplitudes
def calcAmplitudes(dataFile, chirpMasses):

    pc_to_m = (u.pc.to(u.m)) # Parsec to meters ratio

    freq = dataFile["fgw(Hz)"] # Frequency of orbiting binaries is TWICE the normal frequency (two bodies)

    Dist = dataFile["D(pc)"] * pc_to_m # Converts distance given in parsec units to meters
    
    # Calculates the gravitational wave amplitude
    amplitudes = (4 * ((GravConstMcKpcmHz * chirpMasses)**(5/3)) * ((np.pi * (freq/MILHZ_TO_HZ))**(2/3))) / ((cKpcmHz**4) * (Dist/KPC_TO_M))

    return amplitudes

# Creates random inclination angles for binary rotation
def makeRanInclinAngles(dataFile):
    
    N = len(dataFile["D(pc)"])
    cos_iota = np.random.uniform(-1.0, 1.0, N) # Uniform isotropic distribution
    iotas = np.arccos(cos_iota)

    return iotas

# Calculates random initial phase angles for gravitational wave (in freq atm *note commented out line for time func)
def calcPhaseAngles(dataFile):

    phiInitial = np.random.uniform(0.0, 2.0*np.pi, len(dataFile["D(pc)"]))

    #phaseAngles = 2 * np.pi * dataFile["fgw(Hz)"] * dataFile["tphys"] + phiInitial
    
    phaseAngles = phiInitial

    return phaseAngles

# Calculates random psi angles for rotation of gravitational wave
def calcPsiAngles(dataFile):

    psiInitial = np.random.uniform(0.0, 2.0*np.pi, len(dataFile["D(pc)"]))

    return psiInitial

# Calculates polar angles depending on heliocentric position
def calcPolarAngles(dataFile):

    x = dataFile["xHx(kpc)"]
    y = dataFile["yHx(kpc)"]
    z = dataFile["zHx(kpc)"]

    radius = np.sqrt((x ** 2) + (y ** 2) + (z ** 2))

    polar = np.arccos(z/radius)

    return polar

# Calculates azimuthal angles depending on heliocentric position
def calcAzimuthAngles(dataFile):

    x = dataFile["xHx(kpc)"]
    y = dataFile["yHx(kpc)"]

    azimuthal = np.arctan2(y,x)

    return azimuthal

# Calculates h-plus wave in frequency domain
def calcPlusWaveInFreq(amplitudes, inclinationAngles, phaseAngles):
    
    waves = amplitudes * 0.5 * (1 + np.cos(inclinationAngles)**2) * np.exp(1j * phaseAngles)

    return waves

# Calculates h-cross wave in frequency domain
def calcCrossWaveInFreq(amplitudes, inclinationAngles, phaseAngles):

    waves = 1j * amplitudes * np.cos(inclinationAngles) * np.exp(1j * phaseAngles)

    return waves

# Calculates f-plus detector response function in frequency domain
def calcfPlus(polarAngles, azimuthalAngles, psiAngles):
    
    fPlus = np.sqrt(3) / 2 * (0.5 * (1.0 + np.cos(polarAngles)**2) * np.cos(2.0 * azimuthalAngles) * np.cos(2.0 * psiAngles) \
        - np.cos(polarAngles) * np.sin(2.0 * azimuthalAngles) * np.sin(2.0 * psiAngles))
    
    return fPlus

# Calculates f-cross detector response function in frequency domain
def calcfCross(polarAngles, azimuthalAngles, psiAngles):

    fCross = np.sqrt(3) / 2 * (0.5 * (1.0 + np.cos(polarAngles)**2) * np.cos(2.0 * azimuthalAngles) * np.sin(2.0 * psiAngles) \
        + np.cos(polarAngles) * np.sin(2.0 * azimuthalAngles) * np.cos(2.0 * psiAngles))
    
    return fCross

# Simple Cornish-Robson style LISA noise estimate
def simple_noise_psd(freqs):

    L_arm = 2.5 * 10 ** 9 # Length of lisa arms in meters

    freqs = np.asarray(freqs)

    f_safe = np.maximum(freqs, 1e-12)
    f_m = cSI / (2.0 * np.pi * L_arm)

    P_oms = (1.5e-11) ** 2 * (1.0 + (2e-3 / f_safe) ** 4)
    P_acc = ((3e-15) ** 2 * (1.0 + (0.4e-3 / f_safe) ** 2) * (1.0 + (f_safe / 8e-3) ** 4))

    Sn = ((10.0 / (3.0 * (L_arm**2))) * (P_oms + 2 * (1 + (np.cos(f_safe / f_m) ** 2)) * (P_acc / ((2 * np.pi * f_safe) ** 4))) * (1 + (((6 / 10) * ((f_safe / f_m) ** 2)))))

    if Sn.ndim == 0:

        return float(Sn)

    return Sn

# Calculates partial derivatives for fisher matrix
def calcPartials():

    Mc, d, f, incl, polar, azim, phase, psi = sy.symbols('m d f i p a z y')

    # Entire h equation for h = f_plus * h_plus + f_cross * h_cross
    h = np.sqrt(3) / 2 * (0.5 * (1.0 + sy.cos(polar) **2) * sy.cos(2.0 * azim) * sy.cos(2.0 * psi) - sy.cos(polar) * sy.sin(2.0 * azim) * sy.sin(2.0 * psi)) \
    * ((4 * ((GravConstMcKpcmHz * Mc)**(5/3)) * ((sy.pi * f)**(2/3))) / ((cKpcmHz**4) * d)) * 0.5 * (1 + sy.cos(incl)**2) * sy.exp(1j * phase) \
    + np.sqrt(3) / 2 * (0.5 * (1.0 + sy.cos(polar) **2) * sy.cos(2.0 * azim) * sy.sin(2.0 * psi) + sy.cos(polar) * sy.sin(2.0 * azim) * sy.cos(2.0 * psi)) \
    * ((4 * ((GravConstMcKpcmHz * Mc)**(5/3)) * ((sy.pi * f)**(2/3))) / ((cKpcmHz**4) * d)) * 1j * sy.cos(incl) * sy.exp(1j * phase)

    dh_dM = sy.diff(h, Mc)
    print(dh_dM) # Print chirp mass partial

    dh_dd = sy.diff(h, d)
    print(dh_dd) # Print distance partial
    
    dh_df = sy.diff(h, f)
    print(dh_df) # Print frequency partial

    dh_dincl = sy.diff(h, incl)
    print(dh_dincl) # Print inclination angle partial

    dh_dpolar = sy.diff(h, polar)
    print(dh_dpolar) # Print polar angle partial

    dh_dazim = sy.diff(h, azim)
    print(dh_dazim) # Print azimuthal angle partial

    dh_dphase = sy.diff(h, phase)
    print(dh_dphase) # Print phase angle partial

    dh_dpsi = sy.diff(h, psi)
    print(dh_dpsi) # Print psi angle partial

    symbolic_derivs = {
        "dh_dM": sy.diff(h, Mc),
        "dh_dd": sy.diff(h, d),
        "dh_df": sy.diff(h, f),
        "dh_dincl": sy.diff(h, incl),
        "dh_dpolar": sy.diff(h, polar),
        "dh_dazim": sy.diff(h, azim),
        "dh_dphase": sy.diff(h, phase),
        "dh_dpsi": sy.diff(h, psi)
    }

    # Turn symbolic derivatives into fast numerical functions
    deriv_funcs = {
        name: sy.lambdify(
            (Mc, d, f, incl, polar, azim, phase, psi),
            expr,
            modules="numpy"
        )
        for name, expr in symbolic_derivs.items()
    }

    return deriv_funcs

def saveFisherResults(F, covariance, correlation, sigma, paramNames, outFile):

    fisherDF = pd.DataFrame(F,index=paramNames,columns=paramNames)
    covarianceDF = pd.DataFrame(covariance,index=paramNames,columns=paramNames)
    correlationDF = pd.DataFrame(correlation,index=paramNames,columns=paramNames)
    uncertaintyDF = pd.DataFrame({"parameter": paramNames,"sigma": sigma})

    summaryDF = pd.DataFrame({
        "quantity": [
            "condition_number",
            "num_parameters",
            "used_pseudoinverse"
        ],"value": [np.linalg.cond(F),len(paramNames),"yes"]})

    with pd.ExcelWriter(outFile, engine="openpyxl") as writer:
        summaryDF.to_excel(writer, sheet_name="Summary", index=False)
        uncertaintyDF.to_excel(writer, sheet_name="Uncertainties", index=False)
        fisherDF.to_excel(writer, sheet_name="Fisher Matrix")
        covarianceDF.to_excel(writer, sheet_name="Covariance Matrix")
        correlationDF.to_excel(writer, sheet_name="Correlation Matrix")

    print(f"Saved Fisher results to: {outFile}")

# Main Program for plotting graphs of given star parameters
def main():

    # Checks for the desired script inputs
    if len(sys.argv) != 2:
        sys.exit("Wrong inputs! Desired: python fisherMatrix.py input.csv")

    # Assigns file variables from inputs
    csv = pathlib.Path(sys.argv[1])
    dataFile = pd.read_csv(csv)

    # Checks that files exist and are of appropriate type
    checkFiles(csv)

    # Checks file has correct column headers
    columnCheck(dataFile)

    # Runs fisher matrix code and outputs to text file
    chirpMasses = calcChirpMass(dataFile)

    # Makes random inclination angles
    inclinationAngles = makeRanInclinAngles(dataFile)

    # Calculates gravitational wave amplitudes
    amplitudes = calcAmplitudes(dataFile, chirpMasses)

    # Calculates polar angles based off heliocentric positions
    polarAngles = calcPolarAngles(dataFile)

    # Calculates azimuthal angles based off heliocentric positions
    azimuthalAngles = calcAzimuthAngles(dataFile)

    # Makes random phase angles
    phaseAngles = calcPhaseAngles(dataFile)

    # Makes random psi angles
    psiAngles = calcPsiAngles(dataFile)

    # Calculates h_plus and h_cross waves in frequency domain
    hPlusWaves = calcPlusWaveInFreq(amplitudes, inclinationAngles, phaseAngles)
    hCrossWaves = calcCrossWaveInFreq(amplitudes, inclinationAngles, phaseAngles)

    # Calculates f_plus and f_cross detector response in frequency domain
    fPlus = calcfPlus(polarAngles, azimuthalAngles, psiAngles)
    fCross = calcfCross(polarAngles, azimuthalAngles, psiAngles)

    h = fPlus * hPlusWaves + fCross * hCrossWaves # Calculates h (gw strain)
    print(h)

    Sn = simple_noise_psd(dataFile["fgw(Hz)"]) # Calculates noise function

    derivs = calcPartials() # Calculates partial derivatives

    # Converts sympy values back to float array
    Mc_vals = np.asarray(chirpMasses, dtype=float)
    d_vals = np.asarray(dataFile["D(pc)"], dtype=float) / 1000
    f_vals = np.asarray(dataFile["fgw(Hz)"], dtype=float) / MILHZ_TO_HZ
    incl_vals = np.asarray(inclinationAngles, dtype=float)
    polar_vals = np.asarray(polarAngles, dtype=float)
    azim_vals = np.asarray(azimuthalAngles, dtype=float)
    phase_vals = np.asarray(phaseAngles, dtype=float)
    psi_vals = np.asarray(psiAngles, dtype=float)
    Sn_vals = np.asarray(Sn, dtype=float)

    deriv_arrays = {}

    for name, func in derivs.items():
        vals = func(
            Mc_vals,
            d_vals,
            f_vals,
            incl_vals,
            polar_vals,
            azim_vals,
            phase_vals,
            psi_vals)
        
        deriv_arrays[name] = np.asarray(vals, dtype=np.complex128)
        print(name, deriv_arrays[name].shape, deriv_arrays[name].dtype)

    # Frequency step size is 1 / observation time (in sec)
    Tobs =  4 * 365 * 24 * 60 * 60

    paramNames = [
    "dh_dM",
    "dh_dd",
    "dh_df",
    "dh_dincl",
    "dh_dpolar",
    "dh_dazim",
    "dh_dphase",
    "dh_dpsi"]

    numParams = len(paramNames)
    F = np.zeros((numParams, numParams))

    # Fisher matrix loop through derivatives
    for i, pi in enumerate(paramNames):
        for j, pj in enumerate(paramNames):

            dh_i = deriv_arrays[pi]
            dh_j = deriv_arrays[pj]

            integrand = np.conj(dh_i) * dh_j / Sn_vals

            F[i, j] = 4.0 * np.real(np.sum(integrand) * Tobs)

    print("Fisher matrix condition number:")
    print(np.linalg.cond(F))

    covariance = np.linalg.inv(F)
    sigma = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(sigma, sigma)

    saveFisherResults(
    F=F,
    covariance=covariance,
    correlation=correlation,
    sigma=sigma,
    paramNames=paramNames,
    outFile="fisher_matrix_results.xlsx")

    # P-inv correlation code
    F_sym = 0.5 * (F + F.T)
    covariance_pinv = np.linalg.pinv(F_sym)
    sigma_pinv = np.sqrt(np.diag(covariance_pinv))
    correlation_pinv = covariance_pinv / np.outer(sigma_pinv, sigma_pinv)

    saveFisherResults(
    F=F_sym,
    covariance=covariance_pinv,
    correlation=correlation_pinv,
    sigma=sigma_pinv,
    paramNames=paramNames,
    outFile="fisher_matrix_pinv_results.xlsx")


if __name__ == "__main__":
    main()