# CODE    :  makeGalaxy.py
# AUTHOR  :  Shane L. Larson (s.larson@northwestern.edu)
# MODIFIED:  17 Mar 2022 (v7)
# --------------------------------------------------------------------
# REQUIRES:
#    -- MW_maker_lw.py   (kmb, sll modded -- uses legwork, "_lw" name)
#    -- legwork       http://legwork.readthedocs.io
#    -- MC_samp.py   (kmb, from COSMIC 3.3)
#    -- fixed population from COSMIC
#    -- file  "cosmicRuntimeValues.txt" with filenames and filepaths
# --------------------------------------------------------------------
# VERSION NOTES
#  -- v 0.0 [2021]
#  -- first script, following interactive example at
#     https://cosmic-popsynth.github.io/docs/stable/install/index.html
#     for "Generate a population the COSMIC way"
#     with guidance from Julie Malewicz Jupyter Notebook
#     scale-to-MW-COSMIC3.ipynb
#
# MODIFICATIONS from SLL:
# (0) Deleted all plotting routines during dev (clarify debugging)
# (1) addition of new dat_path_2, dat_path_3 so THIN/BULGE/THICK all have
#     separate FIXED POPS (allows different seed, even if have same model
#     params. Makes better LISA populations without line sampling at high N
# (2) (v6) Made MC_samp as file to import, rather than from COSMIC
#     also updated physical constants to be consistent with MC_samp
#     CORRECTED ERROR IN PARSEC CONVERSION --> x10^16 meters (error was x10^19)
#     thought it doesn't look like it is used anywhere...
# (3) (v7) added some directory leads to a standard structure I use so I don't
#     have to sort and move files by hand to run codes farther down the pipeline
# (4) (v8) Commented out all the LISA bits, since I post-process LISA with a
#     C code. Will delete this all once confident it runs.  Also, GW_calcs and
#     MC_samp not used here; they are used in MW_makerS2.py; comment them out?
# (5) (gitCommit) 21 May 2025: debugged version, now runs throgh with WD,NS,BH
#     now under git version control, filenames simplifed.
#-----------------------------------------------------------------------------
#
# 10/31/2022
# Modifications: Natalie Gottschlich (nataliegottschlich2022@u.northwestern.edu)
#  -- Changed variable 'model' to 'gx_component' to avoid degeneracy with mc_samp
#  -- Likewise changed 'metModel' to 'metComp'



'''`makeGalaxy`
'''


# =============================
#  Python imports
# =============================

import time
import numpy as np
import pandas as pd
import MW_maker_lw
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import colors

# -----------------------------
# sll imports
# -----------------------------

from time import time, ctime

# =============================
#  Setup CONSTANTS
# =============================
# sll: these should all be  updated to higher precision; especially the
#      yr_sec, which in principle affects LISA data analysis; don't know
#      how it propagates through the COSMIC production though
#      Follow values I use in my C headers, documented provenance for values

Rsun_au = 1/215.032 # Radius of the sun in AU
day_yr = 1/365.25   # interesting -- what is this and is it consistent with yr_sec? don't see it used
yr_sec = 3.155e7    # seconds in a year  [Sidereal Year = 31558149.7635456 sec]
G = 6.67384e-11     # SI units
c = 299792458.0     # SI units: m/s
Msun = 1.9891e30    # SI units: kg
parsec = 3.0856775814913673e16  # SI units: meters


# =============================
#  USER DEFINED QUANTITIES
#  Update these before each run
# =============================

# useful color hex code/RGB/CMYK/HSL tool: https://htmlcolorcodes.com
#   #7fcdbb = turquoise
#   #225ea8 = med blue
#   #081d58 = dark blue
# Colors used in the plot
# pcolors = ['#7fcdbb',  '#225ea8', '#081d58', 'black']

# Path to the save folder (where the user wishes to save output data)
# ### change at runtime ###
dat_save_path = '../gx_data/'

# Define path to plots folder (where the image files will save)
# ### change at runtime ###
# plot_path = './plots/'


# ===========================================================
#  (sll) READ IN RUNTIME METADATA FILE
#    Looks for the file "cosmicRuntimeValues.txt" which has the
#    info and filenames to the FIXED POPULATIONS for each of the
#    three components of the galaxy. Test development code in
#    file  shane_gxRead.py
# ===========================================================
#    (sll) 23 Sep 2024: combined "gxInputs.txt" and "gxProcess.txt"
#       into a single file that needs edited, to streamline
# ===========================================================

startRun = time()  # get the starting time in seconds for elapsed time calculation
print('==================================\n')
print('Starting run at: ',ctime(time()))

with open('../buildFiles/cosmicRuntimeValues.txt','r') as gxInputs:
    # First, read in the basename for all files in this run
    dataLine = gxInputs.readline()  # read 1 line -- 1st comment line
    dataLine = gxInputs.readline()  # read 1 line -- 2nd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 3rd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 4th comment line
    dataLine = gxInputs.readline()  # read blankline
    dataLine = gxInputs.readline()           # read in dataline
    parseLine = dataLine.split("=")           # parses line around the = sign
    basename = str(parseLine[-1]).strip()     # rstrip removes trailing \n
    print("Galaxy Basename  = ", basename)    # echo to screen

    dataLine = gxInputs.readline()          # read blankline
    dataLine = gxInputs.readline()          # read blankline

    # Next, reading in the inputs for the galaxy model
    dataLine = gxInputs.readline()  # read 1 line -- 1st comment line
    dataLine = gxInputs.readline()  # read 1 line -- 2nd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 3rd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 4th comment line
    dataLine = gxInputs.readline()  # read 1 line -- 5th comment line

    # Star Formation Lines ----------------------------------
    dataLine = gxInputs.readline()          # read blankline
    dataLine = gxInputs.readline()          # read in dataline
    parseLine = dataLine.split("=")         # parses line around the = sign
    SF_start = float(parseLine[-1])         # store the number
    print("SF Start (Myr)   = ", SF_start)  # echo to screen
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    SF_duration = float(parseLine[-1])
    print("SF Duration (Myr)   = ", SF_duration)
   
    starFormation = 'SFstart_' + str(SF_start)  + '_SFduration_' +  str(SF_duration)     # construct the filename bit with Star Formation
    print("Star Formation file tag  = ", starFormation)
    
    dataLine = gxInputs.readline()          # read blankline

    # Metallicities ----------------------------------------
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    metThin = float(parseLine[-1])
    print("Thin Disk Metallicity   = ", metThin)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    metBulge = float(parseLine[-1])
    print("Bulge Metallicity       = ", metBulge)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    metThick = float(parseLine[-1])
    print("Thick Disk Metallicity  = ", metThick)
    
    dataLine = gxInputs.readline()          # read blankline

    # Pathways to Fixed Population Files ---------------------
    dataLine = gxInputs.readline()          # read in dataline
    parseLine = dataLine.split("=")         # parses line around the = sign
    dat_path_1 = str(parseLine[-1]).strip() # rstrip removes trailing \n
    print("thin disk PATH  = ", dat_path_1) # echo to screen
    
    dataLine = gxInputs.readline()          # read blankline
    parseLine = dataLine.split("=")
    dat_path_2 = str(parseLine[-1]).strip()
    print("bulge PATH      = ", dat_path_2)
    
    dataLine = gxInputs.readline()          # read blankline
    parseLine = dataLine.split("=")
    dat_path_3 = str(parseLine[-1]).strip()
    print("thick disk PATH = ", dat_path_3)
    
    dataLine = gxInputs.readline()          # read blankline

    # Galaxy Generation RANDOM SEED ---------------------
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    gxSeed = int(parseLine[-1])
    gxSeed = np.random.randint(1,10000000000)
    print("galaxy generation random seed  = ", gxSeed)
    print("fixed pop info retrieved\n", flush=True) # flush -- force stdout so I know its running

    dataLine = gxInputs.readline()          # read blankline
    dataLine = gxInputs.readline()          # read blankline

    # Next, reading in the inputs for the sub-populations
    dataLine = gxInputs.readline()  # read 1 line -- 1st comment line
    dataLine = gxInputs.readline()  # read 1 line -- 2nd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 3rd comment line
    dataLine = gxInputs.readline()  # read 1 line -- 4th comment line
    dataLine = gxInputs.readline()  # read 1 line -- 5th comment line
    
    dataLine = gxInputs.readline()          # read blankline

    # Common White Dwarf Combos ------------------------
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagHeHe = int(parseLine[-1])
    print("He_He Flag  = ", flagHeHe)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagCoHe = int(parseLine[-1])
    print("CO_He Flag  = ", flagCoHe)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagCoCo = int(parseLine[-1])
    print("CO_CO Flag  = ", flagCoCo)
    
    # Less Common White Dwarf Combos -------------------
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagONeWD = int(parseLine[-1])
    print("ONe_He_ONe Flag  = ", flagONeWD)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagONeNS = int(parseLine[-1])
    print("NS_He_ONe Flag  = ", flagONeNS)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagONeBH = int(parseLine[-1])
    print("BH_He_ONe Flag  = ", flagONeBH)

    # Heavy Remnant Combos ------------------------
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagNSNS = int(parseLine[-1])
    print("NS_NS Flag  = ", flagNSNS)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagBHNS = int(parseLine[-1])
    print("BH_NS Flag  = ", flagBHNS)
    
    dataLine = gxInputs.readline()
    parseLine = dataLine.split("=")
    flagBHBH = int(parseLine[-1])
    print("BH_BH Flag  = ", flagBHBH)
    
    dataLine = gxInputs.readline()          # read blankline
    dataLine = gxInputs.readline()          # read blankline

    # Next, reading in the LISA Processing Specifications
    # No read here; currently not used in python scripts
    # FUTURE PLAN: Make calling Legwork an option?



# ===========================================================
# The different subpopulations
# ===========================================================
gx_components = ['ThinDisk', 'Bulge', 'ThickDisk']
gxc_labels = ['Thin Disk', 'Bulge', 'Thick Disk'] # was gx_component_labels

kstars = []                                  # empty list to start building our populations list
kstarsName = []                              # empty list for population names

if (flagHeHe == 1):
   kstars.append([10,10])
   kstarsName.append('He_He')
if (flagCoHe == 1):
   kstars.append([11,10])
   kstarsName.append('CO_He')
if (flagCoCo == 1):
   kstars.append([11,11])
   kstarsName.append('CO_CO')
if (flagNSNS == 1):
   kstars.append([13,13])
   kstarsName.append('NS_NS')
if (flagBHNS == 1):
   kstars.append([14,13])
   kstarsName.append('BH_NS')
if (flagBHBH == 1):
   kstars.append([14,14])
   kstarsName.append('BH_BH')

kstarsRare = []
kstarsRareName = []

if (flagONeWD == 1):
   kstarsRare.append([12,10,12])
   kstarsRareName.append('ONe_He_ONe')
if (flagONeNS == 1):
   kstarsRare.append([13,10,12])
   kstarsRareName.append('NS_He_ONe')
if (flagONeBH == 1):
   kstarsRare.append([14,10,12])
   kstarsRareName.append('BH_He_ONe')


# ====================================================================
#  CREATE A MILKY WAY POPULATION
#
#  This is where the galaxy is built from the fixed populations. It
#  leans on the MW_makerS2 routines from kmb in the COSMIC distrubution
#  which were imported at the top.  Note the change to the first lines
#  that set the data paths, to accomodate different paths to all three
#  components.
#  Note the data is large
#     ~ 15 GB
# ====================================================================

# --------------------------------------------------------------------------
# First, this does the COMMON (He, CO) combinations ------------------------
# --------------------------------------------------------------------------
for kstarII, kstarNameJJ in zip(kstars, kstarsName):
    for gx_componentII in gx_components:
 
        print('\n COMMON kstarII = ',kstarII,'\n kstarnameJJ = ',kstarNameJJ,'\n gxcomponentII = ',gx_componentII,flush=True)
 
        if gx_componentII == 'ThickDisk':
            dat_path = dat_path_3   # thick disk
            metComp = metThick
        elif gx_componentII == 'Bulge':
            dat_path = dat_path_2   # bulge
            metComp = metBulge
        else:
            dat_path = dat_path_1   # thin disk
            metComp = metThin
        
        print(' metComp = ',metComp,'\n',flush=True)

        # Filename construction:
        #   dat_path      : read in from gxModel.txt, directory with fixed populations
        #   dat_          : COSMIC default leading filename for fixed population files
        #   kstar1_       : COSMIC default filename string indicator for kstar1 value
        #   str(kstar[0]) : make a string from the integer value of kstar1
        #   _kstar2_      : COSMIC default filename string indicator for kstar2 value
        #   str(kstar[1]) : make a string from the integer value of kstar2
        #   str(metModel) : make a string from the metallicty for current component (read from gxModel.txt)
        
        datFileIN = 'dat_' + 'kstar1_' + str(kstarII[0]) + '_kstar2_' + str(kstarII[-1]) + '_' +  starFormation + '_metallicity_' + str(metComp) + '.h5'
        datFileOUT = dat_save_path + basename + '_' + kstarNameJJ + '_' + gx_componentII + '_metallicity_' + str(metComp) +'.h5'
        # datLISApwr = dat_save_path+'dat_power_'+kstar_save+'_'+model+ '_metallicity_' + str(metModel) +'.h5'

        print('datPath IN = ', dat_path,flush=True)
        print('datFile IN = ', datFileIN,flush=True)
        print('datFile OUT = ', datFileOUT,flush=True)

        conv, m_sim, m_sim_tot = MW_maker_lw.read_dat(dat_path=dat_path, dat_file=datFileIN, kstar_types=kstarII)
        
        # ----------------------------------------------------------------------
        # sll 6 May 2025: CONV seems like the obvious place to trap the issue with kstar_1 not being
        #   a NS or BH when we requested that it be so.  insert and test that code here.
        # following this stackOverflow: https://stackoverflow.com/questions/64492382/
        # I should make this a callable function, so I don't have to repeat it in this code
        
        if (kstarII[0] == 13) or (kstarII[0] == 14):
            swap_at = [10,11,12]    # these are the non-NS values that contaminate kstar_1 == all WD types
        
            # make a mask for the kstar_1 column -- these are the entries in the data frame that are WD for kstar1
            swap_at_mask = conv["kstar_1"].isin(swap_at)
        
            # now swap all the columns for the kstar1 = WD entries that need swapped, kstar1 <--> kstar2
            conv.loc[swap_at_mask,["mass_1","mass_2"]] = conv.loc[swap_at_mask,["mass_2","mass_1"]].to_numpy()
            conv.loc[swap_at_mask,["kstar_1","kstar_2"]] = conv.loc[swap_at_mask,["kstar_2","kstar_1"]].to_numpy()
            conv.loc[swap_at_mask,["RRLO_1","RRLO_2"]] = conv.loc[swap_at_mask,["RRLO_2","RRLO_1"]].to_numpy()
            conv.loc[swap_at_mask,["aj_1","aj_2"]] = conv.loc[swap_at_mask,["aj_2","aj_1"]].to_numpy()
            conv.loc[swap_at_mask,["tms_1","tms_2"]] = conv.loc[swap_at_mask,["tms_2","tms_1"]].to_numpy()
            conv.loc[swap_at_mask,["massc_1","massc_2"]] = conv.loc[swap_at_mask,["massc_2","massc_1"]].to_numpy()
            conv.loc[swap_at_mask,["rad_1","rad_2"]] = conv.loc[swap_at_mask,["rad_2","rad_1"]].to_numpy()
            conv.loc[swap_at_mask,["mass0_1","mass0_2"]] = conv.loc[swap_at_mask,["mass0_2","mass0_1"]].to_numpy()
            conv.loc[swap_at_mask,["lum_1","lum_2"]] = conv.loc[swap_at_mask,["lum_2","lum_1"]].to_numpy()
            conv.loc[swap_at_mask,["teff_1","teff_2"]] = conv.loc[swap_at_mask,["teff_2","teff_1"]].to_numpy()
            conv.loc[swap_at_mask,["radc_1","radc_2"]] = conv.loc[swap_at_mask,["radc_2","radc_1"]].to_numpy()
            conv.loc[swap_at_mask,["menv_1","menv_2"]] = conv.loc[swap_at_mask,["menv_2","menv_1"]].to_numpy()
            conv.loc[swap_at_mask,["renv_1","renv_2"]] = conv.loc[swap_at_mask,["renv_2","renv_1"]].to_numpy()
            conv.loc[swap_at_mask,["omega_spin_1","omega_spin_2"]] = conv.loc[swap_at_mask,["omega_spin_2","omega_spin_1"]].to_numpy()
            conv.loc[swap_at_mask,["B_1","B_2"]] = conv.loc[swap_at_mask,["B_2","B_1"]].to_numpy()
            conv.loc[swap_at_mask,["bacc_1","bacc_2"]] = conv.loc[swap_at_mask,["bacc_2","bacc_1"]].to_numpy()
            conv.loc[swap_at_mask,["tacc_1","tacc_2"]] = conv.loc[swap_at_mask,["tacc_2","tacc_1"]].to_numpy()
            conv.loc[swap_at_mask,["epoch_1","epoch_2"]] = conv.loc[swap_at_mask,["epoch_2","epoch_1"]].to_numpy()
            conv.loc[swap_at_mask,["bhspin_1","bhspin_2"]] = conv.loc[swap_at_mask,["bhspin_2","bhspin_1"]].to_numpy()

        print('\nComponent: {}\nkstar: {}'.format(gx_componentII, kstarNameJJ))
        print('Number in converged population: {}'.format(len(conv)))
        print('Total mass required to make converged population: {}'.format(m_sim_tot), flush=True)
        
        # (sll) -- here I modded the line where they construct the filename to output the galaxy
        # Note you can include all stars in the call to MW Maker by specifying 'all stars' as in this commented line:
        # LISA_Gx = MW_makerS2.LISA_Galaxy(conv=conv, model=model, m_sim_tot=m_sim_tot, kstars='all stars')
        # LISA_Gx = MW_makerS2.LISA_Galaxy(conv=conv, model=model, m_sim_tot=m_sim_tot, kstars=kstar_save)
        LISA_Gx = MW_maker_lw.LISA_Galaxy(conv=conv, gx_component = gx_componentII, m_sim_tot=m_sim_tot, kstars=kstarII, dat_write=datFileOUT)
        
        # (sll) -- this is the old call, before MW_Maker wrote its own files
        # LISA_Gx.to_hdf(datFileOUT, key='galaxy')
        print('LISAgx done... ',flush=True)
        print('LISAgx to hdf done... ',flush=True)
        
        conv = []
        LISA_Gx = []
    
    endRun = time()                        # ending time in seconds
    secondsElapsed = endRun - startRun     # total elapsed time in seconds
    print('Time so Far: ', int(secondsElapsed // 3600),' HR ',int(secondsElapsed % 3600 % 60),' MIN',flush=True)
    print('{} done!\n\n'.format(kstarNameJJ), flush=True)


# --------------------------------------------------------------------------
# Now this does the RARE (ONe + Any; NS or BH) combinations ----------------
# --------------------------------------------------------------------------
# for kstar, kstar_save in zip(kstars_1, kstars_save_1):
for kstarII, kstarNameJJ in zip(kstarsRare, kstarsRareName):
    for gx_componentII in gx_components:
 
        print('\n RARE kstarII = ',kstarII,'\n kstarnameJJ = ',kstarNameJJ,'\n gxcomponentII = ',gx_componentII,flush=True)
 
        if gx_componentII == 'ThickDisk':
            dat_path = dat_path_3   # thick disk
            metComp = metThick
        elif gx_componentII == 'Bulge':
            dat_path = dat_path_2   # bulge
            metComp = metBulge
        else:
            dat_path = dat_path_1   # thin disk
            metComp = metThin
    
        print(' metComp = ',metComp,'\n',flush=True)

        # Filename construction:
        #  dat_path      : read in from gxModel.txt, directory with fixed populations
        #  dat_          : COSMIC default leading filename for fixed population files
        #  kstar1_       : COSMIC default filename string indicator for kstar1 value
        #  str(kstar[0]) : make a string from the integer value of kstar1
        #  _kstar2_      : COSMIC default filename string indicator for kstar2 value
        #  str(kstar[1]) : make a string from the integer value of kstar2
        #  str(metModel) : make a string from the metallicty for current component (read from gxModel.txt)
        
        datFileIN = 'dat_' + 'kstar1_' + str(kstarII[0]) + '_kstar2_' + str(kstarII[-2]) + '_' + str(kstarII[-1]) + '_' +  starFormation + '_metallicity_' + str(metComp) + '.h5'
        datFileOUT = dat_save_path + basename + '_' + kstarNameJJ +'_'+ gx_componentII + '_metallicity_' + str(metComp) +'.h5'
        # datLISApwr = dat_save_path+'dat_power_'+kstar_save+'_'+model+ '_metallicity_' + str(metModel) +'.h5'

        print('datPath IN = ', dat_path,flush=True)
        print('datFile IN = ', datFileIN,flush=True)
        print('datFile OUT = ', datFileOUT,flush=True)

        conv, m_sim, m_sim_tot = MW_maker_lw.read_dat(dat_path=dat_path, dat_file=datFileIN, kstar_types=kstarII)
        
        # ----------------------------------------------------------------------
        # sll 6 May 2025: CONV seems like the obvious place to trap the issue with kstar_1 not being
        #   a NS or BH when we requested that it be so.  insert and test that code here.
        # following this stackOverflow: https://stackoverflow.com/questions/64492382/
        # I should make this a callable function, so I don't have to repeat it in this code
        
        if (kstarII[0] == 13) or (kstarII[0] == 14):
            swap_at = [10,11,12]    # these are the non-NS values that contaminate kstar_1 == all WD types
        
            # make a mask for the kstar_1 column -- these are the entries in the data frame that are WD for kstar1
            swap_at_mask = conv["kstar_1"].isin(swap_at)
        
            # now swap all the columns for the kstar1 = WD entries that need swapped, kstar1 <--> kstar2
            conv.loc[swap_at_mask,["mass_1","mass_2"]] = conv.loc[swap_at_mask,["mass_2","mass_1"]].to_numpy()
            conv.loc[swap_at_mask,["kstar_1","kstar_2"]] = conv.loc[swap_at_mask,["kstar_2","kstar_1"]].to_numpy()
            conv.loc[swap_at_mask,["RRLO_1","RRLO_2"]] = conv.loc[swap_at_mask,["RRLO_2","RRLO_1"]].to_numpy()
            conv.loc[swap_at_mask,["aj_1","aj_2"]] = conv.loc[swap_at_mask,["aj_2","aj_1"]].to_numpy()
            conv.loc[swap_at_mask,["tms_1","tms_2"]] = conv.loc[swap_at_mask,["tms_2","tms_1"]].to_numpy()
            conv.loc[swap_at_mask,["massc_1","massc_2"]] = conv.loc[swap_at_mask,["massc_2","massc_1"]].to_numpy()
            conv.loc[swap_at_mask,["rad_1","rad_2"]] = conv.loc[swap_at_mask,["rad_2","rad_1"]].to_numpy()
            conv.loc[swap_at_mask,["mass0_1","mass0_2"]] = conv.loc[swap_at_mask,["mass0_2","mass0_1"]].to_numpy()
            conv.loc[swap_at_mask,["lum_1","lum_2"]] = conv.loc[swap_at_mask,["lum_2","lum_1"]].to_numpy()
            conv.loc[swap_at_mask,["teff_1","teff_2"]] = conv.loc[swap_at_mask,["teff_2","teff_1"]].to_numpy()
            conv.loc[swap_at_mask,["radc_1","radc_2"]] = conv.loc[swap_at_mask,["radc_2","radc_1"]].to_numpy()
            conv.loc[swap_at_mask,["menv_1","menv_2"]] = conv.loc[swap_at_mask,["menv_2","menv_1"]].to_numpy()
            conv.loc[swap_at_mask,["renv_1","renv_2"]] = conv.loc[swap_at_mask,["renv_2","renv_1"]].to_numpy()
            conv.loc[swap_at_mask,["omega_spin_1","omega_spin_2"]] = conv.loc[swap_at_mask,["omega_spin_2","omega_spin_1"]].to_numpy()
            conv.loc[swap_at_mask,["B_1","B_2"]] = conv.loc[swap_at_mask,["B_2","B_1"]].to_numpy()
            conv.loc[swap_at_mask,["bacc_1","bacc_2"]] = conv.loc[swap_at_mask,["bacc_2","bacc_1"]].to_numpy()
            conv.loc[swap_at_mask,["tacc_1","tacc_2"]] = conv.loc[swap_at_mask,["tacc_2","tacc_1"]].to_numpy()
            conv.loc[swap_at_mask,["epoch_1","epoch_2"]] = conv.loc[swap_at_mask,["epoch_2","epoch_1"]].to_numpy()
            conv.loc[swap_at_mask,["bhspin_1","bhspin_2"]] = conv.loc[swap_at_mask,["bhspin_2","bhspin_1"]].to_numpy()

        print('Component: {}\nkstar: {}'.format(gx_componentII, kstarNameJJ))
        print('Number in converged population: {}'.format(len(conv)))
        print('Total mass required to make converged population: {}'.format(m_sim_tot), flush=True)
        
        # (sll) -- here I modded the line where they construct the filename to output the galaxy
        # Note you can include all stars in the call to MW Maker by specifying 'all stars' as in this commented line:
        # LISA_Gx = MW_makerS2.LISA_Galaxy(conv=conv, model=model, m_sim_tot=m_sim_tot, kstars='all stars')
        # LISA_Gx = MW_makerS2.LISA_Galaxy(conv=conv, model=model, m_sim_tot=m_sim_tot, kstars=kstar_save)
        LISA_Gx = MW_maker_lw.LISA_Galaxy(conv=conv, gx_component = gx_componentII, m_sim_tot=m_sim_tot, kstars=kstarNameJJ, dat_write=datFileOUT)
        
        #(sll) -- this is the old call, before MW_Maker wrote its own files
        #LISA_Gx.to_hdf(datFileOUT, key='galaxy')
        print('LISAgx done... ',flush=True)
        print('LISAgx to hdf done... ',flush=True)
        
        conv = []
        LISA_Gx = []
        
    endRun = time()                        # ending time in seconds
    secondsElapsed = endRun - startRun     # total elapsed time in seconds
    print('Time so Far: ', int(secondsElapsed // 3600),' HR ',int(secondsElapsed % 3600 % 60),' MIN',flush=True)
    print('{} done!\n\n'.format(kstarNameJJ), flush=True)


# =============================
#  COMPLETE
# =============================

print('Galaxy complete.\n', flush=True)

endRun = time()                        # ending time in seconds
secondsElapsed = endRun - startRun     # total elapsed time in seconds

print('----------------------------------\n')
print('Ending run at: ',ctime(time()))
print('Elapsed time: ', int(secondsElapsed // 3600),' HR ',int(secondsElapsed % 3600 % 60),' MIN\n')
print('==================================\n')

