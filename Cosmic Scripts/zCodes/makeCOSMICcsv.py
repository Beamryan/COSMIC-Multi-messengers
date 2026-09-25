# CODE     : makeCOSMICcsv.py
# AUTHOR   : Shane L. Larson (s.larson@northwestern.edu)
# DATE     : 06 July 2021
# MODIFIED : 17 Mar 2022
# --------------------------------------------------------------------
# This converts COSMIC HDF5 files to CSV files (suitable to be read
#   and used by other analysis codes (eg. those written in C)
#
# REQUIRES:
#    -- galaxy files from COSMIC (in directory with script)
#    -- file  "gxModel.txt" with filenames and filepaths (same used to make galaxy)
# --------------------------------------------------------------------
#
# This script allows you to decide what columns to include in the output
# CSV file; that means you need to know what the columns in the HDF5 file
# are!  You can probe the columns in an hdf5 file in python via:
#
#   import pandas
#   fileIN = pandas.read_hdf('filename')
#   fileIN.columns
#
# Additionally, you need to know the KEYs in the file. You can also probe
# the keys in an hdf5 file in python via:
#
#   import h5py
#   fileIN = h5py.File('basename_He_He_ThinDisk_metallicity_0.017.h5', 'r')
#   list(fileIN.keys())
#
# --------------------------------------------------------------------
#
# VERSION NOTES
#  -- v 0.0 [2020] convert_COSMICtoCSV_v3.py from Julie Malewicz
#  -- v 1.0 [2021] makeCOSMICcsv.py, in tandem with my COSMIC galaxy building scripts
#  -- (v2) Winter 2022: updated to output spin columns for Minli
#
#  -- (v3) 17 Mar 2022: updated to run with auto script, so points up to
#      needed directories.
#
# --------------------------------------------------------------------

''' `makeCOSMICcsv`

# What does this do ?
#####################

This script allows you to convert full Milky Way COSMIC files from .hdf format to .csv.
The .hdf format is primarily used for the transfer of large datasets as it allows you to save on space, 
but cannot be read as is by the user.

MW COSMIC files are sampled from COSMIC fixed populations in order to go from a statistical representative
to an actual astrophysical realization of a given population.

The script is accompanied by comments so as to make it easier for the user to modify it to fit their purpose.

You can run the script by entering
    python makeCOSMICcsv.py

If you are a python user, you can also simply open the files in your python script with the following:
import pandas as pd
myfile = pd.read_hdf(filename, key='galaxy')

Below, some info about the files themselves:

# file columns
##################################################################
# The current standard output of file columns. For components, the
# same quantity exists for both the PRIMARY (_1) and SECONDARY (_2)
# in the same units (writing both out suppressed here for space)
# A full accounting in the COSMIC documentation here:
#   https://cosmic-popsynth.github.io/docs/v3.3.0/output_info/index.html
bin_num		   Unique binary index	
tphys          Evolution time		            [Myr]
mass_1		   Primary mass		                [Msun]
kstar_1		   Type of primary star	            see kstar section
kstar_2		   Type of secondary star	        see kstar section
sep      	   Separation	                    [Rsun]
porb           Orbital period                   [days]
ecc            Eccentricity	                    [ ]
RRLO_1         Primary radius divided by Roche lobe radius
aj_1           Effective age of the primary     [Myr]
tms_1          Primary main sequence lifetime   [Myr]
massc_1        Primary core mass                [Msun]
rad_1          Primary radius                   [Rsun]
mass0_1        Previous evolutionary stage primary mass    [Msun]
lum_1          Primary luminosity               [Lsun]
teff_1         Primary effective temperature    [K]
radc_1         Primary core radius              [Rsun]
menv_1         Mass of envelope of primary      [Msun]
renv_1         Radius of envelope of primary    [Rsun]
omega_spin_1   Angular speed of primary         [rad/yr]
B_1            Neutron star magnetic field      [G]
bacc_1         dMsun in accretion (pulsars)     [see Eq 7 in COSMIC paper]
tacc_1         Accretion duration               [Myr] (used for B field decay)
epoch_1
bhspin_1       Black hole sping magintude       [Kerr parameter]
t_birth        birth time if gx born at t = 0   [Myr]
t_evol_GW
sep_final
ecc_final
porb_final
RL_1_final
f_gw_peak      Peak GW frequency	            [Hz]
xGx, yGx, zGx  Galacto-centric coord.           [kpc]
dist           Distance to the Sun              [kpc]

# kstar:
10: Helium White Dwarf          [He]
11: Carbon Oxygen White Dwarf   [CO]
12: Oxygen Neon White Dwarf     [ONe]
13: Neutron Star                [NS]
14: Black Hole                  [BH]

[Physical coordinates]
xGx, yGx and zGx are assigned following the McMillan+2002 model
'''
import pandas as pd

# -----------------------------
# sll imports
# -----------------------------

from time import time, ctime


# -----------------------------
# USER INPUT
# -----------------------------

# We offer the user the possibility to only convert NS+BH or DWD populations
# pop = input("Population to convert: all, NS+BH, DWD\n --> ")
#
# Display an error message if the input isn't recognized
# while pop not in ['NS+BH', 'all', 'DWD']:
#     pop = input("Please respond: all, NS+BH, or DWD \n --> ")

#pop = 'DWD'            # ### change at runtime; could also be  'NS+BH' or  'all'

fileformat = 'full'    # ### change at runtime; could also be 'kstar'  'comp'  or  'sep'


# Depending on the user's choice, the list of kstars to convert can be one of 3:
# if pop == 'DWD':
#     kstar_list = ['He_He','CO_He','CO_CO', 'ONe_He_ONe']
# elif pop == 'NS+BH':
#     kstar_list = ['NS_NS','BH_NS','BH_BH']
# elif pop == 'all':
#     kstar_list = ['He_He','CO_He','CO_CO',
#                   'ONe_He_ONe','NS_He_ONe','BH_He_ONe',
#                   'NS_NS','BH_NS','BH_BH']

# print(" ")
# print("Desired Output: \n",
#       "    - full == full galaxy (all components, all types) \n",
#       "    - comp == separate structure components (thinDisk, thickDisk, bulge) \n",
#       "    - kstar == separate files for each binary kstar type (up to 9 in total) \n",
#       "    - sep == separate files for each kstar AND each gx component (up to 27 total)")
# fileformat = input("Enter any other input for a more detailed explanation of those options \n --> ")
#
# while fileformat not in ['full', 'kstar', 'comp', 'sep']:
#     print(" full:  1 file only \n",
#           "       gx_dat_full.csv \n",
#           "comp:  3 files for a full gx\n",
#           "       gx_dat_Bulge.csv \n",
#           "       gx_dat_ThinDisk.csv \n",
#           "       gx_dat_ThickDisk.csv \n",
#           "kstar: 9 files for a full gx (4 for DWD and 3 for NS+BH) \n",
#           "       gx_dat_He_He.csv \n",
#           "       [...]\n",
#           "       gx_dat_BH_BH.csv \n",
#           "sep:   27 files for a full gx (12 for DWD and 9 for NS+BH) \n",
#           "       gx_dat_He_He_ThinDisk.csv \n",
#           "       gx_dat_He_He_ThickDisk.csv \n",
#           "       gx_dat_He_He_Bulge.csv \n",
#           "       [...]\n",
#           "       gx_dat_BH_BH_ThinDisk.csv \n",
#           "       gx_dat_BH_BH_ThickDisk.csv \n",
#           "       gx_dat_BH_BH_Bulge.csv \n")
#
#     fileformat = input("Enter Option: full, kstar, comp, sep? \n --> ")




print('==================================\n')
print('COSMIC hdf5 to CSV conversions\n')
print(" ")
print(" Converting .h5 galaxy files into .csv",
      " and counting number of binary systems in each file. \n",
      "Please be aware this might take up some time (> 2hrs).")


##############
# REAL STUFF #
##############

# Now the real stuff begins:
# We make 4 different loops, for each file format chosen by the user
# to make reading the code a bit easier maybe
#
#
# NB: If you have NS elements, then rejection should include RL_1_final and
#     RL_2_final -- these can always be computed from the stellar parameters
#     later if needed, using Eggleton+1983
#
# Here is the complete column list to edit for column rejection:
#
# reject_col = ['bin_num', 'tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep',
#              'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 'aj_1', 'aj_2', 'tms_1',
#              'tms_2', 'massc_1', 'massc_2', 'rad_1', 'rad_2', 'mass0_1', 'mass0_2',
#              'lum_1', 'lum_2', 'teff_1', 'teff_2', 'radc_1', 'radc_2', 'menv_1',
#              'menv_2', 'renv_1', 'renv_2', 'omega_spin_1', 'omega_spin_2', 'B_1',
#              'B_2', 'bacc_1', 'bacc_2', 'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
#              'bhspin_1', 'bhspin_2', 't_birth', 't_evol_GW', 'sep_final',
#              'ecc_final', 'porb_final', 'RL_1_final', 'RL_2_final', 'f_gw_peak',
#              'xGx', 'yGx', 'zGx', 'dist']

# my usual columns to reject...
# reject_col = ['sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 'aj_1', 'aj_2', 'tms_1',
#              'tms_2', 'massc_1', 'massc_2', 'rad_1', 'rad_2', 'mass0_1', 'mass0_2',
#              'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
#              'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2',
#              'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
#              'bhspin_1', 'bhspin_2', 'RL_1_final', 'RL_2_final', 'f_gw_peak']

# this set keeps the spin of the components
# reject_col = ['sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 'aj_1', 'aj_2', 'tms_1',
#              'tms_2', 'massc_1', 'massc_2', 'rad_1', 'rad_2', 'mass0_1', 'mass0_2',
#              'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
#              'B_1', 'B_2', 'bacc_1', 'bacc_2',
#              'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
#              'bhspin_1', 'bhspin_2', 'RL_1_final', 'RL_2_final', 'f_gw_peak']

# this set keeps the spin of the components; also keeps evol_type for Patti
# reject_col = ['sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'aj_1', 'aj_2', 'tms_1',
#              'tms_2', 'massc_1', 'massc_2', 'rad_1', 'rad_2', 'mass0_1', 'mass0_2',
#              'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
#              'B_1', 'B_2', 'bacc_1', 'bacc_2',
#              'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
#              'bhspin_1', 'bhspin_2', 'RL_1_final', 'RL_2_final', 'f_gw_peak']

# this set keeps the spin of the components; also keeps evol_type for Patti, B_1 for Nycole
reject_col = ['sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'aj_1', 'aj_2', 'tms_1',
              'tms_2', 'massc_1', 'massc_2', 'mass0_1', 'mass0_2',
              'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
              'bacc_1', 'bacc_2',
              'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
              'bhspin_1', 'bhspin_2', 'f_gw_peak']


# ===========================================================
#  (sll) READ IN GALAXY METADATA FILE
#    Looks for the file "gxProcessInputs.txt" which has the
#    basefilename for this galaxy run.
# ===========================================================
#with open('../buildFiles/gxProcessInputs.txt','r') as gxProcess:
#    dataLine = gxProcess.readline()           # read in dataline
#    parseLine = dataLine.split("=")           # parses line around the = sign
#    basename = str(parseLine[-1]).strip()     # rstrip removes trailing \n
#    print("Galaxy Basename  = ", basename)    # echo to screen


# This is a more complicated version of doing the output that Julie's
# original script because I don't loop over the components, I instead have
# a separate section for EACH of the components, which allows me to build
# the filenames in a way consistent with how I do them when producing the
# galaxy from COSMIC

# ===========================================================
#  (sll) READ IN GALAXY METADATA FILE
#    Looks for the file "cosmicRuntimeValues.txt" which has
#    basefilename for this galaxy run and populations used
# ===========================================================

startRun = time()  # get the starting time in seconds for elapsed time calculation
print('\n==================================\n')
print('Starting run at: ',ctime(time()))

# get the data for filename building from the galaxy model file

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




# Now go through the combinations of files

# When making fixed populations, we needed to keep the COMMON and RARE
# cases separate, but here we can merge the two, so make a new master
# list to loop over.

kstarNameLoop = kstarsName + kstarsRareName    # concatenate the two name lists

if fileformat == 'sep':
    for kstarNameJJ in kstarNameLoop:
        print('Now converting ', kstarNameJJ)
        
        for gx_componentII in gx_components:
            if gx_componentII == 'ThickDisk':
                metComp = metThick
            elif gx_componentII == 'Bulge':
                metComp = metBulge
            else:
                metComp = metThin

            # list of columns to delete to save space
            col_to_del = reject_col.copy()

            # Build the Filename and Load the File:
            datFileIN = '../gx_data/' + basename + '_' + kstarNameJJ + '_' + gx_componentII + '_metallicity_' + str(metComp) +'.h5'
            gx = pd.read_hdf(datFileIN, key='LISA_population')
        
            # we drop a column present only for WDs, radius/RL radius which can be computed later by user if needed
            # following Eggleton+1983 e.g.
#            if kstarNameJJ in ['He_He','CO_He','CO_CO','ONe_He_ONe']:
#                col_to_del.extend(['RL_1_final', 'RL_2_final'])
#            elif kstarNameJJ in ['NS_He_ONe','BH_He_ONe']:
#                col_to_del.append('RL_2_final')
            gx = gx.drop(columns=col_to_del)
            
            # sll - 16 Jun 2023 -- added this so create a column with what gxComponent star came from
            gx['gxComponent'] = gx_componentII
            gx.to_csv('../csv/'+ basename + '_' + kstarNameJJ + '_' + gx_componentII + '.csv', index=False)
        
            # some accounting: we add a `#` to the header, making it readable as such by C codes
            gx.rename(columns = {'bin_num':'#bin_num'}, inplace = True)
            gx.to_csv('../csv/'+ basename + '_' + kstarJJ + '_' + gx_componentII + '.csv', index=False)
            print('  {:.2E} {} in the {}'.format(len(gx), kstarJJ, gx_componentII))

elif fileformat == 'comp':
    for gx_componentII in gx_components:
        if gx_componentII == 'ThickDisk':
            metModel = metThick
        elif gx_componentII == 'Bulge':
            metModel = metBulge
        else:
            metModel = metThin
        
        print('Now converting ', gx_componentII)
        # comp_gx = pd.DataFrame()
        comp_gx = []   # make galaxy a list, per Peggy Guo solution to errors; 18 Dec 2023
        for kstarNameJJ in kstarNameLoop:
            col_to_del = reject_col.copy()
            datFileIN = '../gx_data/' + basename + '_' + kstarNameJJ + '_' + gx_componentII + '_metallicity_' + str(metComp) +'.h5'
            gx = pd.read_hdf(datFileIN, key='LISA_population')

            if kstarNameJJ in ['He_He','CO_He','CO_CO','ONe_He_ONe']:
                col_to_del.extend(['RL_1_final', 'RL_2_final'])
            elif kstarNameJJ in ['NS_He_ONe','BH_He_ONe']:
                col_to_del.append('RL_2_final')   
            gx = gx.drop(columns=col_to_del)

            # sll - 16 Jun 2023 -- added this so create a column with what gxComponent star came from
            gx['gxComponent'] = gx_componentII
            print('  {:.2E} {} in the {}'.format(len(gx), kstarNameJJ, gx_componentII))
            comp_gx.append(gx)

        comp_gx = pd.concat(comp_gx) # Convert list to dataframe, per Peggy Guo, 18 Dec 2023
        comp_gx.rename(columns = {'bin_num':'#bin_num'}, inplace = True)
        comp_gx.to_csv('../csv/'+ basename + '_' + gx_componentII + '.csv', index=False)
            
elif fileformat == 'kstar':
    for kstarNameJJ in kstarNameLoop:
        print('Now converting ', kstarNameJJ)
        # kstar_gx = pd.DataFrame()
        kstar_gx = []       # make galaxy a list, per Peggy Guo solution to errors; 18 Dec 2023
        for gx_componentII in gx_components:
            if gx_componentII == 'ThickDisk':
                metModel = metThick
            elif gx_componentII == 'Bulge':
                metModel = metBulge
            else:
                metModel = metThin
            
            col_to_del = reject_col.copy()
            
            datFileIN = '../gx_data/' + basename + '_' + kstarNameJJ + '_' + gx_componentII + '_metallicity_' + str(metComp) +'.h5'
            gx = pd.read_hdf(datFileIN, key='LISA_population')
            
            if kstarNameJJ in ['He_He','CO_He','CO_CO','ONe_He_ONe']:
                col_to_del.extend(['RL_1_final', 'RL_2_final'])
            elif kstarNameJJ in ['NS_He_ONe','BH_He_ONe']:
                col_to_del.append('RL_2_final')   
            gx = gx.drop(columns=col_to_del)
            # sll - 16 Jun 2023 -- added this so create a column with what gxComponent star came from
            gx['gxComponent'] = model
            print('  {:.2E} {} in the {}'.format(len(gx), kstarNameJJ, gx_componentII))
            kstar_gx.append(gx)
        kstar_gx = pd.concat(kstar_gx)      # Convert list to dataframe, per Peggy Guo, 18 Dec 2023
        kstar_gx.rename(columns = {'bin_num':'#bin_num'}, inplace = True)
        kstar_gx.to_csv('../csv/'+ basename + '_' +kstar+'.csv', index=False)
        
        
elif fileformat == 'full':
    # full_gx = pd.DataFrame()
    full_gx = []   # make galaxy a list, per Peggy Guo solution to errors; 18 Dec 2023
    for kstarNameJJ in kstarNameLoop:
        print('Now converting', kstarNameJJ)
        for gx_componentII in gx_components:

            if gx_componentII=='ThickDisk':
                metComp = metThick
            elif gx_componentII=='Bulge':
                metComp = metBulge
            else:
                metComp = metThin
            
            col_to_del = reject_col.copy()
            datFileIN = '../gx_data/' + basename + '_' + kstarNameJJ + '_' + gx_componentII + '_metallicity_' + str(metComp) +'.h5'
            gx = pd.read_hdf(datFileIN, key='LISA_population')

            if kstarNameJJ in ['He_He','CO_He','CO_CO','ONe_He_ONe']:
                col_to_del.extend(['RL_1_final', 'RL_2_final'])
            elif kstarNameJJ in ['NS_He_ONe','BH_He_ONe']:
                col_to_del.append('RL_2_final')   
            gx = gx.drop(columns=col_to_del)

            # sll - 16 Jun 2023 -- added this so create a column with what gxComponent star came from
            gx['gxComponent'] = gx_componentII
            print('  {:.2E} {} in the {}'.format(len(gx), kstarNameJJ, gx_componentII))
            # print(type(full_gx))
            full_gx.append(gx)   # append gx to the list

    full_gx = pd.concat(full_gx)      # Convert list to dataframe, per Peggy Guo, 18 Dec 2023
    full_gx.rename(columns = {'bin_num':'#bin_num'}, inplace = True)
    full_gx.to_csv('../csv/'+ basename + '_' + 'full.csv', index=False)




# =============================
#  COMPLETE
# =============================

print('Galaxy CSV conversion complete.\n', flush=True)

endRun = time()                        # ending time in seconds
secondsElapsed = endRun - startRun     # total elapsed time in seconds

print('----------------------------------\n')
print('Ending run at: ',ctime(time()))
print('Elapsed time: ', int(secondsElapsed // 3600),' HR ',int(secondsElapsed % 3600 % 60),' MIN\n')
print('==================================\n')

