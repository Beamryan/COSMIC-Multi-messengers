#!/bin/bash

# ---------------------------------------------------------------
# -------- TUTORIAL NOTES (sll) ---------------------------------
# ---------------------------------------------------------------
# sll shell script template
# following tutorial: www.linux.com/training-tutorials/writing-simple-bash-script/
#
# bash Scripting Cheatsheet: https://devhints.io/bash
#
# Script needs to be executable, so make sure permissions are set correctly.
# To get rwx only by you:  % chmod 700 scriptname
#
# COMMAND LINE ARGUMENTS
# If use script name followed by arguments at the command line
# then each successive argument can be passed to commands in
# the script, proceeded by a $ and followed by the numerical
# order they appear in.  For example for script fooGalaxy :
#    % fooGalaxy MilkyWay.csv
# in the script then can use $1 anytime it needs MilkyWay.csv
#    % echo $1
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# SCRIPT   : makeFixedPops
# DATE     : 20 May 2021
# MODIFIED : 17 Mar 2022 (v2)
#
# This script takes a set of params.ini type files and
# makes successive fixed populations in COSMIC for
# the subpopulation files by STELLAR TYPE and by
# GALACTIC COMPONENT. These will be the inputs for
# making a FULL GALAXY.
#
# INPUTS: ---------------------
# Input files or command line arguments this script needs
#
# VARIABLES -------------------
# List of variables and what they are goes here for reference
#
# DATADIR = Directory where the fixed population files should be stored
# SECONDS = system var that counts time
# duration = current value of SECONDS
# durationCOHe = SECONDS when COHe calculation ends
# durationCOCO = SECONDS when COCO calculation ends
# durationONe = SECONDS when ONe+Any calculation ends
#
# ---------------------------------------------------------------

NOW=$(date)
echo "START TIME: $NOW"
echo "   Started at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo "------------------------------"

(( SECONDS = 0 ))          ## set time at start of script

# ---------------------------------------------------
# BULGE FIXED POPS
# ---------------------------------------------------

echo "FIXED POPULATIONS FOR BULGE...."
echo

# ---- He + He Systems (kstar1 = 10; kstar2 = 10)
echo "now working on He+He: 10 and 10..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 10 --final-kstar2 10 --inifile ../buildFiles/bulge_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME He+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + He Systems (kstar1 = 11; kstar2 = 10)
echo "now working on CO+He: 11 and 10..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 10 --inifile ../buildFiles/bulge_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + CO Systems (kstar1 = 11; kstar2 = 11)
echo "now working on CO+CO: 11 and 11..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 11 --inifile ../buildFiles/bulge_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+CO = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- ONe + Any Systems (kstar1 = 12; kstar2 = 10,11,12)
echo "now working on ONe+Any: 12 and 10 11 12..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 12 --final-kstar2 10 12 --inifile ../buildFiles/bulge_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME ONe+any = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + WD Systems (kstar1 = 13; kstar2 = 10,11,12)
echo "now working on NS+WD: 13 and 10,11,12..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 10 12 --inifile ../buildFiles/bulge_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + WD Systems (kstar1 = 14; kstar2 = 10,11,12)
# echo "now working on BH+WD: 14 and 10,11,12..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 10 12 --inifile ../buildFiles/bulge_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo

# ---- NS + NS Systems (kstar1 = 13; kstar2 = 13)
echo "now working on NS+NS: 13 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 13 --inifile ../buildFiles/bulge_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + BH Systems (kstar1 = 14; kstar2 = 13)
echo "now working on BH+NS: 14 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 13 --inifile ../buildFiles/bulge_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + BH Systems (kstar1 = 14; kstar2 = 14)
# echo "now working on BH+BH: 14 and 14..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 14 --inifile ../buildFiles/bulge_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

# (( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+BH = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
# echo


# --- note use of MV instead of mv -- how I have my unix set up.
#     for standard unix use:
mv ./dat_* ../fixed/fixed_bulge/.
mv ./log_* ../fixed/fixed_bulge/.

# MV ./dat_* ../fixed/fixed_bulge/.
# MV ./log_* ../fixed/fixed_bulge/.

echo "Bulge complete."
echo
echo "------------------------------"
echo

# ---------------------------------------------------
# THICK DISK FIXED POPS
# ---------------------------------------------------

echo "FIXED POPULATIONS FOR THICK DISK...."
echo

# ---- He + He Systems (kstar1 = 10; kstar2 = 10)
echo "now working on He+He: 10 and 10..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 10 --final-kstar2 10 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME He+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + He Systems (kstar1 = 11; kstar2 = 10)
echo "now working on CO+He: 11 and 10..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 10 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + CO Systems (kstar1 = 11; kstar2 = 11)
echo "now working on CO+CO: 11 and 11..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 11 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+CO = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- ONe + Any Systems (kstar1 = 12; kstar2 = 10,11,12)
echo "now working on ONe+Any: 12 and 10 11 12..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 12 --final-kstar2 10 12 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME ONe+any = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + WD Systems (kstar1 = 13; kstar2 = 10,11,12)
echo "now working on NS+WD: 13 and 10,11,12..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 10 12 --inifile ../buildFiles/thick_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + WD Systems (kstar1 = 14; kstar2 = 10,11,12)
# echo "now working on BH+WD: 14 and 10,11,12..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 10 12 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo

# ---- NS + NS Systems (kstar1 = 13; kstar2 = 13)
echo "now working on NS+NS: 13 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 13 --inifile ../buildFiles/thick_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + BH Systems (kstar1 = 14; kstar2 = 13)
echo "now working on BH+NS: 14 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 13 --inifile ../buildFiles/thick_Params.ini --Nstep 100000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + BH Systems (kstar1 = 14; kstar2 = 14)
# echo "now working on BH+BH: 14 and 14..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 14 --inifile ../buildFiles/thick_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

# (( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+BH = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
# echo


# --- note use of MV instead of mv -- how I have my unix set up.
#     for standard unix use:
mv ./dat_* ../fixed/fixed_thickDisk/.
mv ./log_* ../fixed/fixed_thickDisk/.

# MV ./dat_* ../fixed/fixed_thickDisk/.
# MV ./log_* ../fixed/fixed_thickDisk/.

echo "Thick Disk complete."
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# THIN DISK FIXED POPS
# ---------------------------------------------------

echo "FIXED POPULATIONS FOR THIN DISK...."
echo

# ---- He + He Systems (kstar1 = 10; kstar2 = 10)
echo "now working on He+He: 10 and 10..."
echo "   Started at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"

cosmic-pop --final-kstar1 10 --final-kstar2 10 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME He+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + He Systems (kstar1 = 11; kstar2 = 10)
echo "now working on CO+He: 11 and 10..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 10 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+He = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- CO + CO Systems (kstar1 = 11; kstar2 = 11)
echo "now working on CO+CO: 11 and 11..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 11 --final-kstar2 11 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME CO+CO = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo


# ---- ONe + Any Systems (kstar1 = 12; kstar2 = 10,11,12)
echo "now working on ONe+Any: 12 and 10 11 12..."
(( duration1 = $SECONDS ))            ## time at start of population run

cosmic-pop --final-kstar1 12 --final-kstar2 10 12 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME ONe+any = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + WD Systems (kstar1 = 13; kstar2 = 10,11,12)
echo "now working on NS+WD: 13 and 10,11,12..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 10 12 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + WD Systems (kstar1 = 14; kstar2 = 10,11,12)
# echo "now working on BH+WD: 14 and 10,11,12..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 10 12 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+WD = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo

# ---- NS + NS Systems (kstar1 = 13; kstar2 = 13)
echo "now working on NS+NS: 13 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 13 --final-kstar2 13 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- NS + BH Systems (kstar1 = 14; kstar2 = 13)
echo "now working on BH+NS: 14 and 13..."
(( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 13 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

(( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
echo "   RUNTIME NS+NS = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo

# ---- BH + BH Systems (kstar1 = 14; kstar2 = 14)
# echo "now working on BH+BH: 14 and 14..."
# (( duration1 = $SECONDS ))            ## time at start of population run

# cosmic-pop --final-kstar1 14 --final-kstar2 14 --inifile ../buildFiles/thin_Params.ini --Nstep 10000 --Niter 1000000000 -n 12 >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

# (( duration2 = $SECONDS - duration1))            ## time at end of population run and echo to screen
# echo "   RUNTIME BH+BH = $(($duration2/86400)) DAY $((($duration2 % 86400)/3600)) HR $(((($duration2 % 86400)%3600)/60)) MIN $(((($duration2 % 86400)%3600)%60)) SEC"
# echo "   Finished at: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
# echo


# --- note use of MV instead of mv -- how I have my unix set up.
#     for standard unix use:
mv ./dat_* ../fixed/fixed_thinDisk/.
mv ./log_* ../fixed/fixed_thinDisk/.

# MV ./dat_* ../fixed/fixed_thinDisk/.
# MV ./log_* ../fixed/fixed_thinDisk/.

echo "Thin Disk complete."
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# END OF SCRIPT
# ---------------------------------------------------

NOW=$(date)
echo "FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"

(( duration = $SECONDS ))            ## set time at end of script
echo "TOTAL RUNTIME = $(($duration/86400)) DAY $((($duration % 86400)/3600)) HR $(((($duration % 86400)%3600)/60)) MIN $(((($duration % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo "All done!"

