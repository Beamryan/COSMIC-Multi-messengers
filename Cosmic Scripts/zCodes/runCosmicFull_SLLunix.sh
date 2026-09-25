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
# SCRIPT   : runCosmicFull
# DATE     : 17 Mar 2022
# MODIFIED : 17 Mar 2022 (v1)
#
# This script is a full pipeline script for making a COSMIC galaxy and doing
# an initial LISA analysis for me.  It combines versions of several other
# scripts and codes -- it sets up directories, then runs the codes in
# sequence so the output for each feeds the successive codes:
#    -- makeFixedPops.py    Takes params.ini files and makes gx components fixed pop files
#    -- makeGalaxy.py       Takes fixed pop files then makes a MC galaxy file
#    -- makeCOSMICcsv.py    Takes hdf5 gx files and makes CSV files for C codes
#    -- gxProcess_COSMIC_2021   LISA analysis of galaxy by mono/chirp/confused
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
echo "COSMIC START TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
echo "------------------------------"

(( SECONDS = 0 ))          ## set time at start of script

# ---------------------------------------------------
# Set Up Directories
# This sets up all the directories I use for my full
# galaxy build, on the assumption that making the
# fixed population is the first step. All codes that
# follow this one assume this directory structure.
# This assumes the script exists and is being run from
# a directory called buildCodes/, and makes these
# directories in the parent directory
# ---------------------------------------------------

mkdir ../fixed

mkdir ../gx_data

mkdir ../csv

mkdir ../fixed/fixed_bulge

mkdir ../fixed/fixed_thickDisk

mkdir ../fixed/fixed_thinDisk

echo "Directories created."
echo "------------------------------"
echo

# ---------------------------------------------------
# MAKE FIXED POPS
# ---------------------------------------------------

echo "Making Fixed Populations."
echo

./makeFixedPops_SLLunix.sh >> ../errLogs/log01_Fixed.log 2>> ../errLogs/err01_Fixed.log < /dev/null

echo "Fixed Populations complete."
NOW=$(date)
echo "FIXED FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( durationFixed = $SECONDS ))       ## set time at end of fixed populations
echo "FIXED RUNTIME = $(($durationFixed/86400)) DAY $((($durationFixed % 86400)/3600)) HR $(((($durationFixed % 86400)%3600)/60)) MIN $(((($durationFixed % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# MAKE GALAXY
# Note the files here are saved to the directory
# ../gx_dat/  which is specified in the python code
# makeGalaxy_v9.py  in the variable dat_save_path
# That code reads in a file "../buildFiles/gxModel.txt
# which points to the directories with the fixed
# populations in them.
# ---------------------------------------------------

echo "Making Full Galaxy from Fixed Populations."
echo

python makeGalaxy.py >> ../errLogs/log02_Galaxy.log 2>> ../errLogs/err02_Galaxy.log < /dev/null

echo "Full Galaxy complete."
NOW=$(date)
echo "GALAXY FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( durationGalaxy = $SECONDS - $durationFixed))  ## galaxy runtime in seconds
echo "GALAXY RUNTIME = $(($durationGalaxy/86400)) DAY $((($durationGalaxy % 86400)/3600)) HR $(((($durationGalaxy % 86400)%3600)/60)) MIN $(((($durationGalaxy % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# MAKE CSV
# ---------------------------------------------------

echo "Making CSV files from Galaxy."
echo

python makeCOSMICcsv.py >> ../errLogs/log03_CSV.log 2>> ../errLogs/err03_CSV.log < /dev/null

echo "CSV generation complete."
NOW=$(date)
echo "CSV FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( durationCSV = $SECONDS - $durationFixed - $durationGalaxy)) ## csv runtime in seconds
echo "CSV RUNTIME = $(($durationCSV/86400)) DAY $((($durationCSV % 86400)/3600)) HR $(((($durationCSV % 86400)%3600)/60)) MIN $(((($durationCSV % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# MAKE LISA ANALYSIS
# ---------------------------------------------------

echo "Making LISA sort on CSV Galaxy."
echo

gcc -lm ./gxProcess_COSMIC.c

./a.out >> ../errLogs/log04_LISA.log 2>> ../errLogs/err04_LISA.log < /dev/null

echo "LISA analysis complete."
NOW=$(date)
echo "LISA FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( durationLISA = $SECONDS - $durationFixed - $durationGalaxy - $durationCSV))      ## runtime LISA analysis in seconds
echo "LISA RUNTIME = $(($durationLISA/86400)) DAY $((($durationLISA % 86400)/3600)) HR $(((($durationLISA % 86400)%3600)/60)) MIN $(((($durationLISA % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo


# ---------------------------------------------------
# END OF SCRIPT
# ---------------------------------------------------

echo "COSMIC FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( duration = $SECONDS ))            ## set time at end of script
echo "TOTAL RUNTIME = $(($duration/86400)) DAY $((($duration % 86400)/3600)) HR $(((($duration % 86400)%3600)/60)) MIN $(((($duration % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo "All done!"

