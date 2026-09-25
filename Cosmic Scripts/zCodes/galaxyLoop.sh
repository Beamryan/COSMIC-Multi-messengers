mkdir -p ../loopFiles
cd ../loopFiles

for i in {9..20};
do
    # School PC
    #mkdir -p "$i/gx_data" "$i/csv" "$i/errLogs"

    # Home PC
    mkdir -p "../../../../../../../e/loopFiles/$i/gx_data" "../../../../../../../e/loopFiles/$i/csv" "../../../../../../../e/loopFiles/$i/errLogs"

    NOW=$(date)
    echo "COSMIC START TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
    echo "------------------------------"

    (( SECONDS = 0 ))          ## set time at start of script

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

    python ../zCodes/makeGalaxy.py >> ../errLogs/log02_Galaxy.log 2>> ../errLogs/err02_Galaxy.log < /dev/null

    echo "Full Galaxy complete."
    NOW=$(date)
    echo "GALAXY FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
    (( durationGalaxy = $SECONDS))  ## galaxy runtime in seconds
    echo "GALAXY RUNTIME = $(($durationGalaxy/86400)) DAY $((($durationGalaxy % 86400)/3600)) HR $(((($durationGalaxy % 86400)%3600)/60)) MIN $(((($durationGalaxy % 86400)%3600)%60)) SEC"
    echo
    echo "------------------------------"
    echo


    # ---------------------------------------------------
    # MAKE CSV
    # ---------------------------------------------------

    echo "Making CSV files from Galaxy."
    echo

    python ../zCodes/makeCOSMICcsv.py >> ../errLogs/log03_CSV.log 2>> ../errLogs/err03_CSV.log < /dev/null

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
   
    gcc -g ../zCodes/gxProcess_COSMIC_forResolved.c -o ../zCodes/gxProcess_COSMIC_forResolved -lm

    ../zCodes/gxProcess_COSMIC_forResolved >> ../errLogs/log04_LISA.log 2>> ../errLogs/err04_LISA.log < /dev/null

    echo "LISA analysis complete."
    NOW=$(date)
    echo "LISA FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
    (( durationLISA = $SECONDS - $durationFixed - $durationGalaxy - $durationCSV))      ## runtime LISA analysis in seconds
    echo "LISA RUNTIME = $(($durationLISA/86400)) DAY $((($durationLISA % 86400)/3600)) HR $(((($durationLISA % 86400)%3600)/60)) MIN $(((($durationLISA % 86400)%3600)%60)) SEC"
    echo
    echo "------------------------------"
    echo

    # School PC
    #mv ../gx_data/* "$i/gx_data/"
    #mv ../csv/* "$i/csv/"
    #mv ../errLogs/* "$i/errLogs/"    

    # Home PC 
    mv ../gx_data/* "../../../../../../../e/loopFiles/$i/gx_data"
    mv ../csv/* "../../../../../../../e/loopFiles/$i/csv/"
    mv ../errLogs/* "../../../../../../../e/loopFiles/$i/errLogs/"
    
done
    

# ---------------------------------------------------
# END OF SCRIPT
# ---------------------------------------------------

echo "COSMIC FINISH TIME: $(date +"%a %d %b %Y  %I:%M:%S %p %Z")"
(( duration = $SECONDS ))            ## set time at end of script
echo "TOTAL RUNTIME = $(($duration/86400)) DAY $((($duration % 86400)/3600)) HR $(((($duration % 86400)%3600)/60)) MIN $(((($duration % 86400)%3600)%60)) SEC"
echo
echo "------------------------------"
echo "All done!"

