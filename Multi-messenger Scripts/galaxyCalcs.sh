for i in {1..5};
do
    # School PC
    #mkdir -p "../CosmicScripts/loopFiles/$i/magnitudeFiles"
    
    # Home PC
    mkdir -p "../../../../../../e/loopFiles/$i/magnitudeFiles"

    echo "Converting CSV Parameters"
    echo
    python csvParameterFix.py ../../../../../../e/loopFiles/$i/csv/gx_petoskey_resolved.csv
    echo "Parameters appended"
    echo
    

    echo "Calculating Apparent Magnitudes and Saving Visible Star Files"
    echo
    python appMagCalc.py ../../../../../../e/loopFiles/$i/csv/gx_petoskey_resolved.csv ../globalDensity.fits
    echo "Visible Star Files Created"
    echo


    echo "Making Galaxy Scatter Plots"
    echo
    python galaxyPlot.py visible_stars50.csv Plots.ini
    echo "Scatter Plots Completed"
    echo


    echo "Making Galaxy 3d Plots"
    echo
    python 3dPlotPng.py visible_stars20.csv 
    python 3dPlotPng.py visible_stars50.csv
    echo "3d Plots Completed"
    echo

    # School PC
    #mv *.png "../CosmicScripts/loopFiles/$i/magnitudeFiles/"
    #mv *.csv "../CosmicScripts/loopFiles/$i/magnitudeFiles/"

    # Home PC
    mv *.png "../../../../../../e/loopFiles/$i/magnitudeFiles/"
    mv *.csv "../../../../../../e/loopFiles/$i/magnitudeFiles/"

done
    
echo
echo "------------------------------"
echo "All done!"

