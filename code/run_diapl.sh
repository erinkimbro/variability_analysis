#!/bin/bash


while read -r RA DEC; do

    cd /Users/erinkimbro/Projects/merian_variable/code

    #need to add a cut on minimum number of observations
    #possibly do this in xmatching for selection of targets
    #ADD: throw out obs that are not 300x300
    #ADD: exit analysis if # of frames drops below 20 
    python ztf_data_download.py "$RA" "$DEC"

    #run diapl

    cd ../WorkingDir

    num="$(ls -1 | wc -l)" 
    echo "$num"

    if ((num>42))
    then
        cd ../code
        python get_image_list.py

        cd ../WorkingDir
        bash fwhms.bash

        #insert astrometric reference in diapl setup par
        cd ../code
        python set_ref_image.py

        cd ../WorkingDir
        bash mktpllist.bash

        bash shifts.bash

        bash template.bash

        bash subtract1.bash

        cd ../code

        #run photometry program
        #create frame work for tpl data
        python photometry.py "$RA" "$DEC"

        #run qso fit
        #create frame work for fit data
        python run_qsofit.py "$RA" "$DEC"

        #produce lightcurve
        python lightcurve.py "$RA" "$DEC"

        #pipeline mostly working line by line 
        #need to fine tune plot

        #flush data
        rm -r ../WorkingDir

        cp -R ../WorkingDirectoryTemplate ../WorkingDir
        
    else
        echo "not enough exposures"
        #flush data
        rm -r ../WorkingDir

        cp -R ../WorkingDirectoryTemplate ../WorkingDir

        fi


done < pos_var_sample.txt