
while read -r FIELD CCDID QID FILTER; do

    python ztf_download_field.py "$FIELD" "$CCDID" "$QID" "$FILTER"

    while read line; do

        NUM=15

        echo "$NUM"

        echo "$line"

        if [[ $line -ge $NUM ]]; then

            echo "Enough Data"

            cd ../../../results_"$FILTER"

            mkdir "$FIELD"_"$CCDID"_"$QID"

            cd "$FIELD"_"$CCDID"_"$QID"

            mkdir tables
            mkdir plots

            cd tables

            mkdir measurements
            mkdir models

            cd ../plots

            mkdir cutouts
            mkdir light_curves

            cd ../../../variability_analysis/code/field

            python ztf_dwnld_ims_field.py

            cd ..

            python get_image_list.py

            cd field
            python diapl_setup.py

            cd ../../../diapl2/WorkingDir
            bash fwhms.bash
            cd ../../variability_analysis/code  
            python set_ref_image.py

            cd ../../diapl2/WorkingDir   
            bash mktpllist.bash
            bash shifts.bash   
            bash template.bash 
            bash subtract1.bash

            if test -f 1_1/simages.txt
            then
                cd ../../variability_analysis/code/field

                python photometry_field.py "$FIELD" "$CCDID" "$QID" "$FILTER"
                python qsofit_field.py "$FIELD" "$CCDID" "$QID" "$FILTER"
                python lightcurve_field.py "$FIELD" "$CCDID" "$QID" "$FILTER"

                rm -rf ../../../diapl2/WorkingDir
                cp -R ../../../diapl2/WorkingDirectoryTemplate ../../../diapl2/WorkingDir

                cd ../../../results_"$FILTER"
                echo "$FIELD" "$CCDID" "$QID" "Successful" >> field_log.txt
                cd ../variability_analysis/code/field

            else
                cd ../../results_"$FILTER"
                rm -rf "$FIELD"_"$CCDID"_"$QID"

                cd ../variability_analysis/code/field
                rm -rf ../../../diapl2/WorkingDir
                cp -R ../../../diapl2/WorkingDirectoryTemplate ../../../diapl2/WorkingDir
                cd ../../../results_"$FILTER"
                echo "$FIELD" "$CCDID" "$QID" "$FILTER" "DIAPL FAILED" >> field_log.txt
                cd ../variability_analysis/code/field
            fi


        else
            echo "Not Enough Data"
            cd ../../results_"$FILTER"
            rm -rf "$FIELD"_"$CCDID"_"$QID"

            cd ../variability_analysis/code/field
            rm -rf ../../../diapl2/WorkingDir
            cp -R ../../../diapl2/WorkingDirectoryTemplate ../../../diapl2/WorkingDir
            cd ../../../results_"$FILTER"
            echo "$FIELD" "$CCDID" "$QID" "NOT ENOUGH DATA" >> field_log.txt
            cd ../variability_analysis/code/field
            

        fi

    done < test.txt

done < ../../../workspace/data/fields.dat