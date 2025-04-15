cd /Users/erinkimbro/Projects/merian_variable/code/field

while read -r FIELD CCDID QID; do

    python ztf_download_field.py "$FIELD" "$CCDID" "$QID"

    while read line; do

        NUM=15

        echo "$NUM"

        echo "$line"

        if [[ $line -ge $NUM ]]; then

            echo "Enough Data"

            cd ../../results

            mkdir "$FIELD"_"$CCDID"_"$QID"

            cd "$FIELD"_"$CCDID"_"$QID"

            mkdir tables
            mkdir plot

            cd tables

            mkdir measurements
            mkdir models

            cd ../plot

            mkdir cutouts
            mkdir light_curves

            cd ../../../code/field

            python ztf_dwnld_ims_field.py

            cd ..

            python get_image_list.py

            cd field
            python diapl_setup.py

            cd ../../WorkingDir
            bash fwhms.bash
            cd ../code  
            python set_ref_image.py

            cd ../WorkingDir   
            bash mktpllist.bash
            bash shifts.bash   
            bash template.bash 
            bash subtract1.bash

            if test -f 1_1/simages.txt
            then
                cd ../code/field

                python photometry_field.py "$FIELD" "$CCDID" "$QID"
                python qsofit_field.py "$FIELD" "$CCDID" "$QID"
                python lightcurve_field.py "$FIELD" "$CCDID" "$QID"

                rm -rf ../../WorkingDir
                cp -R ../../WorkingDirectoryTemplate ../../WorkingDir

                echo "$FIELD" "$CCDID" "$QID" "Successful" >> field_log.txt

            else
                cd ../results
                rm -rf "$FIELD"_"$CCDID"_"$QID"

                cd ../code/field
                rm -rf ../../WorkingDir
                cp -R ../../WorkingDirectoryTemplate ../../WorkingDir
                echo "$FIELD" "$CCDID" "$QID" "DIAPL FAILED" >> field_log.txt
            fi


        else
            echo "Not Enough Data"
            rm -rf ../../WorkingDir
            cp -R ../../WorkingDirectoryTemplate ../../WorkingDir
            echo "$FIELD" "$CCDID" "$QID" "NOT ENOUGH DATA" >> field_log.txt
            

        fi

    done < test.txt

done < fields.dat