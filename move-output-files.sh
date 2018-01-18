#!/bin/bash

# basic function move syntax with error check / results print
move_function()
{
    if [ "$#" -ne "2" ]
    then
        echo "Usage $0 from to"
    else

        FROM="$1"
        TO="$2"

        echo '$from = '"$FROM"
        echo '$to = '"$TO"

        if [ -d "$FROM" ]
        then
            if [ -d "$TO" ]
            then
                
                COUNT_GOOD=0
                COUNT_BAD=0
                for FILE in $(cat file-list.txt)
                do
                    #echo "$FILE"
                    mv "$FROM/$FILE" "$TO/$FILE"
                    RET="$?"

                    if [ "$RET" -eq "0" ]
                    then
                        (( COUNT_GOOD += 1 ))
                    else
                        (( COUNT_BAD += 1 ))
                    fi
                done

                echo "$COUNT_GOOD files moved, $COUNT_BAD errors"
        
            else
                echo "$TO does not exist or is not a directory"    
            fi
        else
            echo "$FROM does not exist or is not a directory"
        fi

    fi
}

# main function
main()
{
    if [ "$#" -eq 2 ]
    then
        # syntax is move-output-files.sh [before/after] [mc/data]

        FROM="./bin" # from dir is always default output dir
        TO="./bin" # to dir is changed later
        
        TIME="$1" # before/after (changing the PPT)
        TYPE="$2" # data/mc

        #DATAFILE=""

        if [ "$TIME" == "before" ] || [ "$TIME" == "after" ]
        then
            if [ "$TYPE" == "data" ] || [ "$TYPE" == "mc" ]
            then
                GOOD="GOOD"
                
                # switch datafile depending on mc/data
                if [ "$TYPE" == "data" ]
                then
                    DATAFILE="cell8_out.root"
                elif [ "$TYPE" == "mc" ]
                then
                    DATAFILE="timestamp_cpp_out.root"
                fi
            else
                echo "Error incorrect TYPE $TYPE"
            fi

        else
            echo "Error incorrect TIME $TIME"
        fi

        if [ "$GOOD" == "GOOD" ]
        then
        
            # set before/after dir
            if [ "$TIME" == "before" ]
            then
                TO="$TO/before-change-ppt"
                
            elif [ "$TIME" == "after" ]
            then
                TO="$TO/after-change-ppt"
            # can't be anything else, already performed this check
            fi
            
            # set data/mc dir
            if [ "$TYPE" == "data" ]
            then
                TO="$TO/data"
            elif [ "$TYPE" == "mc" ]
            then
                TO="$TO/mc"
            # can't be anything else, already performed this check
            fi
        
            mv "$FROM/$DATAFILE" "$TO/$DATAFILE"
            move_function "$FROM" "$TO"
        fi

        # move the output data file
        #if [ "$TYPE" == "data" ]
        #then
        #    mv "$FROM/cell8_out.root" "$TO/cell8_out.root"
        #    move_function $FROM $TO
        #elif [ "$TYPE" == "mc" ]
        #then
        #    mv "$FROM/timestamp_cpp_out.root" "$TO/timestamp_cpp_out.root"
        #    move_function $FROM $TO
        #else
        #    echo "Error incorrect TYPE $TYPE"
        #fi
    else
        echo "Usage $0 TIME TYPE"
    fi
}

main $@
