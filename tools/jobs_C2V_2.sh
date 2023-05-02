#!/bin/bash
declare -a BKGDIRS=("/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_2oslep_2ak4_1ak8_v2"
                    "/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_2oslep_2ak4_1ak8"
                    "/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_1lep_1ak8_2ak4_v1")

declare -a SIGDIRS=("/ceph/cms/store/user/jguiang/VBSVHSkim/sig_2oslep_2ak4_1ak8_v2"
                    "/ceph/cms/store/user/jguiang/VBSVHSkim/sig_2oslep_2ak4_1ak8"
                    "/ceph/cms/store/user/jguiang/VBSVHSkim/sig_1lep_1ak8_2ak4_v1")

declare -a FILENAMES=("2oslep_2ak4_1ak8_v2"
                      "2oslep_2ak4_1ak8"
                      "1lep_1ak8_2ak4_v1")

BKGYEARS="20UL16NanoAODv9 \
20UL16NanoAODAPVv9 \
20UL17NanoAODv9 \
20UL18NanoAODv9"


SIGYEARS="20UL16APV-106X \
20UL16-106X \
20UL17-106X \
20UL18-106X"

# Desired backgrounds
BKGPROCESSES=$( jq 'keys' /home/users/joytzphysics/Analysis/data/xsections.json | tr -d '[],"' )

SIGPROCESSES="VBSWZH_Inclusive_4f \
VBSWWH_Inclusive_4f \
VBSZZH_Inclusive_4f \
VBSOSWWH_Inclusive_4f"

for i in "${FILENAMES[@]}"; do
    rm -f /home/users/joytzphysics/Analysis/.semiMerge_jobs_"$i".txt
    rm -f /home/users/joytzphysics/Analysis/.fullyMerge_jobs_"$i".txt
    rm -f /home/users/joytzphysics/Analysis/data/weights_sig_"$i".txt
    rm -f /home/users/joytzphysics/Analysis/data/weights_bkg_"$i".txt
done

for i in {0,1,2}; do
    BKGDIR=${BKGDIRS[$i]}
    echo ${BKGDIR}
    SIGDIR=${SIGDIRS[$i]}
    # all the xsections are in pb
    for BKGPROCESS in ${BKGPROCESSES}; do
        echo $BKGPROCESS
        # the following command can also deal with special characters such as dash
        XSECTION=$( jq --arg keyvar "$BKGPROCESS" '.[$keyvar]' /home/users/joytzphysics/Analysis/data/xsections.json )
        for YEAR in ${BKGYEARS}; do
            if [[ ${YEAR} == "20UL16NanoAODv9" ]]; then LUMINOSITY=16.81; fi
            if [[ ${YEAR} == "20UL16NanoAODAPVv9" ]]; then LUMINOSITY=19.52; fi
            if [[ ${YEAR} == "20UL17NanoAODv9" ]]; then LUMINOSITY=41.48; fi
            if [[ ${YEAR} == "20UL18NanoAODv9" ]]; then LUMINOSITY=59.83; fi
            SUMMEDWGT=0
            RAWEVENTS=0
            for INPUTDIR in $(ls ${BKGDIR} | grep ^${BKGPROCESS} | grep ${YEAR}); do
                echo ${BKGDIR}/${INPUTDIR}
                WGT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR} -c 'W' )
                RAWEVENT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR} -c 'R' )
                SUMMEDWGT=$( echo "$SUMMEDWGT + $WGT" | bc )
                RAWEVENTS=$( echo "$RAWEVENTS + $RAWEVENT" | bc )
            done
            if [[ $SUMMEDWGT == 0 ]]; then echo Directory of the process not found; fi
            if [[ $SUMMEDWGT != 0 ]]
            then
                SCALE1FB=$( echo "scale=15; ${XSECTION} * ${LUMINOSITY} * 1000 / ${SUMMEDWGT}" | bc )
                echo ${BKGDIR}/${INPUTDIR},${BKGPROCESS},${YEAR},${LUMINOSITY},${XSECTION},${SUMMEDWGT},${RAWEVENTS},${SCALE1FB} >> /home/users/joytzphysics/Analysis/data/weights_bkg_${FILENAMES[$i]}.txt
            fi
            IFILE=0
            for INPUTDIR in $(ls ${BKGDIR} | grep ${BKGPROCESS} | grep ${YEAR}); do
                for FILE in $(ls ${BKGDIR}/${INPUTDIR}); do
                    IFILE=$(( IFILE+1 ))
                    echo "./semiMerge -t Events -d outputs/output_semiMerge_${FILENAMES[$i]} -s ${SCALE1FB} -n ${BKGPROCESS}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/outputs/semiMerge_logfiles_${FILENAMES[$i]}/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> /home/users/joytzphysics/Analysis/.semiMerge_jobs_${FILENAMES[$i]}.txt
                    echo "./fullyMerge -t Events -d outputs/output_fullyMerge_${FILENAMES[$i]} -s ${SCALE1FB} -n ${BKGPROCESS}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/outputs/fullyMerge_logfiles_${FILENAMES[$i]}/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> /home/users/joytzphysics/Analysis/.fullyMerge_jobs_${FILENAMES[$i]}.txt
                done
            done
        done
    done

    echo ${SIGDIR}
    for SIGPROCESS in ${SIGPROCESSES}; do
        if [[ ${SIGPROCESS} == "VBSWZH_Inclusive_4f" ]]
        then 
            XSECTION=1.67403
            KEY="WZH"
        fi
        if [[ ${SIGPROCESS} == "VBSWWH_Inclusive_4f" ]]
        then
            XSECTION=1.50253
            KEY="WWH"
        fi
        if [[ ${SIGPROCESS} == "VBSZZH_Inclusive_4f" ]]
        then
            XSECTION=1.06028
            KEY="ZZH"
        fi
        if [[ ${SIGPROCESS} == "VBSOSWWH_Inclusive_4f" ]]
        then
            XSECTION=2.52444
            KEY="OSWWH"
        fi
        for YEAR in ${SIGYEARS}; do
            if [[ ${YEAR} == "20UL16APV-106X" ]]; then LUMINOSITY=19.52; fi
            if [[ ${YEAR} == "20UL16-106X" ]]; then LUMINOSITY=16.81; fi
            if [[ ${YEAR} == "20UL17-106X" ]]; then LUMINOSITY=41.48; fi
            if [[ ${YEAR} == "20UL18-106X" ]]; then LUMINOSITY=59.83; fi

            SUMMEDWGT=0
            RAWEVENTS=0
            # one signal sample one year loop
            # calculate weights
            for INPUTDIR in $(ls ${SIGDIR} | grep ^${SIGPROCESS} | grep ${YEAR}); do
                echo ${SIGDIR}/${INPUTDIR}/${FILE}
                WGT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'W' )
                RAWEVENT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'R' )
                SUMMEDWGT=$( echo "$SUMMEDWGT + $WGT" | bc )
                RAWEVENTS=$( echo "$RAWEVENTS + $RAWEVENT" | bc )
            done
            SCALE1FB=$( echo "scale=10; ${XSECTION} * ${LUMINOSITY} / ${SUMMEDWGT}" | bc )
            echo ${SIGDIR}/${INPUTDIR},${KEY},${YEAR},${LUMINOSITY},${XSECTION},${SUMMEDWGT},${RAWEVENTS},${SCALE1FB}>> /home/users/joytzphysics/Analysis/data/weights_sig_${FILENAMES[$i]}.txt
            
            IFILE=0
            for INPUTDIR in $(ls ${SIGDIR} | grep ^${SIGPROCESS} | grep ${YEAR}); do
                for FILE in $(ls ${SIGDIR}/${INPUTDIR}); do
                    IFILE=$(( IFILE+1 ))
                    echo "./semiMerge -t Events -d outputs/output_semiMerge_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/outputs/semiMerge_logfiles_${FILENAMES[$i]}/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> /home/users/joytzphysics/Analysis/.semiMerge_jobs_${FILENAMES[$i]}.txt
                    echo "./fullyMerge -t Events -d outputs/output_fullyMerge_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/outputs/fullyMerge_logfiles_${FILENAMES[$i]}/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> /home/users/joytzphysics/Analysis/.fullyMerge_jobs_${FILENAMES[$i]}.txt
                done
            done
        done
    done
    python3 /home/users/joytzphysics/Analysis/tools/tocsv.py -f /home/users/joytzphysics/Analysis/data/weights_sig_${FILENAMES[$i]}
    python3 /home/users/joytzphysics/Analysis/tools/tocsv.py -f /home/users/joytzphysics/Analysis/data/weights_bkg_${FILENAMES[$i]}
done

