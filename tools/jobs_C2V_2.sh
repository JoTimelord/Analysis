#!/bin/bash
BKGDIRS=(/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_2oslep_2ak4_1ak8_v2 \
/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_2oslep_2ak4_1ak8 \
/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_1lep_1ak8_2ak4_v1)

SIGDIRS=(/ceph/cms/store/user/jguiang/VBSVHSkim/sig_2oslep_2ak4_1ak8_v2 \
/ceph/cms/store/user/jguiang/VBSVHSkim/sig_2oslep_2ak4_1ak8 \
/ceph/cms/store/user/jguiang/VBSVHSkim/sig_1lep_1ak8_2ak4_v1)

FILENAMES=(2oslep_2ak4_1ak8_v2 \
2oslep_2ak4_1ak8 \
1lep_1ak8_2ak4_v1)

BKGYEARS="20UL16NanoAODv9 \
20UL16NanoAODAPVv9 \
20UL17NanoAODv9 \
20UL18NanoAODv9"

SIGYEARS="20UL16APV-106X \
20UL16-106X \
20UL17-106X \
20UL18-106X"

# Desired backgrounds
BKGPROCESSES="DYJETSbkg \
WWdilep \
WWinclusive \
ttdilep"

SIGPROCESSES="VBSWZH_Inclusive_4f \
VBSWWH_Inclusive_4f \
VBSZZH_Inclusive_4f \
VBSOSWWH_Inclusive_4f"

for NAME in FILENAMES; do
    rm -f .jobs1_${NAME}.txt
    rm -f .jobs2_${NAME}.txt
    rm -f weights_sig_${NAME}.txt
    rm -f weights_bkg_${NAME}.txt
done

for i in {0,1,2}; do
    BKGDIR=${BKGDIRS[$i]}
    SIGDIR=${SIGDIRS[$i]}
    # all the xsections are in pb
    for BKGPROCESS in ${BKGPROCESSES}; do
        if [[ ${BKGPROCESS} == *"DYJETSbkg" ]]
        then
            XSECTION=6197900
            KEY="DYJetsToLL_M-50"
        fi
        if [[ ${BKGPROCESS} == *"WWinclusive" ]]
        then
            XSECTION=118710
            KEY="WW_TuneCP5"
        fi
        if [[ ${BKGPROCESS} == *"WWdilep" ]]
        then
            XSECTION=12178
            KEY="WWTo2L2Nu_TuneCP5"
        fi
        if [[ ${BKGPROCESS} == *"ttdilep" ]]
        then
            XSECTION=72100
            KEY="TTTo2L2Nu_TuneCP5"
        fi
        for YEAR in ${BKGYEARS}; do
            if [[ ${YEAR} == "20UL16NanoAODv9" ]]; then LUMINOSITY=16.81; fi
            if [[ ${YEAR} == "20UL16NanoAODAPVv9" ]]; then LUMINOSITY=19.52; fi
            if [[ ${YEAR} == "20UL17NanoAODv9" ]]; then LUMINOSITY=41.48; fi
            if [[ ${YEAR} == "20UL18NanoAODv9" ]]; then LUMINOSITY=59.83; fi
            IFILE=0
            for INPUTDIR in $(ls ${BKGDIR} | grep ^${KEY} | grep ${YEAR}); do
                for FILE in $(ls ${BKGDIR}/${INPUTDIR}); do
                    echo ${BKGDIR}/${INPUTDIR}/${FILE}
                    IFILE=$(( IFILE+1 ))
                    SUMMEDWGT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR}/${FILE} -c 'W' )
                    RAWEVENTS=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR}/${FILE} -c 'R' )
                    SCALE1FB=$( echo "scale=10; ${XSECTION} * ${LUMINOSITY} / ${SUMMEDWGT}" | bc )
                    echo ${BKGDIR}/${INPUTDIR},${KEY},${YEAR},${LUMINOSITY},${XSECTION},${SUMMEDWGT},${RAWEVENTS},${SCALE1FB} >> weights_bkg_${FILENAMES[$i]}.txt
                    echo "./runLooper2 -t Events -d output2_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> .jobs2_${FILENAMES[$i]}.txt
                    echo "./runLooper1 -t Events -d output1_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> .jobs1_${FILENAMES[$i]}.txt
                done
            done
        done
    done


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
            IFILE=0
            for INPUTDIR in $(ls ${SIGDIR} | grep ^${SIGPROCESS} | grep ${YEAR}); do
                for FILE in $(ls ${SIGDIR}/${INPUTDIR}); do
                    IFILE=$(( IFILE+1 ))
                    echo ${SIGDIR}/${INPUTDIR}/${FILE}
                    SUMMEDWGT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'W' )
                    RAWEVENTS=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'R' )
                    SCALE1FB=$( echo "scale=10; ${XSECTION} * ${LUMINOSITY} / ${SUMMEDWGT}" | bc )
                    echo ${SIGDIR}/${INPUTDIR},${KEY},${YEAR},${LUMINOSITY},${XSECTION},${SUMMEDWGT},${RAWEVENTS},${SCALE1FB}>> weights_sig_${FILENAMES[$i]}.txt
                    echo "./runLooper2 -t Events -d output2_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> .jobs2_${FILENAMES[$i]}.txt
                    echo "./runLooper1 -t Events -d output1_${FILENAMES[$i]} -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> .jobs1_${FILENAMES[$i]}.txt
                done
            done
        done
    done
    python3 tocsv.py -f weights_sig_${FILENAMES[$i]}
    python3 tocsv.py -f weights_bkg_${FILENAMES[$i]}
done

# Convert the signal and background file samples to csv files for import


