#!/bin/bash

BKGDIR=/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_1lep_1ak8_2ak4_v1
SIGDIR=/ceph/cms/store/user/jguiang/VBSVHSkim/sig_1lep_1ak8_2ak4_v1

BKGYEARS="20UL16NanoAODv9 \
20UL16NanoAODAPVv9 \
20UL17NanoAODv9 \
20UL18NanoAODv9"

SIGYEARS="20UL16APV \
20UL16 \
20UL17 \
20UL18"

# Desired backgrounds
BKGPROCESSES="DYJETSbkg \
WWdilep \
WWinclusive \
ttdilep"

SIGPROCESSES="VBSWZH_C2V_3 \
VBSZZH_C2V_3 \
VBSOSWWH_C2V_3 \
VBSOSWWH_C2V_4 \
VBSZZH_C2V_4 \
VBSWWH_C2V_4"

rm -f .jobs1_C2V_3and4.txt
rm -f .jobs2_C2V_3and4.txt
rm -f weights_sig_C2V_3and4.txt
rm -f weights_bkg_C2V_3and4.txt


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
        let SUMMEDWGT=0.00
        let RAWEVENTS=0.00
        let WGT=0.00
        let RAWEVENT=0.00
        for INPUTDIR in $(ls ${BKGDIR} | grep ^${KEY} | grep ${YEAR}); do
            for FILE in $(ls ${BKGDIR}/${INPUTDIR}); do
                WGT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR}/${FILE} -c 'W' )
                RAWEVENT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR}/${FILE} -c 'R' )
                SUMMEDWGT=`expr $SUMMEDWGT + $WGT`
                RAWEVENTS=`expr $RAWEVENTS + $RAWEVENT`
            done
        done
        echo ${SUMMEDWGT}
        SCALE1FB=$( echo "scale=10; ${XSECTION}/${SUMMEDWGT}*${LUMINOSITY}" | bc )
        echo ${BKGPROCESS} ${KEY} ${YEAR} ${LUMINOSITY} ${XSECTION} ${SUMMEDWGT} ${SCALE1FB} >> weights_bkg_C2V_3and4.txt
        
        IFILE=0
        for INPUTDIR in $(ls ${BKGDIR} | grep ^${KEY} | grep ${YEAR}); do
            for FILE in $(ls ${BKGDIR}/${INPUTDIR}); do
                echo ${BKGDIR}/${INPUTDIR}/${FILE}
                IFILE=$(( IFILE+1 ))                
                echo "./runLooper1 -t Events -d output1_C2V_3and4 -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> .jobs1_C2V_3and4.txt
                echo "./runLooper2 -t Events -d output2_C2V_3and4 -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${BKGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${BKGPROCESS}_${YEAR}_${IFILE}.log 2>&1" >> .jobs2_C2V_3and4.txt
            done
        done
    done
done


for SIGPROCESS in ${SIGPROCESSES}; do
    if [[ ${SIGPROCESS} == *"VBSOSWWH_C2V_3" ]]
    then
        XSECTION=5.652
        KEY="VBSOSWWH_incl_C2V_3"
    fi
    if [[ ${SIGPROCESS} == *"VBSOSWWH_C2V_4" ]]
    then
        XSECTION=11.8
        KEY="VBSOSWWH_incl_C2V_4"
    fi
    if [[ ${SIGPROCESS} == *"VBSWZH_C2V_3" ]]
    then
        XSECTION=3.742
        KEY="VBSWZH_incl_C2V_3"
    fi
    if [[ ${SIGPROCESS} == *"VBSWZH_C2V_4" ]]
    then
        XSECTION=7.865
        KEY="VBSWZH_incl_C2V_4"
    fi
    if [[ ${SIGPROCESS} == *"VBSZZH_C2V_4" ]]
    then
        XSECTION=6.592
        KEY="VBSZZH_incl_C2V_4"
    fi
    if [[ ${SIGPROCESS} == *"VBSZZH_C2V_3" ]]
    then
        XSECTION=2.994
        KEY="VBSZZH_incl_C2V_3"
    fi
    
    for YEAR in ${SIGYEARS}; do
        if [[ ${YEAR} == "20UL16APV" ]]; then LUMINOSITY=19.52; fi
        if [[ ${YEAR} == "20UL16" ]]; then LUMINOSITY=16.81; fi
        if [[ ${YEAR} == "20UL17" ]]; then LUMINOSITY=41.48; fi
        if [[ ${YEAR} == "20UL18" ]]; then LUMINOSITY=59.83; fi
        SUMMEDWGT=0
        RAWEVENTS=0
        for INPUTDIR in $(ls ${SIGDIR} | grep ${KEY} | grep ${YEAR}); do
            for FILE in $(ls ${SIGDIR}/${INPUTDIR}); do
                WGT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'W' )
                RAWEVENT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR}/${FILE} -c 'R' )
                SUMMEDWGT=$(( $SUMMEDWGT + $WGT ))
                RAWEVENTS=$(( $RAWEVENTS + $RAWEVENT))
            done
        done
        SCALE1FB=$( echo "scale=10; ${XSECTION}/${SUMMEDWGT}*${LUMINOSITY}" | bc )
        echo ${SIGPROCESS} ${KEY} ${YEAR} ${LUMINOSITY} ${XSECTION} ${SUMMEDWGT} ${SCALE1FB} >> weights_sig_C2V_3and4.txt
        
        IFILE=0
        for INPUTDIR in $(ls ${SIGDIR} | grep ${KEY} | grep ${YEAR}); do
            for FILE in $(ls ${SIGDIR}/${INPUTDIR}); do
                echo ${SIGDIR}/${INPUTDIR}/${FILE}
                IFILE=$(( IFILE+1 ))
                echo "./runLooper2 -t Events -d output2_C2V_3and4 -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> .jobs2_C2V_3and4.txt
                echo "./runLooper1 -t Events -d output1_C2V_3and4 -s ${SCALE1FB} -n ${KEY}_${YEAR}_${IFILE} -T tree ${SIGDIR}/${INPUTDIR}/${FILE} > ~/Analysis/logfiles/${SIGPROCESS}_${YEAR}_${IFILE}.log 2>&1">> .jobs1_C2V_3and4.txt
            done
        done
    done
done


