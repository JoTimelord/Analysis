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

SIGPROCESSES="VBSWZH_Inclusive_4f"

rm -f .jobs.txt
rm -f weightssig.txt
rm -f weightsbkg.txt

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
    IFILE=0
    for YEAR in ${BKGYEARS}; do
        if [[ ${YEAR} == "20UL16NanoAODv9" ]]; then LUMINOSITY=16.81; fi
        if [[ ${YEAR} == "20UL16NanoAODAPVv9" ]]; then LUMINOSITY=19.52; fi
        if [[ ${YEAR} == "20UL17NanoAODv9" ]]; then LUMINOSITY=41.48; fi
        if [[ ${YEAR} == "20UL18NanoAODv9" ]]; then LUMINOSITY=59.83; fi
        for INPUTDIR in $(ls ${BKGDIR} | grep ^${KEY} | grep ${YEAR}); do
            echo ${INPUTDIR}
            IFILE=$(( IFILE+1 ))
            SUMMEDWGT=$( python3 scale1fb.py -f ${BKGDIR}/${INPUTDIR} )
            SCALE1FB=$( echo "scale=10; ${XSECTION} * ${LUMINOSITY} / ${SUMMEDWGT}" | bc )
            echo ${BKGPROCESS} ${KEY} ${YEAR} ${LUMINOSITY} ${XSECTION} ${SUMMEDWGT} ${SCALE1FB} >> weightsbkg.txt
            echo "./runLooper -t Events -d output -s ${SCALE1FB} -n ${KEY}_${YEAR} -T tree ${SIGDIR}/${INPUTDIR}/merged.root> ~/Analysis/logfiles/${BKGPROCESS}_${IFILE}.log 2>&1" >> .jobs.txt
        done
    done
done


for SIGPROCESS in ${SIGPROCESSES}; do
    IFILE=0
    XSECTION=1.67403
    for YEAR in ${SIGYEARS}; do
        if [[ ${YEAR} == "20UL16APV" ]]; then LUMINOSITY=19.52; fi
        if [[ ${YEAR} == "20UL16" ]]; then LUMINOSITY=16.81; fi
        if [[ ${YEAR} == "20UL17" ]]; then LUMINOSITY=41.48; fi
        if [[ ${YEAR} == "20UL18" ]]; then LUMINOSITY=59.83; fi
        for INPUTDIR in $(ls ${SIGDIR} | grep ${SIGPROCESS} | grep ${YEAR}); do
            echo ${INPUTDIR}
            IFILE=$(( IFILE+1 ))
            SUMMEDWGT=$( python3 scale1fb.py -f ${SIGDIR}/${INPUTDIR} )
            SCALE1FB=$( echo "scale=10; ${XSECTION} * ${LUMINOSITY} / ${SUMMEDWGT}" | bc )
            echo ${SIGPROCESS} ${KEY} ${YEAR} ${LUMINOSITY} ${XSECTION} ${SUMMEDWGT} ${SCALE1FB} >> weightssig.txt
            echo "./runLooper -t Events -d output -s ${SCALE1FB} -n WZH_${YEAR} -T tree ${SIGDIR}/${INPUTDIR}/merged.root > ~/Analysis/logfiles/${SIGPROCESS}_${IFILE}.log 2>&1">> .jobs.txt
        done
    done
done


