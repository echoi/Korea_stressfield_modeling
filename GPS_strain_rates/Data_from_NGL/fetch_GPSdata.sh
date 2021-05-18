#!/usr/bin/bash

stations=('DAEJ' 'JEJU' 'KOHG' 'MKPO' 'MLYN' 'SBAO' 'SKCH' 'SKMA' 'ANSG' 'BOEN' 'CHCN' 'CHEN' 'CHJU' 'CHLW' 'CHNG' 'CHSG' 'CHWN' 'CHYG' 'CNJU' 'DOKD' 'DOND' 'EOCH' 'GOCH' 'GSAN' 'HADG' 'HOMI' 'HONC' 'INCH' 'INJE' 'JAHG' 'JEOJ' 'JINJ' 'JUNG' 'JUNJ' 'KANR' 'KIMC' 'KUNW' 'KWNJ' 'MARA' 'MUJU' 'NAMW' 'NONS' 'PAJU' 'PUSN' 'SEOS' 'SNJU' 'SONC' 'SOUL' 'SUWN' 'TABK' 'TEGN' 'ULLE' 'WNJU' 'WOLS' 'WPWN' 'WULJ' 'YANP' 'YECH' 'YONK' 'YOWL' 'SEJO' 'GOJE' 'KUSN' 'YODK' 'JIND' 'YOIN' 'SEJN' 'GANH' 'BONH' 'DONH' 'CHUL' 'DANJ' 'GOSG' 'HCHN' 'SMAN')

for s in "${stations[@]}";
do
    fname=${s}.EU.tenv3
    #echo $s $fname
    cmd="wget http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/${fname}"
    echo $cmd
    $cmd
done

