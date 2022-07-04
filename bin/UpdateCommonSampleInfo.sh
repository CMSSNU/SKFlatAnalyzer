#!/bin/bash
while read line <&3
do
    array=(${line//\// })
    YEAR=${array[6]}
    if [ "${array[7]}" = "DATA" ];then
	NAME=$(echo ${array[8]/.root/}|sed 's/GetEffLumi_//'|grep -o "^[^_]*")
	PERIOD=$(echo ${array[8]/.root/}|sed "s/GetEffLumi_${NAME}_//")
	NUMS=($(echo 'cout<<Form("%d\t",(int)sumW->GetEntries())<<endl;'|root -b -l ${line}|tail -n1))
	SUMMARYFILE=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/SampleSummary_DATA.txt
	if cat $SUMMARYFILE|grep -q "$NAME.*$PERIOD"; then
	    OLDLINE=$(cat $SUMMARYFILE|grep "$NAME.*$PERIOD")
	    if [ "$(echo $OLDLINE|awk '{print $3}')" != "${NUMS[0]}" ];then
		echo -e "$SUMMARYFILE: $OLDLINE\t-->\t$NAME\t$PERIOD\t${NUMS[0]}"
		read -p "(y/n): " YES
		[ "$YES" = "y" ] && sed -i "s/.*$NAME.*$PERIOD.*/$NAME\t$PERIOD\t${NUMS[0]}/" $SUMMARYFILE
	    fi
	else
	    echo "echo -e $NAME\t$PERIOD\t${NUMS[0]} >> $SUMMARYFILE"
	    read -p "(y/n): " YES
	    [ "$YES" = "y" ] && echo -e "$NAME\t$PERIOD\t${NUMS[0]}" >> $SUMMARYFILE
	fi
    else
	## collect info
	NAME=$(echo ${array[7]/.root/}|sed 's/GetEffLumi_//')
	DASNAME=$(head -n1 $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/ForSNU/${NAME}.txt|sed 's@/@ @g'|awk '{print $7}')
	CROSSSECTIONS=($(cat $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/SampleSummary_MC.txt|grep "^$NAME[^a-zA-Z0-9_]"|awk '{print $3}'|grep -v FIXMECROSSSECTION|uniq))
	if [ ${#CROSSSECTIONS[@]} -eq 1 ]; then 
	    CROSSSECTION=${CROSSSECTIONS[0]}
	elif [ ${#CROSSSECTIONS[@]} -eq 0 ]; then
	    CROSSSECTION=""
	else
	    echo "candidate cross sections= ${CROSSSECTIONS[@]}"
	    read -p "select cross section: " CROSSSECTION
	    [ -z "$CROSSSECTION" ] && CROSSSECTION=FIXMECROSSSECTION
	fi
	if [ "$CROSSSECTION" = "" ]; then
	    CROSSSECTIONS=($(cat $SKFlat_WD/data/$SKFlatV/*/Sample/SampleSummary_MC.txt|grep "^$NAME[^a-zA-Z0-9_]"|awk '{print $3}'|grep -v FIXMECROSSSECTION|uniq))
	    if [ ${#CROSSSECTIONS[@]} -eq 1 ];
	    then 
		CROSSSECTION=${CROSSSECTIONS[0]}
	    fi
	fi
	if [ "$CROSSSECTION" = "" ]; then
	    CROSSSECTIONS=($(cat $SKFlat_WD/data/*/*/Sample/SampleSummary_MC.txt|grep "^$NAME[^a-zA-Z0-9_]"|awk '{print $3}'|grep -v FIXMECROSSSECTION|uniq))
	    CROSSSECTIONS=($(echo "${CROSSSECTIONS[@]}"|tr ' ' '\n'|sort -nu))
	    if [ ${#CROSSSECTIONS[@]} -eq 1 ];
	    then 
		CROSSSECTION=${CROSSSECTIONS[0]}
	    else
		echo "candidate cross sections= ${CROSSSECTIONS[@]}"
		read -p "select cross section: " CROSSSECTION
		[ -z "$CROSSSECTION" ] && CROSSSECTION=FIXMECROSSSECTION
	    fi
	fi

	NUMS=($(echo 'cout<<Form("%d\t%.1f\t%.12e",(int)sumW->GetEntries(),sumSign->GetSum(),sumW->GetSum())<<endl;'|root -b -l ${line}|tail -n1))

	## for CommonSampleInfo file
	OUT=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/CommonSampleInfo/${NAME}.txt
	NEWLINE="$NAME\t$DASNAME\t$CROSSSECTION\t${NUMS[0]}\t${NUMS[1]}\t${NUMS[2]}"
	if [ -f "$OUT" ];then
	    DIFF=$(diff "$OUT" <(echo -e "# alias PD xsec nmc sumsign sumw\n$NEWLINE"))
	    if [ -n "$DIFF" ];then
		echo "$DIFF"
		echo " > Update $OUT with $line"
		echo -e "$NEWLINE"
		read -p "(y/n): " YES
		if [ "$YES" = "y" ];then
		    echo -e "# alias PD xsec nmc sumsign sumw\n$NEWLINE" > $OUT
		fi
	    fi
	else
	    echo " > Update $OUT with $line"
	    echo -e "$NEWLINE"
	    read -p "(y/n): " YES
	    if [ "$YES" = "y" ];then
		echo -e "# alias PD xsec nmc sumsign sumw\n$NEWLINE" > $OUT
	    fi
	fi

	## for summary file
	SUMMARYFILE=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/SampleSummary_MC.txt
	EXISTSUMMARY=$(grep "^$NAME[^a-zA-Z0-9_]*$DASSNAME[^a-zA-Z0-9_].*$" $SUMMARYFILE)
	if [ -z "$EXISTSUMMARY" ]
	then
	    echo "SampleSummary_MC.txt: add $NEWLINE"
	    echo -e "$NEWLINE" >> $SUMMARYFILE
	elif [ "$EXISTSUMMARY" != "$(echo -e "$NEWLINE")" ];then
	    echo -e "SampleSummary_MC.txt: $EXISTSUMMARY --> $NEWLINE"
	    read -p "(y/n): " YES
	    if [ "$YES" = "y" ];then
		sed -i "s/^$NAME[^a-zA-Z0-9_]*$DASSNAME[^a-zA-Z0-9_].*$/$NEWLINE/" $SUMMARYFILE
	    fi
	fi
    fi
done 3< <(find $SKFlatOutputDir$SKFlatV/GetEffLumi -type f|sort -V)



#find data/Run2Legacy_v3/*/Sample/CommonSampleInfo -type f -mtime -1|while read line;do array=(${line//\// });NAME=${array[5]/.txt/}; DASNAME=$(head -n1 ${line/CommonSampleInfo/ForSNU/}|awk 'BEGIN{FS="/"}{print $9}'); NUMS=($(echo 'cout<<sumW->GetSum()<<"\t"<<sumW->GetEntries()<<"\t"'|root -b -l $SKFlatOutputDir/$SKFlatV/GetEffLumi/${array[2]}/GetEffLumi_${array[5]/.txt/.root}|tail -n1|awk '{print $1"\t"$2}')); echo $NAME $DASNAME CROSSSECTION ${NUMS[@]} > $line;done
