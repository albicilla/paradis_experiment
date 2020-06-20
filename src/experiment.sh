for i in `seq 0 7`; do
    for j in `seq 1 5`; do
	((thnum = 2 ** i))
	#echo ${i}
	echo ${thnum}
	./paradis_ompf_repair_read ${thnum} ./1_1_10^9 1000000000
    done
done
