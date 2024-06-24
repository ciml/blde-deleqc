

mkdir results
i=1

for seed in {1..10}; do
	for fl in 0.6 0.7 0.8; do
		for genL in 50 100 200 500; do
		
			qsub run.sh $seed $fl $genL $i
			
			i=$(($i+1))
			
		done
	done
done

