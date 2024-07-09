#$ -S /bin/bash
#$ -q linux.q
#----------------------------------------------------------
# Parameters:
# 1- 
# 2- 
# 3- 
# 4- 
#----------------------------------------------------------

# static parameters
cr=0.9
#tol=0.0001
#eps=0.0001
tol=0
eps=0

p=1
neval=1000

varL=3
popL=30
popF=30
f=0.7
genL=200


for seed in {1..10}; do
	
	for genF in 50 100 200 500; do

		# blde+deleqc-iii
		echo "blde+deleqc-iii"
		g++ -O2 main.cpp funcoes.cpp LU.c -o cilamce2024
		
		for fobj in 20 21 22 23 24; do
			echo "./cilamce2024 -genL $genL -popL $popL -popF $popF -genF $genF -F $f -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $seed | bzip2 -9 -c > results/cilamce2024_blde-deleqc-iii_${genL}_${popL}_${genF}_${popF}_${f}_${cr}_${varL}_${fobj}_${tol}_${eps}_${seed}.out.bz2"
			./cilamce2024 -genL $genL -popL $popL -popF $popF -genF $genF -F $f -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $seed | bzip2 -9 -c > results/cilamce2024_blde-deleqc-iii_${genL}_${popL}_${genF}_${popF}_${f}_${cr}_${varL}_${fobj}_${tol}_${eps}_${seed}.out.bz2
		done

	done

done


