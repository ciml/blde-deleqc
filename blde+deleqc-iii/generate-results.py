import bz2
from statistics import median

# Specify the path to the bz2 file
file_path = "/home/hedersb/workspace/blde-deleqc/blde+deleqc-iii/results/"

# static parameters
cr=0.9
tol=0.0001
eps=0.0001
#tol=0
#eps=0

varL=3
popL=30
popF=30
f=0.7
genL=200

for fobj in range(20, 24+1):
	
    for genF in (50, 100, 200, 500):

        results = []
        objectiveFunctionEvaluations = []

        for seed in range(1, 10+1):

            # Open the bz2 file in read mode
            filename = f"cilamce2024_blde-deleqc-iii_{genL}_{popL}_{genF}_{popF}_{f}_{cr}_{varL}_{fobj}_{tol}_{eps}_{seed}.out.bz2"
            with bz2.open(file_path + filename, "rt") as file:
                # Read lines from the file
                for line in file:
                    # Process each line
                    eachLine = line.strip()
                    if eachLine.startswith(f"G-{genL-1}"):
                        lineElements = eachLine.split(" ")
                        i = 0
                        while i < len(lineElements) and lineElements[i] != "Fit:":
                            i += 1
                        i += 1
                        #print(lineElements[i])
                        result = float(lineElements[i])
                        while i < len(lineElements) and lineElements[i] != "Const:":
                            i += 1
                        i += 1
                        constraintViolation = float(lineElements[i])
                        #print(lineElements[i])
                        #print(i)
                        objectiveFunctionEvaluationsInThisRun = int(lineElements[len(lineElements)-1])
                        if (constraintViolation < 0.0001):
                            results.append(result)
                            objectiveFunctionEvaluations.append(objectiveFunctionEvaluationsInThisRun)

        #print(f"&   {fobj}  & {genF} & {sum(results)/len(results):.2f}  & {median(results):.2f} & {round(sum(objectiveFunctionEvaluations)/len(objectiveFunctionEvaluations))} & {round(median(objectiveFunctionEvaluations))} & {len(results)} \\\\")
        if len(results) > 0:
            print(f"& {genF} & {sum(results)/len(results):.2f}  & {median(results):.2f} & {round(sum(objectiveFunctionEvaluations)/len(objectiveFunctionEvaluations))} & {round(median(objectiveFunctionEvaluations))} & {len(results)} \\\\")
        else:
            print(f"& {genF} & -- & -- & -- & -- & -- \\\\")