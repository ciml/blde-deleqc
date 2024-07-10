import bz2
from statistics import median

# Specify the path to the bz2 file
file_path = "/home/hedersb/workspace/blde-deleqc/blde+deleqc-iii/results/"

# static parameters
tol=0.0001
eps=0.0001
popL=30
popF=30
genL=200

cr=0.9
varL=3
f=0.7

for f in (0.6, 0.7, 0.8, 0.9, 1):

    for cr in (0.5, 0.6, 0.7, 0.8, 0.9):

        for varL in (1, 2, 3, 4):

            for fobj in range(20, 24+1):
                
                for genF in (50, 100, 200, 500):

                    results = []

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
                                    if (constraintViolation < 0.0001):
                                        results.append(result)

                    #print(f"&   {fobj}  & {genF} & {sum(results)/len(results):.2f}  & {median(results):.2f} & {round(sum(objectiveFunctionEvaluations)/len(objectiveFunctionEvaluations))} & {round(median(objectiveFunctionEvaluations))} & {len(results)} \\\\")
                    # write file
                    with open(file_path + f"{f}-{cr}-{varL}.dat", "a") as file:
                        if len(results) > 0:
                            file.write(f"{fobj}-{genF}\t{median(results):.16f}\n")