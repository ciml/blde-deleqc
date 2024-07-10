list_of = "algorithms" #"problems"


# static parameters
tol=0.0001
eps=0.0001
popL=30
popF=30
genL=200

cr=0.9
varL=3
f=0.7

if list_of == "problems":

    for fobj in range(20, 24+1):
        
        for genF in (50, 100, 200, 500):

            print(f"{fobj}-{genF}")

else:
    for f in (0.6, 0.7, 0.8, 0.9, 1):

        for cr in (0.5, 0.6, 0.7, 0.8, 0.9):

            for varL in (1, 2, 3, 4):

                print(f"{f}-{cr}-{varL}")