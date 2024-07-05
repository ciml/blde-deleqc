#include "math.h"
#include "stdlib.h"

void calculaFuncao(double *ind, int d, int nivel, double *leader, double *follower, int funcao);
int getDimensao(int funcao, int nivel);
int getTipo(int funcao, int nivel);
double getLower(int nivel, int funcao, int indice);
double getUpper(int nivel, int funcao, int indice);
int getNumberOfConstraints(int funcao, int nivel);

int getNEval(int nivel);

extern double EPS;
