/* 
 * File:   main.cpp
 * Author: Jaque
 *
 * Created on 26 de Abril de 2012, 15:06
 */

#include <cstdlib>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

int fim = 0;
int nobj = 2;
const int totalSolucoes = 1000000;

struct naodominada{
   double obj[2];
   char dominada;  /* 0 = nao dominada, 1 = dominada */
}ndom[totalSolucoes], objs;


/*####################################################
              VERIFICAR DOMINANCIA
####################################################*/
int checar_dominancia (struct naodominada a, struct naodominada b ){
    int i, flag1, flag2;
    flag1 = 0;
    flag2 = 0;

    for (i = 0; i < nobj ; i++ ){
        if (a.obj[i] < b.obj[i])
           flag1 = 1;
        else {
            if (a.obj[i] > b.obj[i])
               flag2 = 1;
        }
    }

    if (flag1==1 && flag2==0) {
       return (1);        /* a domina b */
    }
    else {
         if (flag1==0 && flag2==1) {
            return (-1);  /* d domina a */
         }
         else{
            return (0);   /* nao dominados */
         }
    }
}

/*####################################################
        INCLUIR ELEMENTO NO CONJ. DE PARETO
####################################################*/
void incluir_conj_pareto (struct naodominada op ){

     for (int i=0 ; i<nobj ; i++){
          ndom[fim].obj[i] = op.obj[i];
     }
     fim = fim + 1;

}


/*####################################################
                        PRINCIPAL
####################################################*/
int main( int argc, char** argv ){
    FILE *entrada, *saida;

	cout << "Entrada: " << argv[1] << "\nSaida: " << argv[2] << endl;

    entrada = fopen(argv[1], "r"); ;
    saida = fopen( argv[2] , "w+"); ;

    int ponteiro = 0;
    int val, pont;

    /* INICIAR FLAG PARA NAO DOMINADA */
    for (int i=0 ; i < totalSolucoes ;i++)
        ndom[i].dominada = 0; // NAO DOMINADA

    while ( !feof(entrada) ) {

        //fscanf(entrada, "%lf %lf", &objs.obj[0], &objs.obj[1]);
        fscanf(entrada, "%lf\t%lf", &objs.obj[0], &objs.obj[1]);
        incluir_conj_pareto(objs);

    }

    /* VERIFICAR DOMINANCIA */
    for (ponteiro = 0; ponteiro < fim; ponteiro++){
        for (pont = 0; pont < fim; pont++){
            if (ndom[ponteiro].dominada != 1){
               val = checar_dominancia (ndom[pont], ndom[ponteiro]);
               switch (val){
                  case 1: // a domina b
                       ndom[ponteiro].dominada = 1;
                       pont = fim;
                  break;
               }
            } else {
                   pont = fim;
            }
        }
    }

    printf("\n NAO DOMINADAS FINAL \n");
    for (int k = 0; k < fim ; k++){
        if (ndom[k].dominada == 0){
  //         printf("\n NDOM[%d]\t OBJ1:\t%f \t OBJ2:\t%f \tDOM:%d", k, ndom[k].obj[0], ndom[k].obj[1], ndom[k].dominada);
           fprintf(saida, "%f\t %f\n", ndom[k].obj[0], ndom[k].obj[1]);
        }
    }

    fclose (entrada);
    fclose (saida);
//    getchar();

}
