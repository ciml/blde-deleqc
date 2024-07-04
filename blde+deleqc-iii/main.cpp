#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include <iostream>
using namespace std;

#define PENALTY 1

int SIZEL, SIZEF; // tamanho da população (L = leader; F = follower)
int GENL, GENF;   // número de gerações
int SEED;         // semente
int FUNCAO;       // função teste
int DIML, DIMF;   // dimensão do problema (variáveis X-Lider e Y-Seguidor)
double F, CR;     // parâmetros do DE
int VARIANTE;     // variante do DE
double TOLERANCE; // tolerancia para o criterio de parada
double EPS;       // tolerancia da restricao

#include "funcoes.h"

double **popL;
double **popLNova;
// Populacao para armazenar os valores do Follower correspondentes ao Leader em popL
double **popLValoresF;
double residuo;
double epsilon = pow(10, -4);

void inicializaFollower(double **&pop, double *leader, int n, int d);
void inicializa(double **&pop, int n, int d, int nivel = 2);

//=============================================================================
//=============================================================================
void calculaVariancia(double **pop, double *var, int dim, int size)
{
    // variancias e medias de cada variavel
    double *med = new double[dim];

    for (int d = 0; d < dim; d++)
    {
        double soma = 0;
        for (int n = 0; n < size; n++)
        {
            soma += pop[n][d];
        }
        med[d] = soma / size;
    }

    for (int d = 0; d < dim; d++)
    {
        double soma_Pvar = 0;
        for (int n = 0; n < size; n++)
        {
            soma_Pvar += (pop[n][d] - med[d]) * (pop[n][d] - med[d]);
        }
        var[d] = soma_Pvar / size;
    }

    delete[] med;
}

//==========================================================
void calculaAptidao(double *ind, int d, int nivel, double *leader, double *follower)
{

    if (follower == NULL)
    {
        ind[d] = RAND_MAX;
        ind[d + 1] = RAND_MAX;
        ind[d + 2] = RAND_MAX;
        return;
    }

    calculaFuncao(ind, d, nivel, leader, follower, FUNCAO);
}

void imprimePopulacao(double **pop, int n, int d)
{
    for (int i = 0; i < n; i++)
    {
        cout << i << ") ";
        for (int j = 0; j < d; j++)
        {
            cout << pop[i][j] << " ";
        }
        cout << " Fit: " << pop[i][d] << " Const: " << pop[i][d + 1];

        cout << " Foll.: ";
        for (int j = 0; j < DIMF; j++)
        {
            cout << popLValoresF[i][j] << " ";
        }
        cout << " Fit: " << popLValoresF[i][DIMF] << " Const: " << popLValoresF[i][DIMF + 1] << endl;
    }
}

void selecionaIndividuos(int &ind1, int &ind2, int &ind3, int i, int n)
{
    do
    {
        ind1 = rand() % n;
    } while (ind1 == i);

    do
    {
        ind2 = rand() % n;
    } while (ind2 == i || ind2 == ind1);

    do
    {
        ind3 = rand() % n;
    } while (ind3 == i || ind3 == ind1 || ind3 == ind2);
}

int compara(double *ind1, double *ind2, int d, int nivel)
{

    // Critério de seleção do Deb para tratamento de restrições
    // Esta parte considera o caso que alguma restrição foi violado
    if (nivel == 1)
    {
        if (ind1[d + 1] <= EPS && ind2[d + 1] > EPS)
        {
            return 1; // ind1 nao violou
        }
        else if (ind1[d + 1] > EPS && ind2[d + 1] <= EPS)
        {
            return 0; // ind2 nao violou
        }
        else if (ind1[d + 1] > EPS && ind2[d + 1] > EPS)
        { // ambos violam as rest.
            if (ind1[d + 1] <= ind2[d + 1])
            { // verifica quem viola menos
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }
    else if (nivel == 2)
    {
        if (ind1[d + 2] <= 0 && ind2[d + 2] > 0)
        {
            return 1; // ind1 nao violou
        }
        else if (ind1[d + 2] > 0 && ind2[d + 2] <= 0)
        {
            return 0; // ind2 nao violou
        }
        else if (ind1[d + 2] > 0 && ind2[d + 2] > 0)
        { // ambos violam as rest.
            if (ind1[d + 2] <= ind2[d + 2])
            { // verifica quem viola menos
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }

    // Esta parte considera o caso que nenhuma restrição foi violado;
    // neste caso o procedimento retorna o de melhor aptidao
    //-----------------------------------
    // Se a função é de maximicao compara com >=
    // Caso contrário usa <=
    if (getTipo(FUNCAO, nivel) == 1)
    {
        if (ind1[d] <= ind2[d])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (ind1[d] >= ind2[d])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

int selecionaMelhor(double *ind, double **pop, int n, int d, int nivel)
{
    /*
        int m = 0;

        // Se for leader
        if (nivel == 1){
          for(int j = 0; j < d + 3; j++){
                ind[j] = pop[0][j];
          }

          for (int i = 1; i < n; i++){
                if( compara(pop[i], ind, d, nivel) > 0 ){
                    for(int j = 0; j < d + 3; j++){ //======!!!!!!!!!!!!
                        ind[j] = pop[i][j];
                        m = i;
                    }
                }
          }
          return m;
        } else {  // se for follower
                  for(int j = 0; j < d + 3; j++){
                      ind[j] = pop[0][j];
                   }

                  for (int i = 1; i < n; i++){
                     if( compara(pop[i], ind, d, nivel) > 0 ){
                        for(int j = 0; j < d + 3; j++){
                           ind[j] = pop[i][j];
                          m = i;
                        }
                     }
                  }
                return m;
               }
      * */

    int m = 0;
    for (int j = 0; j < d + 3; j++)
    {
        ind[j] = pop[0][j];
    }

    for (int i = 1; i < n; i++)
    {
        if (compara(pop[i], ind, d, nivel) > 0)
        {
            for (int j = 0; j < d + 3; j++)
            {
                ind[j] = pop[i][j];
                m = i;
            }
        }
    }
    return m;
}

int iguais(double *ind1, double *ind2, int d)
{

    for (int i = 0; i < d; i++)
    {
        if (ind1[i] != ind2[i])
        {
            return 0;
        }
    }
    return 1;
}

// ===== Algoritmo do Seguidor =====
void deFollower(double *uL, double *uF)
{

    // uL é a variável que vem do Lider (variável X);
    // Neste código vamos obter o Y associado ao X que veio do Líder.

    double **popF;
    double **popFNova;
    inicializaFollower(popF, uL, SIZEF, DIMF);
    inicializa(popFNova, SIZEF, DIMF, 2);

    // variancias e medias de cada variavel
    double *var_inicial = new double[DIMF];
    double *var_atual = new double[DIMF];

    calculaVariancia(popF, var_inicial, DIMF, SIZEF);

    for (int gF = 0; gF < GENF; gF++)
    {

        selecionaMelhor(uF, popF, SIZEF, DIMF, 2);
        // uF = melhor da populacao do seguidor; utilizado somente na
        // variante DE/best/1/bin

        for (int i = 0; i < SIZEL; i++)
        {

            int ind1, ind2, ind3;
            selecionaIndividuos(ind1, ind2, ind3, i, SIZEL);

            double *u = new double[DIMF + 3];
            //            int jRand = rand()%DIMF;
            for (int j = 0; j < DIMF; j++)
            {
                //                if (j == jRand || rand()/(float)RAND_MAX < CR){

                if (VARIANTE == 1)
                {
                    // DE/rand/1/bin
                    u[j] = popF[ind1][j] + F * (popF[ind2][j] - popF[ind3][j]);
                }
                else if (VARIANTE == 2)
                {
                    // DE/best/1/bin
                    u[j] = uF[j] + F * (popF[ind2][j] - popF[ind3][j]);
                }
                else if (VARIANTE == 3)
                {
                    // DE/target-to-rand/1/bin
                    u[j] = popF[i][j] + F * (popF[ind1][j] - popF[i][j]) + F * (popF[ind2][j] - popF[ind3][j]);
                }
                else if (VARIANTE == 4)
                {
                    // DE/target-to-best/1/bin
                    u[j] = popF[i][j] + F * (uF[j] - popF[i][j]) + F * (popF[ind2][j] - popF[ind3][j]);
                }

                //cout << " " << u[j];
            }
            //cout << endl;

            calculaAptidao(u, DIMF, 2, uL, u);

            if (compara(u, popF[i], DIMF, 2) > 0)
            {
                for (int j = 0; j < DIMF + 3; j++)
                {
                    popFNova[i][j] = u[j];
                }
            }
            else
            {
                for (int j = 0; j < DIMF + 3; j++)
                {
                    popFNova[i][j] = popF[i][j];
                }
            }
            delete[] u;
        }

        // copia a populacao nova
        for (int i = 0; i < SIZEF; i++)
        {
            for (int j = 0; j < DIMF + 3; j++)
            {
                popF[i][j] = popFNova[i][j];
            }
        }

        // Criterio de parada Sinha&Deb
        calculaVariancia(popF, var_atual, DIMF, SIZEF);
        double soma_total = 0;
        for (int d = 0; d < DIMF; d++)
        {
            soma_total += var_atual[d] / var_inicial[d];
        }
        //        cout << endl << "alfa: " << soma_total;
        if (soma_total < TOLERANCE)
        {
            break;
        }
        //---------------------------------------------
    }

    selecionaMelhor(uF, popF, SIZEF, DIMF, 2);

    //    cout << endl << " BEST: " ;
    //    for (int i = 0 ; i < DIMF+3 ; i++ )
    //        cout << " " << uF[i];
    //    getchar();
    // uF = melhor da populacao do seguidor (variável Y),
    // que irá retornar para o Líder.

    for (int i = 0; i < SIZEF; i++)
    {
        delete[] popF[i];
        delete[] popFNova[i];
    }
    delete[] popF;
    delete[] popFNova;
}

void deLeader()
{

    /////////////////////////////////////////////
    double *uBest = new double[DIML + 3];
    selecionaMelhor(uBest, popL, SIZEL, DIML, 1);
    /////////////////////////////////////////////

    for (int i = 0; i < SIZEL; i++)
    {

        int ind1, ind2, ind3;
        selecionaIndividuos(ind1, ind2, ind3, i, SIZEL);

        double *u = new double[DIML + 3];
        int jRand = rand() % DIML;
        for (int j = 0; j < DIML; j++)
        {
            if (j == jRand || rand() / (float)RAND_MAX < CR)
            {

                if (VARIANTE == 1)
                {
                    //---------- DE/rand/1/bin
                    u[j] = popL[ind1][j] + F * (popL[ind2][j] - popL[ind3][j]);
                }
                else if (VARIANTE == 2)
                {
                    //--------- DE/best/1/bin
                    u[j] = uBest[j] + F * (popL[ind2][j] - popL[ind3][j]);
                }
                else if (VARIANTE = 3)
                {
                    //--------- DE/target-to-rand/1/bin
                    u[j] = popL[i][j] + F * (popL[ind1][j] - popL[i][j]) + F * (popL[ind2][j] - popL[ind3][j]);
                }
                else if (VARIANTE = 4)
                {
                    //--------- DE/target-to-best/1/bin
                    u[j] = popL[i][j] + F * (uBest[j] - popL[i][j]) + F * (popL[ind2][j] - popL[ind3][j]);
                }
            }
            else
            {
                u[j] = popL[i][j];
            }
            // verifica se ultrapassou o limite superior ou inferior das variáveis
            if (u[j] < getLower(1, FUNCAO, j)) // LOWER
                u[j] = getLower(1, FUNCAO, j);
            else if (u[j] > getUpper(1, FUNCAO, j)) // UPPER
                u[j] = getUpper(1, FUNCAO, j);
        }

        double *uF = new double[DIMF + 3];
        // Para cada X do líder, executamos o código deFollower para
        // obter a variável Y do seguidor, associado ao X

        //        cout << "  " << u[0] << "  " << u[1] << endl;

        deFollower(u, uF);

        // se a solução do líder atual (u) for igual a anterior (popL) e
        // a solução do seguidor atual (uF) for diferente da anteior (popLValoresF).
        if ((iguais(u, popL[i], DIML) == 1) && (iguais(uF, popLValoresF[i], DIMF) == 0))
        {

            // verifica se a solução atual (uF) do seguidor  é melhor
            // que a anterior (popLValoresF)
            // compara(...) -> critério de selecao do Deb
            if (compara(uF, popLValoresF[i], DIMF, 2) > 0)
            { // uF melhor

                for (int j = 0; j < DIMF + 3; j++)
                {
                    popLValoresF[i][j] = uF[j]; // atualiza solução com valor de uF
                }
                calculaAptidao(popL[i], DIML, 1, popL[i], uF); // calcula nova aptidao para popL

                // verifica se seguidor viola restrição
                if (uF[DIMF + 1] > 0 || uF[DIMF + 2] > 0)
                {
                    // acumula valor da restricao do seguidor
                    // nas restricoes do lider em popL
                    popL[i][DIML + 1] = popL[i][DIML + 1] + PENALTY * uF[DIMF + 1];
                    popL[i][DIML + 1] = popL[i][DIML + 1] + PENALTY * uF[DIMF + 2];
                }
            }
            else
            {
                for (int j = 0; j < DIMF + 3; j++)
                {
                    uF[j] = popLValoresF[i][j]; // atualiza solução com valor de popLValoresF
                }
            }
        }
        // calcula aptidao de u
        calculaAptidao(u, DIML, 1, u, uF);

        // verifica se seguidor viola restrição
        if (uF[DIMF + 1] > 0 || uF[DIMF + 2] > 0)
        {
            // acumula valor da restricao do seguidor
            // nas restricoes do lider em u
            u[DIML + 1] = u[DIML + 1] + PENALTY * uF[DIMF + 1];
            u[DIML + 1] = u[DIML + 1] + PENALTY * uF[DIMF + 2];
        }

        // compara(...) -> critério de selecao do Deb
        if (compara(u, popL[i], DIML, 1) > 0)
        {

            for (int j = 0; j < DIML + 3; j++)
            {
                popLNova[i][j] = u[j];
            }
            for (int j = 0; j < DIMF + 3; j++)
            {
                popLValoresF[i][j] = uF[j];
            }
        }
        else
        {
            for (int j = 0; j < DIML + 3; j++)
            {
                popLNova[i][j] = popL[i][j];
            }
        }
        delete[] uF;
        delete[] u;
    }

    // copia a populacao nova
    for (int i = 0; i < SIZEL; i++)
    {
        for (int j = 0; j < DIML + 3; j++)
        {
            popL[i][j] = popLNova[i][j];
        }
    }

    delete[] uBest;
}

void imprimeCabecalho()
{
    cout << "g leader ";
    for (int i = 0; i < DIML; i++)
    {
        cout << "x" << i << " ";
    }
    cout << "fitLeader fitLeaderValue constLeader constLeaderValue follower ";
    for (int i = 0; i < DIMF; i++)
    {
        cout << "y" << i << " ";
    }
    cout << "fitFollower fitFollowerValue constFollower constFollowerValue" << " nEvalL nEvalF" << endl;
}

// void  MatrizE(double **E, double *c, int restigual,double *leader)
//{

//}
//*********************************************************************
// Decompõe a matriz M = LU

void decomposicao_LU(double **M, int *p, unsigned int n)
{

    double max; // atual valor do pivo
    int imax;   // atual valor do índice do pivo
    int temp;   // variavel auxiliar
    double m;   // fator multiplicativo

    int k, i, j;
    for (k = 0; k < n - 1; k++)
    {
        // Encontrar o pivo--------
        max = fabs(M[p[k]][k]);
        imax = k;
        for (i = k + 1; i < n; i++)
        {
            if (max < fabs(M[p[i]][k]))
            {
                max = fabs(M[p[i]][k]);
                imax = i;
            }
        }
        //------------------------

        // Troca 'vitual das linhas da matriz A'---
        temp = p[k];
        p[k] = p[imax];
        p[imax] = temp;
        //----------------------------------------

        // Anular elementos--------
        for (i = k + 1; i < n; i++)
        {
            m = M[p[i]][k] / M[p[k]][k];
            for (j = k; j < n; j++)
            {
                M[p[i]][j] = M[p[i]][j] - m * M[p[k]][j];
            }
            M[p[i]][k] = m;
        }
        //------------------------
    }
}

// Resolve o sistema  My= c

double *resolveSistema(double **M, int *p, double *c, unsigned int n)
{

    double *aux2 = (double *)malloc(n * sizeof(double)); // variável para amazenar valores parciais da resposta final
    double *resF = (double *)malloc(n * sizeof(double)); // variável para amazenar valores finais

    int k, j;
    // Substituição progressiva-------------------------
    for (k = 0; k < n; k++)
    {
        aux2[k] = c[p[k]];
        for (j = 0; j < k; j++)
        {
            aux2[k] = aux2[k] - M[p[k]][j] * aux2[j];
        }
        aux2[k] = aux2[k] / 1.0;
    }
    //-------------------------------------------------

    // Substituição regressiva--------------------------
    // obs: como usamos tipo 'unsigned int', para que
    // o comando 'for' funcionasse corretamente, foi necessário
    // uma adaptação na variável k (-1 == 1)
    for (k = n; k > 0; k--)
    {
        resF[k - 1] = aux2[k - 1];
        for (j = k; j < n; j++)
        {
            resF[k - 1] = resF[k - 1] - M[p[k - 1]][j] * resF[j];
        }
        resF[k - 1] = resF[k - 1] / M[p[k - 1]][k - 1];
    }
    //-------------------------------------------------

    free(aux2);

    return resF;
}

void Funcao_Ajuste_1(double **populacao, int restigual, int dim, int n, double **E, double *c, double **Residuo)
{
    int i, j, k;

    for (k = 0; k < n; k++)
    {
        for (i = 0; i < restigual; i++)
        {
            Residuo[k][i] = 0.0;
            for (j = 0; j < dim; j++)
            {
                Residuo[k][i] = Residuo[k][i] + E[i][j] * populacao[k][j];
            }
            Residuo[k][i] = Residuo[k][i] - c[i];
        }
    }
}

void Funcao_Ajuste_2(double **M, double **Residuo, int restigual, int *p, double *c, double *y, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        y = resolveSistema(M, p, Residuo[i], restigual); // Encontra o resultado do sistema Mp=r
    }
}

double *Funcao_Ajuste_3(double *y, double **Et, double *id, int restigual, int dim)
{
    int i, j;
    double *vetor_correto;

    vetor_correto = (double *)malloc(dim * sizeof(double));

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < restigual; j++)
        {
            vetor_correto[i] = id[i] - Et[i][j] * y[i];
        }
    }
    return vetor_correto;
}

void Calcula_Residuo(int n, double *id, double **populacao, int restigual, int dim, double **E, double *c, double **Residuo, double **M, int *p, double *y, double **Et, double *vetor_correto)
{
    int i, j, k;
    double *maior;

    maior = (double *)malloc((n) * sizeof(double));

    Funcao_Ajuste_1(populacao, restigual, dim, n, E, c, Residuo);

    // residuo = 0.0;
    for (i = 0; i < n; i++)
    {
        maior[i] = 0.0;
        for (j = 0; j < restigual; j++)
        {
            if (maior[i] < fabs(Residuo[i][j]))
            {
                maior[i] = fabs(Residuo[i][j]);
                residuo = residuo + Residuo[i][j];
            }
            if (maior[i] > epsilon)
            {
                Funcao_Ajuste_2(M, Residuo, restigual, p, c, y, n);
                for (k = 0; k < dim; k++)
                {
                    vetor_correto = Funcao_Ajuste_3(y, Et, id, restigual, dim);
                    populacao[i][k] = vetor_correto[k];
                }
            }
        }
    }
    free(maior);
}

// Inicializa a população do Seguidor
void inicializaFollower(double **&pop, double *leader, int n, int dim)
{

    pop = new double *[n];
    for (int i = 0; i < n; i++)
    {
        pop[i] = new double[dim + 3];
    }

    //    cout << " " << leader[0] << " " << leader[1] << endl;

    int restigual = 4;
    double *c, *vetor_correto;
    double *u, *y, *z;
    double **E, **Et, **M, **Residuo;
    double x0[dim], v[dim], x[dim];
    int *p;
    int i, j, k; // contadores
    double d[dim];

    p = (int *)malloc(restigual * sizeof(int));
    y = (double *)malloc(restigual * sizeof(double));
    z = (double *)malloc(restigual * sizeof(double));
    E = (double **)malloc(restigual * sizeof(double *));
    Et = (double **)malloc(dim * sizeof(double *));
    M = (double **)malloc(restigual * sizeof(double));
    c = (double *)malloc(restigual * sizeof(double *));
    Residuo = (double **)malloc(n * sizeof(double *));
    vetor_correto = (double *)malloc(dim * sizeof(double));

    for (i = 0; i < restigual; i++)
    {
        M[i] = (double *)malloc(restigual * sizeof(double));
        E[i] = (double *)malloc(dim * sizeof(double));
    }

    for (i = 0; i < dim; i++)
    {
        Et[i] = (double *)malloc(restigual * sizeof(double));
    }

    for (j = 0; j < n; j++)
    {
        Residuo[j] = (double *)malloc(restigual * sizeof(double));
    }

    // Pega a matriz E, o vetor c, a quantidade de restrições de igualdade

    if (FUNCAO == 20)
    {

        restigual = 3;

        E[0][0] = -1.0;
        E[0][1] = 1.0;
        E[0][2] = 1.0;
        E[0][3] = 1.0;
        E[0][4] = 0.0;
        E[0][5] = 0.0;
        E[1][0] = -1.0;
        E[1][1] = 2.0;
        E[1][2] = -0.5;
        E[1][3] = 0.0;
        E[1][4] = 1.0;
        E[1][5] = 0.0;
        E[2][0] = 2.0;
        E[2][1] = -1.0;
        E[2][2] = -0.5;
        E[2][3] = 0.0;
        E[2][4] = 0.0;
        E[2][5] = 1.0;

        c[0] = 1.0;
        c[1] = 1.0 - 2.0 * leader[0];
        c[2] = 1.0 - 2.0 * leader[1];
    }
    else if (FUNCAO == 21)
    {

        restigual = 3;

        E[0][0] = 1.0;
        E[0][1] = -1.0;
        E[0][2] = 1.0;
        E[0][3] = 0.0;
        E[0][4] = 0.0;
        E[1][0] = 0.0;
        E[1][1] = 1.0;
        E[1][2] = 0.0;
        E[1][3] = 1.0;
        E[1][4] = 0.0;
        E[2][0] = 0.0;
        E[2][1] = 0.0;
        E[2][2] = 0.0;
        E[2][3] = 0.0;
        E[2][4] = 1.0;

        c[0] = 2 * leader[0] - 2.5;
        c[1] = 2 - leader[0] + 3 * leader[1];
        c[2] = 2 - leader[0] - leader[1];
    }
    else if (FUNCAO == 22)
    {

        restigual = 4;

        E[0][0] = 5.0;
        E[0][1] = 4.0;
        E[0][2] = 1.0;
        E[0][3] = 0.0;
        E[0][4] = 0.0;
        E[0][5] = 0.0;
        E[1][0] = -5.0;
        E[1][1] = 4.0;
        E[1][2] = 0.0;
        E[1][3] = 1.0;
        E[1][4] = 0.0;
        E[1][5] = 0.0;
        E[2][0] = -4.0;
        E[2][1] = 5.0;
        E[2][2] = 0.0;
        E[2][3] = 0.0;
        E[2][4] = 1.0;
        E[2][5] = 0.0;
        E[3][0] = 4.0;
        E[3][1] = 5.0;
        E[3][2] = 0.0;
        E[3][3] = 0.0;
        E[3][4] = 0.0;
        E[3][5] = 1.0;

        c[0] = 12 - 4 * leader[0];
        c[1] = -4 + 4 * leader[0];
        c[2] = 4 - 4 * leader[0];
        c[3] = 4 + 4 * leader[0];
    }
    else if (FUNCAO == 23)
    {

        restigual = 3;

        E[0][0] = -1.0;
        E[0][1] = 1.0;
        E[0][2] = 1.0;
        E[0][3] = 1.0;
        E[0][4] = 0.0;
        E[0][5] = 0.0;
        E[1][0] = -1.0;
        E[1][1] = 2.0;
        E[1][2] = -0.5;
        E[1][3] = 0.0;
        E[1][4] = 1.0;
        E[1][5] = 0.0;
        E[2][0] = 2.0;
        E[2][1] = -1.0;
        E[2][2] = -0.5;
        E[2][3] = 0.0;
        E[2][4] = 0.0;
        E[2][5] = 1.0;

        c[0] = 1;
        c[1] = 1 - 2 * leader[0];
        c[2] = 1 - 2 * leader[1];
    }
    else if (FUNCAO == 24)
    {

        restigual = 3;

        E[0][0] = -1.0;
        E[0][1] = 1.0;
        E[0][2] = 1.0;
        E[0][3] = 1.0;
        E[0][4] = 0.0;
        E[0][5] = 0.0;
        E[1][0] = -1.0;
        E[1][1] = 2.0;
        E[1][2] = -0.5;
        E[1][3] = 0.0;
        E[1][4] = 1.0;
        E[1][5] = 0.0;
        E[2][0] = 2.0;
        E[2][1] = -1.0;
        E[2][2] = -0.5;
        E[2][3] = 0.0;
        E[2][4] = 0.0;
        E[2][5] = 1.0;

        c[0] = 1;
        c[1] = 1 - 2 * leader[0];
        c[2] = 1 - 2 * leader[1];
    }

    // Calcula a matriz transposta de E
    for (i = 0; i < restigual; i++)
    {
        for (j = 0; j < dim; j++)
        {
            Et[j][i] = E[i][j];
        }
    }

    // Calcula o produto de E por Et e encontra a matriz M
    for (i = 0; i < restigual; i++)
    {
        for (j = 0; j < restigual; j++)
        {
            M[i][j] = 0.0;
            for (k = 0; k < dim; k++)
            {
                M[i][j] = M[i][j] + E[i][k] * Et[k][j];
            }
        }
    }

    // Inicializa Pivo
    for (i = 0; i < restigual; i++)
    {
        p[i] = i;
    }

    // Faz a decomposição da Matriz M = LU
    decomposicao_LU(M, &p[0], restigual);   // Encontra as matrizes 'L' e 'U'
    y = resolveSistema(M, p, c, restigual); // Encontra o resultado do sistema Mx=b

    // Calcula o valor de x0
    for (i = 0; i < dim; i++)
    {
        x0[i] = 0;
        for (j = 0; j < restigual; j++)
        {
            x0[i] = x0[i] + Et[i][j] * y[j];
        }
    }

    i = 0;
    while (i < n)
    { // Enquanto número de indivíduos total da população
        // gera d randomicamente
        for (j = 0; j < dim; j++)
        {
            d[j] = getLower(2, FUNCAO, j) + (rand() % 1000000 / 1000000.f) * (getUpper(2, FUNCAO, j) - getLower(2, FUNCAO, j)); //( (rand()%1000000) / 1000000.f ) * ( bounds[j][1]-bounds[j][0] ) + bounds[j][0]; /// arrumar quem é bounds?
        }

        // calcula o vetor z
        for (j = 0; j < restigual; j++)
        {
            z[j] = 0.0;
            for (k = 0; k < dim; k++)
            {
                z[j] = z[j] + E[j][k] * d[k];
            }
        }

        // Resolve o sistema Mu=z
        u = resolveSistema(M, p, z, restigual);

        // Calcula o vetor v = Et*u
        for (j = 0; j < dim; j++)
        {
            v[j] = 0.0;
            for (k = 0; k < restigual; k++)
            {
                v[j] = v[j] + Et[j][k] * u[k];
            }
        }
        // Calcula a solução xi
        for (j = 0; j < dim; j++)
        {
            x[j] = x0[j] + d[j] - v[j];
            // verifica se ultrapassou o limite superior ou inferior das variáveis
            //   if (x[j] < getLower(2, FUNCAO, j)) //LOWER
            //       	x[j] = getLower(2, FUNCAO, j);
            //       else if (x[j] > getUpper(2, FUNCAO, j)) //UPPER
            //       	x[j] = getUpper(2, FUNCAO, j);
        }

        for (j = 0; j < dim; j++)
        {
            pop[i][j] = x[j];
        }
        calculaAptidao(pop[i], dim, 2, leader, pop[i]);
        // Calcula_Residuo(n, pop[i], pop, restigual, dim, E, c, Residuo, M, p, y, Et, vetor_correto);

        i = i + 1;
    }

    //    for (int i = 0 ; i < n ; i++){
    //        cout << "u[" << i << "]: " ;
    //        for (int j = 0 ; j < dim ; j++){
    //            cout << " " << pop[i][j];
    //        }
    //        cout << "\t " << pop[i][dim] << "\t " << pop[i][dim+1]  << "\t " << pop[i][dim+2];
    //        cout << endl;
    //    }
    //    getchar();

    // Liberando a Memória
    for (i = 0; i < n; ++i)
    {
        free(Residuo[i]);
    }
    for (i = 0; i < dim; i++)
    {
        free(Et[i]);
    }

    for (i = 0; i < restigual; i++)
    {
        free(E[i]);
        free(M[i]);
    }

    free(p);
    free(y);
    free(Residuo);
    free(E);
    free(Et);
    free(M);
    free(c);
    free(vetor_correto);
}

void inicializa(double **&pop, int n, int d, int nivel)
{

    pop = new double *[n];

    for (int i = 0; i < n; i++)
    {

        pop[i] = new double[d + 3];
        // da posição 0 a d-1 => valores das variáveis
        // posição d => aptidao
        // posição d+1 => valor de violação das restrições

        for (int j = 0; j < d; j++)
        {
            pop[i][j] = getLower(nivel, FUNCAO, j) + (rand() / (double)RAND_MAX) * (getUpper(nivel, FUNCAO, j) - getLower(nivel, FUNCAO, j)); // UPPER - LOWER
        }

        if (nivel == 1)
        { // se lider, determina valores do seguidor

            deFollower(pop[i], popLValoresF[i]);

            calculaAptidao(pop[i], DIML, 1, pop[i], popLValoresF[i]);

            // verifica se seguidor viola restrição
            if (popLValoresF[i][DIMF + 1] > 0 || popLValoresF[i][DIMF + 2] > 0)
            {
                // acumula valor da restricao do seguidor
                // nas restricoes do lider
                pop[i][DIML + 1] = pop[i][DIML + 1] + PENALTY * popLValoresF[i][DIMF + 1];
                pop[i][DIML + 1] = pop[i][DIML + 1] + PENALTY * popLValoresF[i][DIMF + 2];
            }
        }
        else
        {
            calculaAptidao(pop[i], d, 2, NULL, NULL);
        }
    }
}

void BlDE()
{
    imprimeCabecalho();

    for (int g = 0; g < GENL; g++)
    {

        deLeader();

        double *uL = new double[DIML + 3];
        int m = selecionaMelhor(uL, popL, SIZEL, DIML, 2);
        cout << "G-" << g << " [Leader] ";
        for (int j = 0; j < DIML; j++)
        {
            cout << uL[j] << " ";
        }
        cout << "Fit: " << uL[DIML] << " Const: " << uL[DIML + 1] << " [Follower] ";
        for (int j = 0; j < DIMF; j++)
        {
            cout << popLValoresF[m][j] << " ";
        }
        cout << "Fit: " << popLValoresF[m][DIMF] << " Const: " << popLValoresF[m][DIMF + 1] << " " << popLValoresF[m][DIMF + 2] << " " << getNEval(1) << " " << getNEval(2) << endl;
        delete[] uL;
    }
}

int main(int argc, char *argv[])
{

    /* Parametros de entrada
    int SIZEL, SIZEF;
    int GENL, GENF;
    int SEED;
    int FUNCAO;
    double F, CR;
    int VARIANTE;*/

    for (int i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-genL") == 0)
        {
            GENL = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-popL") == 0)
        {
            SIZEL = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-genF") == 0)
        {
            GENF = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-popF") == 0)
        {
            SIZEF = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-seed") == 0)
        {
            SEED = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-func") == 0)
        {
            FUNCAO = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-F") == 0)
        {
            F = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-CR") == 0)
        {
            CR = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-tol") == 0)
        {
            TOLERANCE = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-eps") == 0)
        {
            EPS = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-Var") == 0)
        {
            VARIANTE = atoi(argv[++i]);
        }
    }

    srand(SEED);

    SIZEF = SIZEL;

    DIML = getDimensao(FUNCAO, 1);
    DIMF = getDimensao(FUNCAO, 2);

    inicializa(popLValoresF, SIZEL, DIMF, 2);
    inicializa(popL, SIZEL, DIML, 1);
    inicializa(popLNova, SIZEL, DIML, 2);

    BlDE();
}
