#include <SU3_ops.h>
#include <geometry.h>

void actDiracOperator(Scalar *f, Scalar *g) {
    //	Calcula g=D.f, ou seja, a matriz de Dirac agindo no vetor f, e coloca o resultado em g.

    PosVec position;
    PosVec positionl1;
    PosVec positionl2;

    int posicaol1[d];
    int posicaol2[d];

    int versormu[d];
    int versormenosmu[d];

    double complex Elo1[3][3];
    double complex produto1;

    double complex Elo2[3][3];
    double complex produto2;

    for (int mu = 0; mu < d; mu++) {
        versormu[mu] = 0;
        versormenosmu[mu] = 0;
    }

    for (posicao[0] = 0; posicao[0] < Nt; posicao[0]++)
        for (posicao[1] = 0; posicao[1] < Nxyz; posicao[1]++)
            for (posicao[2] = 0; posicao[2] < Nxyz; posicao[2]++)
                for (posicao[3] = 0; posicao[3] < Nxyz; posicao[3]++)
                    for (int alfa = 0; alfa < 4; alfa++)
                        for (int a = 0; a < 3; a++) {
                            // Termo diagonal, simplesmente copia o conteúdo de f para g

                            (*(ELEM_VEC_INV(g, posicao, alfa, a))) = (*(ELEM_VEC_INV(f, posicao, alfa, a)));

                            // Termo de hopping com vizinhos da "frente" (mu positivo, 1) e de "tras" (mu negativo, 2)

                            for (int mu = 0; mu < d; mu++) {
                                versormu[mu] = 1;
                                versormenosmu[mu] = -1;

                                EloVizinhoSU3(posicao, mu, FRENTE, Elo1);
                                EloVizinhoSU3(posicao, mu, TRAS, Elo2);

                                SomaVetoresPosicao(posicao, versormu, posicaol1);
                                SomaVetoresPosicao(posicao, versormenosmu, posicaol2);

                                for (int alfal = 0; alfal < 4; alfal++)
                                    for (int al = 0; al < 3; al++) {
                                        produto1 = menoskappa * IdentidadeDiracMenosGama[mu][alfa][alfal] * Elo1[a][al] * (*(ELEM_VEC_INV(f, posicaol1, alfal, al)));
                                        produto2 = menoskappa * IdentidadeDiracMaisGama[mu][alfa][alfal] * Elo2[a][al] * (*(ELEM_VEC_INV(f, posicaol2, alfal, al)));

                                        //	Condição de contorno antiperiódica na direção temporal
                                        (*(ELEM_VEC_INV(g, posicao, alfa, a))) += (posicao[0] != (Nt - 1) || mu != 0) ? produto1 : -produto1;
                                        (*(ELEM_VEC_INV(g, posicao, alfa, a))) += (posicao[0] != 0 || mu != 0) ? produto2 : -produto2;
                                    }

                                versormu[mu] = 0;
                                versormenosmu[mu] = 0;
                            }

                            // Melhoria de trevo /* Vou ignorar improvements por enquanto */

                            // for(int alfal=0;alfal<4;alfal++)
                            // 	for(int al=0;al<3;al++){
                            // 		(*(ELEM_VEC_INV(g, posicao, alfa, a))) += cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al] * (*(ELEM_VEC_INV(f, posicao, alfal, al)));
                            // }
                        }
}