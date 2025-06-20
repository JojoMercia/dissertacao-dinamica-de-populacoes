#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

const int NumbPopulations = 2;

#define NumbNeurons 500

const int NumbNeuronsExc = (int)(0.8 * NumbNeurons); //escrever com memória dinâmica pop[i].
const int NumbNeuronsInh = NumbNeurons - NumbNeuronsExc;

#define NumbSynInPop 50 //0.1 * NumbNeurons //10% of connectivity
#define NumbSynBetweenPop 20

const double dt = 0.05;

#include "../auxiliares/nrecipes.c"
#include "populacaobib.c"

#define ALLOC(size, type) (type *)malloc((size_t)size * sizeof(type));

int main(int argc, char const *argv[]){
    clock_t start_time, end_time;
    double elapsed_time;
    start_time = clock();
    char namefile[100];
    FILE *file[NumbPopulations], *filem;
    if(argc < 8){
        fputs("Recommend parameters: 53408123 5000.0 -5.0 0.5 4.0 5.0 0.5 0.0\n", stderr);
        return 1;
    }
    long int seed = -atol(argv[1]);
    //----------------Open files------------------//
    for(int i = 0; i < NumbPopulations; i++){
        sprintf(namefile, "Resultados/excitatorio_TxX_seed_53408123/Vxt[Pop=%d][X=%d][GABAInR=%d].dat", i + 1, (int)atof(argv[3]), (int)atof(argv[6]));
        file[i] = fopen(namefile, "w");
        if (file[i] == NULL) {
            perror("Error opening file\n");
            return 2; 
        }
    }
    /////////////////////////////INITS/////////////////////////////
    //--------------------Neuronal Parameters-------------------//
    double X = atof(argv[3]), Y = 2.0 * X / 5.0;
    double Xi = atof(argv[8]);
    double sigma[2]; //random number [0,1]
    Populacao pop[NumbPopulations];

    BehavIZ randombhvexc = {
        .a = 0.02,
        .b = 0.2
    };
    BehavIZ randombhvinh = {
        .c = -65.0,
        .d = 2.0
    };

    for(int i = 0; i < NumbNeuronsExc; i++){
        // excitatory pop Sender
        sigma[0] = ran2(&seed);
        sigma[0] *= sigma[0];

        randombhvexc.c = -65.0 + (15.0 * sigma[0]);
        randombhvexc.d = 8.0 - (6.0 * sigma[0]);

        pop[0].neuron[i] = InitIZ(randombhvexc); 
    }    
    for(int i = NumbNeuronsExc; i < NumbNeurons; i++){
        // inhibitory pop Sender
        sigma[0] = ran2(&seed);

        randombhvinh.a = 0.02 + (0.08 * sigma[0]);
        randombhvinh.b = 0.25 - (0.05 * sigma[0]);

        pop[0].neuron[i] = InitIZ(randombhvinh); 
    }
    for(int i = 0; i < NumbNeuronsExc; i++){
        // excitatory pop Receiver
        sigma[0] = ran2(&seed);
        sigma[1] = ran2(&seed);

        sigma[0] *= sigma[0];
        sigma[1] *= sigma[1];

        randombhvexc.c = -55.0 - X + (5.0 + X) * sigma[0] - (10.0 - X) * sigma[1];
        randombhvexc.d = 4.0 + Y - (2.0 + Y) * sigma[0] + (4.0 - Y) * sigma[1];
        
        pop[1].neuron[i] = InitIZ(randombhvexc);
    }
    for(int i = NumbNeuronsExc; i < NumbNeurons; i++){
        // inhibitory pop Receiver
        sigma[0] = ran2(&seed);

        randombhvinh.a = 0.02 + (0.08 * sigma[0]);
        randombhvinh.b = 0.25 - (0.05 * sigma[0]);

        pop[1].neuron[i] = InitIZ(randombhvinh);
    }

    for(int i = 0; i < NumbPopulations; i++){
        // external current
        for(int j = 0; j < NumbNeuronsExc; j++)
            pop[i].neuron[j].Iext = 0.0;
        for(int j = NumbNeuronsExc; j < NumbNeurons; j++)
            pop[i].neuron[j].Iext = 0.0;
    }
    //---------------------Synapses parameters------------------------//
    const SYNType AMPA_IZ_POP = {.E = 0.0, .tau = 5.26, .G = 0.0};
    const SYNType GABAa_IZ_POP = {.E = -65.0, .tau = 5.6, .G = 0.0};
    //----------------Synapses INTO the population----------------//
    SYNType AMPAparamInPop[NumbPopulations], GABAaParamInPop[NumbPopulations];

    for(int i = 0; i < NumbPopulations; i++){
        AMPAparamInPop[i] = AMPA_IZ_POP;
        GABAaParamInPop[i] = GABAa_IZ_POP;
        AMPAparamInPop[i].G = atof(argv[4]); //0.5
    }

    GABAaParamInPop[0].G = atof(argv[5]); //4.0
    GABAaParamInPop[1].G = atof(argv[6]);
    
    for(int i = 0; i < NumbPopulations; i++)
        InitSynPop(&pop[i], &pop[i], AMPAparamInPop[i], GABAaParamInPop[i], seed);

    //---------------Synapses BETWEEN populations----------------//
    SYNType AMPAparamBetPop = AMPA_IZ_POP;
    SYNType GABAaParamBetPop = GABAa_IZ_POP;

    AMPAparamBetPop.G = atof(argv[7]); //0.5

    InitSynPop(&pop[0], &pop[1], AMPAparamBetPop, GABAaParamBetPop, seed); 
    //---------------------Poisson parameters--------------------//
    double const PoissonRate = 2400.0; //Hz
    double tpoisson[NumbPopulations][NumbNeurons];
    double const rPoisson = PoissonRate * (dt / 1000.0); // dt in miliseconds
    int noise[NumbPopulations][NumbNeurons];
    double GPoisson = 0.5;

    for(int i = 0; i < NumbPopulations; i++){
        for(int j = 0; j < NumbNeurons; j++)
            tpoisson[i][j] = 30.0 * 8.0 + 600.0 * ran2(&seed); //?
    }
    //--------------------Filter parameters---------------------//
    int iT = 0, iDT = (int)(0.05/dt);

    //////////////////NUMERICAL INTEGRATION/////////////////////////
    double maxtime = atof(argv[2]);
    double Vavarage;

    for(double t = 0.0; t <= maxtime; t += dt){
        //-------------------------Noise------------------------//
        for(int i = 0; i < NumbPopulations; i++)
            for(int j = 0; j < NumbNeurons; j++)
                // integer representation of if(t > tpoisson[i][j]) and if((1.0 - exp(-rPoisson)) > ran2(&seed))
                noise[i][j] = (int)(t > tpoisson[i][j]) * (int)((1.0 - exp(-rPoisson)) > ran2(&seed));
        //-------------------Synapse evolution------------------//
        for(int j = 0; j < NumbNeurons; j++){
            EvolvePOPEULERSYN(&pop[0].neuron[j], noise[0][j], GPoisson, NumbSynInPop, AMPAparamInPop[0], GABAaParamInPop[0], AMPA_IZ_POP, GABAa_IZ_POP);
            EvolvePOPEULERSYN(&pop[1].neuron[j], noise[1][j], GPoisson, (NumbSynInPop + NumbSynBetweenPop), AMPAparamInPop[1], GABAaParamInPop[1], AMPAparamBetPop, GABAa_IZ_POP);
        }
        //-------------------Neuron evolution-------------------//
        EvolvePOPEULERIZ(pop[0].neuron, NumbNeurons);
        EvolvePOPEULERIZ(pop[1].neuron, NumbNeurons);
        /////////////////////PRINT V IN FILES////////////////////////
        //---------------------Avarage of V--------------------//
        iT++;
        if (iT > iDT){
            iT = 0;
            for (int i = 0; i < NumbPopulations; i++){
                Vavarage = 0.0;
                for (int j = 0; j < NumbNeurons; j++)
                    Vavarage += pop[i].neuron[j].V;
                Vavarage /= NumbNeurons;
                fprintf(file[i], "%lf %lf\n", t, Vavarage);
            }

        }
        /////////////////////END OF PRINTS///////////////////////////
    }

    for(int i = 0; i < NumbPopulations; i++)
        fclose(file[i]);

    end_time = clock();
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Simulation Time: %f minutes\n", elapsed_time/60.0);

    return 0;
}