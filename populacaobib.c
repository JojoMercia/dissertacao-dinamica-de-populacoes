#include "populacaobib.h"

NeuronIZ InitIZ(BehavIZ bhv){
    NeuronIZ n;

    n.a = bhv.a;
    n.b = bhv.b;
    n.c = bhv.c;
    n.d = bhv.d;
    n.V = -60.0;
    n.u = 0.2 * n.V;
    for (int i = 0; i < NumbSynBetweenPop + NumbSynInPop; i++)
        n.pre[i] = NULL;
    n.r = rcptrs;

    return n;
}

inline void InitSynPop(Populacao *pop1, Populacao *pop2, SYNType AMPA, SYNType GABAa, long seed){
    int n = 0, N = NumbSynInPop, NumbMax = NumbNeurons;
    int nrand;

    if(pop1 != pop2){
    //synapse between the populations
        n = NumbSynInPop;
        N += NumbSynBetweenPop;
        NumbMax = NumbNeuronsExc;
    }

    for (int k = 0; k < NumbNeurons; k++){
        for (int j = n; j < N; j++){
            do nrand = (int)floor(NumbMax * ran2(&seed));
            while (k == nrand);

            pop2->neuron[k].pre[j] = &(pop1->neuron[nrand]);

            if (nrand < NumbNeuronsExc) pop2->neuron[k].pre[j]->st = AMPA;
            if (nrand >= NumbNeuronsExc) pop2->neuron[k].pre[j]->st = GABAa;
        }
    }

    return;
}

inline double dotVIZ(double V, double u, double Isyn){
    return 0.04 * V * V + 5.0 * V + 140.0 - u + Isyn;
}
inline double dotu(double V, double u, double a, double b){
    return a * (b * V - u);
}
inline void EvolvePOPEULERIZ(NeuronIZ neuroniz[], int n){
    double kV, ku, Vreset = 30.0;

    for (int i = 0; i < n; i++){
        kV = dotVIZ(neuroniz[i].V, neuroniz[i].u, neuroniz[i].Isyn + neuroniz[i].Iext);
        ku = dotu(neuroniz[i].V, neuroniz[i].u, neuroniz[i].a, neuroniz[i].b);

        neuroniz[i].V += kV * dt;
        neuroniz[i].u += ku * dt;

        //-------------------Reset Izhikevich-------------------//
        if (neuroniz[i].V >= Vreset){
            neuroniz[i].V = neuroniz[i].c;
            neuroniz[i].u += neuroniz[i].d;
            neuroniz[i].spiketime = true;
        }
        else neuroniz[i].spiketime = false;
    }

    return;
}

inline double popdotr(double r, double tau, int n){
    return (-r + n) / tau;
}
inline double Isyn(double V, double r, double G, double E){
    return r * G * (E - V);
}

inline void EvolvePOPEULERSYN(NeuronIZ *n, int noise, double GPoisson, int N, SYNType AMPA, SYNType GABAa, SYNType AMPABet, SYNType GABAaBet){
    int NumbSpikeExc = 0, NumbSpikeInh = 0;
    //For a neuron

    n->Isyn = 0.0;
    
    for (int i = 0; i < N; i++){ //looks for all pre
        if (n->pre[i]->spiketime){ // if its potential spiked
            if(n->pre[i]->st.tau == AMPA.tau || n->pre[i]->st.tau == AMPABet.tau) NumbSpikeExc++;
            if(n->pre[i]->st.tau == GABAa.tau || n->pre[i]->st.tau == GABAaBet.tau) NumbSpikeInh++;
        }
    }   

    //Excitatory synapses
    n->r.Exc += popdotr(n->r.Exc, AMPA.tau, NumbSpikeExc) * dt;
    n->Isyn += Isyn(n->V, n->r.Exc, AMPA.G, AMPA.E);
    //Inhibitory synapses
    n->r.Inh += popdotr(n->r.Inh, GABAa.tau, NumbSpikeInh) * dt;
    n->Isyn += Isyn(n->V, n->r.Inh, GABAa.G, GABAa.E);
    //Poisson synapses
    n->r.Poisson += popdotr(n->r.Poisson, AMPA.tau, noise) * dt;
    n->Isyn += Isyn(n->V, n->r.Poisson, GPoisson, AMPA.E);

    return;
}