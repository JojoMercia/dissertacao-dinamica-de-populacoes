#if !defined(POP_HEADER)
#define POP_HEADER

typedef struct BehavIZ{
    double a, b, c, d;
} BehavIZ;

const BehavIZ IZ_RS = {.a = 0.02, .b = 0.2, .c = -65.0, .d = 8.0};
const BehavIZ IZ_IB = {.a = 0.02, .b = 0.2, .c = -55.0, .d = 4.0};
const BehavIZ IZ_CH = {.a = 0.02, .b = 0.2, .c = -50.0, .d = 2.0};

//INHIBITORY
const BehavIZ IZ_FS = {.a = 0.1, .b = 0.2, .c = -65.0, .d = 2.0};
const BehavIZ IZ_LTS = {.a = 0.02, .b = 0.25, .c = -65.0, .d = 2.0};

typedef struct SYNType{
    double G, E, tau;
} SYNType;

typedef struct Receptors{
    double Exc, Inh, Poisson;
} Receptors;

Receptors rcptrs = {
    .Exc = 0.0,
    .Inh = 0.0,
    .Poisson = 0.0};

typedef struct NeuronIZ{
    double V, u;
    double a, b, c, d;
    double Isyn, Iext;
    bool spiketime;
    struct Receptors r;
    struct NeuronIZ *pre[NumbSynBetweenPop + NumbSynInPop];
    struct SYNType st;
} NeuronIZ;

typedef struct Populacao{
    int NumbSynAMPA, NumbSynGABAa;
    NeuronIZ neuron[NumbNeurons];
} Populacao;

extern inline double dotVIZ(double, double, double);
extern inline double dotu(double, double, double, double);
extern inline double popdotr(double, double, int);
extern inline double Isyn(double, double, double, double);
extern inline void InitSynPop(Populacao*, Populacao*, SYNType, SYNType, long);
extern inline void EvolvePOPEULERIZ(NeuronIZ*, int);
extern inline void EvolvePOPEULERSYN(NeuronIZ*, int, double, int, SYNType, SYNType, SYNType, SYNType);
#endif