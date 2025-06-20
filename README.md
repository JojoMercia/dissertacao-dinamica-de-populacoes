# Simulação de Populações Neuronais

Este repositório reúne o código desenvolvido no contexto da dissertação de mestrado em Física (UFAL, 2025), voltado à simulação numérica da dinâmica de populações neuronais acopladas. Utilizando o modelo de Izhikevich, o programa permite investigar como diferentes regimes de heterogeneidade e conectividade afetam a sincronização entre duas populações de neurônios.

A implementação considera sinapses excitatórias (AMPA) e inibitórias (GABA<sub>A</sub>), além de ruído externo modelado por processos estocásticos de Poisson. O código foi projetado de forma modular, com possibilidade de expansão para incluir outros tipos de plasticidade ou topologias mais complexas.

O programa fornece como saída a evolução temporal do potencial de membrana médio de cada população, permitindo posterior análise da atividade coletiva. O tratamento e análise dos dados simulados não estão incluídos neste repositório.

## 📁 Estrutura do Repositório

- [`populacao.c`](./populacao.c): código principal da simulação.
- [`populacaobib.c`](./populacaobib.c): biblioteca com funções auxiliares (inicialização, evolução, cálculo de correntes etc.).
- [`populacaobib.h`](./populacaobib.h): cabeçalho com as definições das estruturas e funções utilizadas.

## 🔐 Licença

Este repositório sob a licença [MIT](./LICENSE).

## ⚙️ Compilação

O código pode ser compilado utilizando o compilador `gcc`. No terminal, execute:

```bash
gcc populacao.c -o populacao.out -lm -O3
```

## 🚀 Execução
Após a compilação, a simulação é executada com:

```bash
./populacao.out seed tmax X gE(IN) gI(S) gI(R) gE(BET) Xi
```
Exemplo:

```bash
./populacao.out 53408123 5000.0 -5.0 0.5 4.0 5.0 0.5 0.0
```
Parâmetros:

- `seed`: semente aleatória
- `tmax`: tempo total de simulação (ms)
- `X`, `Xi`: controlam heterogeneidade excitatória/inibitória em R
- `gE_IN`, `gE_BET`: condutâncias AMPA intra e entre populações
- `gI_S`, `gI_R`: condutâncias GABA<sub>A</sub> nas populações S e R

## 🧠 Modelo

A simulação envolve duas populações de neurônios acopladas. Cada população é composta por:

- 400 neurônios excitatórios
- 100 neurônios inibitórios
- 10% de conectividade sináptica aleatória

As sinapses podem ser:
- Excitatórias: AMPA
- Inibitórias: GABA<sub>A</sub>

Cada neurônio recebe estímulo externo via ruído de Poisson. A dinâmica é resolvida numericamente pelo método de Euler com passo fixo.

A heterogeneidade neuronal é incorporada variando os parâmetros do modelo de Izhikevich a partir de distribuições dependentes dos parâmetros `X` (para neurônios excitatórios) e `Xi` (para neurônios inibitórios) na população **R**. Já os neurônios da população **S** mantêm distribuições fixas ao longo das simulações.

### Variação dos neurônios na população R

#### Neurônios excitatórios:
```c
for(int i = 0; i < NumbNeuronsExc; i++){
    sigma[0] = ran2(&seed);
    sigma[1] = ran2(&seed);

    sigma[0] *= sigma[0];
    sigma[1] *= sigma[1];

    randombhvexc.c = -55.0 - X + (5.0 + X) * sigma[0] - (10.0 - X) * sigma[1];
    randombhvexc.d = 4.0 + Y - (2.0 + Y) * sigma[0] + (4.0 - Y) * sigma[1];

    pop[1].neuron[i] = InitIZ(randombhvexc);
}
```

#### Neurônios inibitórios (distribuição fixa):
```c
for(int i = NumbNeuronsExc; i < NumbNeurons; i++){
    sigma[0] = ran2(&seed);

    randombhvinh.a = 0.02 + (0.08 * sigma[0]);
    randombhvinh.b = 0.25 - (0.05 * sigma[0]);

    pop[1].neuron[i] = InitIZ(randombhvinh);
}
```

#### Neurônios inibitórios (com heterogeneidade via Xi):
```c
for(int i = NumbNeuronsExc; i < NumbNeurons; i++){
    sigma[0] = ran2(&seed);
    sigma[1] = ran2(&seed);

    sigma[0] *= sigma[0];
    sigma[1] *= sigma[1];

    randombhvinh.a = 0.06 - Xi + (0.04 + Xi) * sigma[0] - (0.04 - Xi) * sigma[1];
    randombhvinh.b = -0.625 * randombhvinh.a + 0.262;

    pop[1].neuron[i] = InitIZ(randombhvinh);
}
```

Para simular redes homogêneas, basta substituir `randombhv` por um tipo fixo `IZ_NEURONTYPE`, como `IZ_FS` ou `IZ_LTS`, na função `InitIZ()`.


<!--
## 🔍 Análise
As análises dos dados de saída podem ser feitas com os scripts em Python fornecidos no diretório scripts/.

Exemplo: Detecção de picos de atividade
```python
def get_tspike(filename):
    # Lê o arquivo com dados do potencial médio e detecta os instantes de pico
```
Exemplo: Cálculo de defasagem (tau)
```python
def tau_calculate(tspike, NumbSpike):
    # Calcula o atraso entre as populações baseado nos tempos de disparo
```
-->
## 📎 Dependências
Numerical Recipes — função ran2 para geração de números pseudoaleatórios.

## 🔗 Referência

Se utilizar este código, por favor cite uma das referências abaixo:

- **J. M. G. L. Silva**. *Diversidade na Sincronização de Fase em Redes Neurais*. Dissertação de Mestrado – Programa de Pós-Graduação em Física, Universidade Federal de Alagoas, Maceió, 2025.  
  Disponível em: [http://www.repositorio.ufal.br/jspui/handle/123456789/12670](http://www.repositorio.ufal.br/jspui/handle/123456789/12670)  
  <!-- Atualize o link do repositório após a publicação oficial -->

ou

- **K. V. Brito**, **J. M. Silva**, **C. R. Mirasso** e **F. S. Matias**.  
  *The role of inhibitory neuronal variability in modulating phase diversity between coupled networks*.  
  *arXiv preprint* [arXiv:2411.19664](https://arxiv.org/abs/2411.19664), 2024.

