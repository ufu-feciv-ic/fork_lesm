# Princípios de Montagem e Conectividade no LESM

Este documento descreve como o LESM gerencia a criação de nós e elementos, garantindo a integridade da matriz de rigidez global independentemente da ordem de criação dos componentes.

## 1. Conectividade Baseada em Objetos (Handles)
Diferente de sistemas que dependem exclusivamente de índices numéricos, o LESM utiliza a natureza de **objetos (handle classes)** do MATLAB.
- Quando um elemento (barra) é criado, ele armazena referências diretas aos objetos dos nós inicial e final em sua propriedade `nodes`.
- Isso significa que, internamente, a barra "sabe" exatamente a quais nós está conectada, mesmo que os IDs numéricos dos nós sejam alterados ou criados em ordem não sequencial.

## 2. O Mapa de Graus de Liberdade (Matriz ID)
Antes de iniciar o cálculo da matriz de rigidez, a classe `Drv` executa o método `assembleDOFNum()`. Este método é o "recenseamento" do modelo:
1. Ele percorre todos os nós presentes no vetor `drv.nodes`.
2. Para cada nó, ele verifica quais graus de liberdade (DOF) estão livres ou restritos (EBC - Essential Boundary Conditions).
3. Ele preenche a matriz **`ID(k, n)`**, onde:
   - `k`: Índice do grau de liberdade (ex: 1 para deslocamento em X).
   - `n`: ID do nó.
   - **Valor**: O número da equação global correspondente na matriz de rigidez.

## 3. O Vetor de Coleta (`gle` - Gather Vector)
Cada elemento possui um vetor chamado `gle` que serve como uma "ponte" entre o local e o global. O método `assembleGle()` da classe `Drv` preenche este vetor para cada elemento:
- O elemento consulta a matriz `ID` usando os IDs de seus nós inicial e final.
- Ele armazena esses números de equações globais em sequência no `gle`.
- *Exemplo:* Se uma barra 2D conecta o Nó A (equações 1, 2, 3) ao Nó B (equações 7, 8, 9), seu `gle` será `[1, 2, 3, 7, 8, 9]`.

## 4. Processo de Montagem (Assembly)
A montagem da matriz de rigidez global `K` é feita através de um processo de "espalhamento" (scattering) no método `assembleElemMtx()`:

```matlab
% keg: matriz de rigidez do elemento no sistema global
% gle: vetor de coleta do elemento
drv.K(gle, gle) = drv.K(gle, gle) + keg;
```

O MATLAB utiliza o vetor `gle` para identificar simultaneamente as linhas e colunas em `K` onde os valores da matriz de rigidez do elemento (`keg`) devem ser somados.

## Resumo do Fluxo de Trabalho
1. **Modelagem**: Usuário cria nós e elementos (ordem irrelevante).
2. **Numeração**: `assembleDOFNum()` organiza as equações globais com base na lista de nós.
3. **Mapeamento**: `assembleGle()` vincula as equações de cada nó aos elementos que neles se encontram.
4. **Integração**: `gblMtx()` percorre os elementos e soma suas contribuições em `K` usando os endereços mapeados no `gle`.

Este desacoplamento entre a **criação** (interface) e a **numeração** (processamento) garante que o sistema seja robusto a qualquer sequência de entrada do usuário.

---

## 5. Partição de Matrizes vs. Método da Penalização

Atualmente, o LESM utiliza a **Partição de Matrizes** para resolver o sistema de equações. Abaixo, detalhamos a diferença entre este método e uma possível implementação futura via **Método da Penalização**.

### O Cenário Atual: Partição de Matrizes
Para possibilitar a solução por partição, a numeração dos graus de liberdade no LESM não é apenas sequencial, mas **estratégica**:
1.  **Numeração Prioritária**: O método `assembleDOFNum` numera primeiro todos os graus de liberdade **livres** (`countF`) e somente depois os **presos** (`countS`).
2.  **Objetivo**: Isso permite que a matriz global $K$ seja visualizada em blocos:
    $$ \begin{bmatrix} K_{ff} & K_{fs} \\ K_{sf} & K_{ss} \end{bmatrix} \begin{Bmatrix} D_f \\ D_s \end{Bmatrix} = \begin{Bmatrix} F_f \\ F_s \end{Bmatrix} $$
3.  **Vantagem**: O sistema resolvido é reduzido apenas à parte livre ($K_{ff} \cdot D_f = F_f$), o que é computacionalmente mais eficiente e numericamente mais estável.

### O Cenário Alternativo: Método da Penalização
Se o LESM migrasse para o método de penalização, os seguintes aspectos mudariam:

1.  **Ordem de Numeração Irrelevante**: Não haveria mais necessidade de separar "livres" de "presos" durante a numeração. Os graus de liberdade poderiam seguir a ordem natural dos nós, simplificando o método `assembleDOFNum`.
2.  **Tamanho do Sistema**: O sistema de equações não seria reduzido. Se o modelo tem 1000 DOFs, a matriz resolvida será sempre 1000x1000, independentemente de quantos apoios existam.
3.  **Aplicação de Vínculos**: Em vez de "cortar" a matriz, adiciona-se uma rigidez fictícia extremamente elevada ($C \approx 10^{15}$) nos termos da diagonal de $K$ correspondentes aos graus presos.
    - $K_{ii, novo} = K_{ii} + C$
    - $F_{i, novo} = F_i + C \cdot \delta_{prescrito}$
4.  **Reações de Apoio**: As reações seriam obtidas diretamente pelo erro residual da penalização: $R_i = C \cdot (D_{calculado} - \delta_{prescrito})$.

### Comparativo Resumido

| Característica | Partição (Atual) | Penalização (Opcional) |
| :--- | :--- | :--- |
| **Lógica de Numeração** | Rígida (Livres antes de Presos) | Flexível (Ordem arbitrária) |
| **Tamanho da Matriz Solvida** | Reduzida (apenas DOFs livres) | Total (todos os DOFs) |
| **Estabilidade Numérica** | Alta (Exata) | Depende da escolha da constante $C$ |
| **Cálculo de Reações** | Via álgebra de blocos ($K_{sf} \cdot D_f$) | Proporcional à penalização |

---

## 6. O Impacto do Uso de Matrizes Esparsas

A decisão do LESM de usar Partição de Matrizes reflete uma abordagem tradicional para matrizes densas. No entanto, ao utilizar **Matrizes Esparsas**, o Método da Penalização torna-se muito mais vantajoso:

1.  **Integridade da Estrutura de Dados**: Matrizes esparsas dependem do "padrão de esparsidade" (sparsity pattern). Particionar a matriz exige a criação de novas estruturas esparsas, o que envolve realocação de memória e reindexação cara. A penalização permite manter a mesma estrutura do início ao fim.
2.  **Custo Computacional**: Em sistemas esparsos, o benefício de reduzir o tamanho da matriz (ex: de 1000x1000 para 900x900) é frequentemente superado pelo custo de processamento necessário para realizar a partição. Solvers esparsos modernos lidam com o sistema completo de forma extremamente eficiente.
3.  **Simplicidade de Código**: O uso de penalização com matrizes esparsas elimina a necessidade de numerar graus de liberdade livres e presos separadamente, simplificando drasticamente os algoritmos de montagem (`Assembly`).

**Conclusão para Novos Desenvolvimentos**: Se o foco do projeto for o uso de matrizes esparsas, o Método da Penalização é a escolha técnica mais coerente, unindo simplicidade de implementação com alta performance.

---

## 7. Mapeamento de Índices e Estratégia de Tripletos (I, J, V)

Mesmo utilizando matrizes esparsas e o método da penalização, ainda é necessário um "tradutor" para gerenciar a relação entre o que o usuário vê e como o computador armazena os dados.

### 7.1. IDs de Usuário vs. Índices de Memória
O usuário pode criar nós com IDs arbitrários (ex: Nó 1 e Nó 50). 
- **O Problema**: Criar uma matriz baseada diretamente no ID "50" geraria um enorme vazio de memória entre os índices 2 e 49.
- **A Solução**: O sistema mantém um vetor sequencial de objetos. O "Nó 50" pode ser o segundo objeto no vetor (`nodes[1]`). O tradutor converte o ID do usuário para o **índice sequencial** da lista de objetos.

### 7.2. Mapeamento de Graus de Liberdade (DOF)
O tradutor também converte o índice do nó nos índices das equações globais:
- Se cada nó tem 3 DOFs, o nó no índice sequencial `n` terá suas equações em:
  - $Eq_1 = (n-1) \times 3 + 1$
  - $Eq_2 = (n-1) \times 3 + 2$
  - $Eq_3 = (n-1) \times 3 + 3$

### 7.3. Montagem Eficiente: Estratégia de Tripletos
Em matrizes esparsas, inserir valores um a um diretamente na matriz final é ineficiente. A técnica recomendada é o uso de **Tripletos**:
1.  **Coleta**: Durante o loop de elementos, em vez de tocar na matriz global, os dados são armazenados em três vetores auxiliares:
    - `I`: Vetor com os índices das linhas.
    - `J`: Vetor com os índices das colunas.
    - `V`: Vetor com os valores numéricos.
2.  **Expansão**: Para cada barra, adicionamos todas as contribuições de sua matriz de rigidez local ($k_{local}$) nesses vetores.
3.  **Compressão Final**: Após percorrer todos os elementos, uma única função (como `sparse(I, J, V)` no MATLAB ou `setFromTriplets` em C++) constrói a matriz esparsa final. Se houver índices repetidos em `I` e `J`, o sistema automaticamente soma os valores correspondentes em `V`, realizando o **Assembly** de forma extremamente veloz.

**Conclusão**: A escolha atual do LESM pela partição justifica a complexidade extra no `assembleDOFNum` em troca de precisão numérica e menor custo computacional na resolução do sistema linear.
