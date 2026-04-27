# Estratégias de Montagem e Solução: Particionamento vs. Matrizes Esparsas

Este documento compara a abordagem atual do LESM (Particionamento de Matrizes) com a alternativa de alto desempenho (Matrizes Esparsas e Numeração Sequencial).

---

## 1. Abordagem Atual: Particionamento Manual
O LESM organiza a numeração das equações de forma que todos os graus de liberdade (DOFs) livres venham primeiro, seguidos pelos DOFs fixos (apoios).

### Como funciona:
- **Numeração (ID):** O vetor `ID` é gerado para que índices $1$ a $neqfree$ sejam livres.
- **Montagem:** O `gle` mapeia os termos locais para esses índices "separados".
- **Solução:** A matriz global $K$ é fatiada em submatrizes ($K_{ff}, K_{fs}, K_{sf}, K_{ss}$).

### Vantagens:
- **Didática:** É a transcrição direta da formulação teórica do Método da Rigidez Direta.
- **Extração de Reações:** As reações de apoio ($F_s$) são calculadas explicitamente através da operação $K_{sf} \cdot D_f$.

---

## 2. Abordagem Alternativa: Matrizes Esparsas e Numeração Sequencial
Nesta abordagem, a numeração segue a ordem física dos nós, e a matriz é tratada como um objeto esparso (armazenando apenas valores não-zero).

### Como funciona:
- **Numeração:** Não há distinção entre livre e fixo na contagem. O nó 1 tem DOFs 1, 2, 3; o nó 2 tem 4, 5, 6, e assim por diante.
- **Montagem (Triplets):** Em vez de somar na matriz global dentro de um loop (lento), coletamos as coordenadas $(I, J)$ e os valores $V$ de cada elemento e criamos a matriz de uma vez: `K = sparse(I, J, V)`.
- **Condições de Contorno:** Como a matriz não é particionada, os apoios são tratados "zerando" a linha/coluna correspondente e colocando `1` na diagonal (ou via método da penalidade).

### Vantagens:
- **Performance:** Reduz drasticamente a "Largura de Banda" (Bandwidth), o que torna a solução muito mais rápida para modelos grandes.
- **Memória:** Economiza RAM ao não armazenar os zeros entre os termos de rigidez.

---

## 3. O Papel do `gle` em Ambos os Cenários

Uma dúvida comum é se o **`gle` (Gather Vector)** ainda seria necessário com matrizes esparsas. **A resposta é SIM.**

O `gle` não existe apenas para particionar; sua função principal é o **Mapeamento Topológico**.

| Função do `gle` | No Particionamento (Atual) | Na Matriz Esparsa (Global) |
| :--- | :--- | :--- |
| **Mapeamento** | Diz onde o termo local vai na partição livre ou fixa. | Diz onde o termo local vai na matriz global sequencial. |
| **Uso na Montagem** | `K(gle, gle) = K(gle, gle) + ke` | Fornece os índices para os vetores de coordenadas `I` e `J`. |
| **Dependência** | Depende da classificação do nó (livre/fixo). | Depende apenas da conectividade entre os nós. |

---

## 4. Comparação Técnica

| Característica | Particionamento (LESM) | Matrizes Esparsas (Indústria) |
| :--- | :--- | :--- |
| **Estrutura de $K$** | Blocos lógicos ($K_{ff}, K_{fs}$) | Matriz única com padrão de "faixa" |
| **Armazenamento** | Geralmente densa (Matlab `zeros`) | Formato comprimido (CCS/CRS) |
| **Solver** | Inversão/Fatoração de $K_{ff}$ | Solvers diretos (Backslash `\`, SparseLU) |
| **Reações de Apoio** | Calculadas por multiplicação de blocos | Extraídas do resíduo $F - K \cdot D$ |
| **Escalabilidade** | Limitada por memória/largura de banda | Alta (milhares/milhões de DOFs) |

## 5. Conclusão
O **particionamento** é excelente para entender a mecânica do método e para pequenos problemas de engenharia. No entanto, para evoluir o LESM para um solver de alta performance (similar ao que bibliotecas como a **Eigen** permitem em C++), a migração para **matrizes esparsas com numeração sequencial** é o caminho natural, mantendo o `gle` como a peça fundamental de integração.
