#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

/**
 * Exemplo de Protótipo C++ para o LESM
 * Modelo com 20 nós para demonstrar Eficiência de Matriz Esparsa
 */

struct Node {
    int id;
    double x, y;
    bool fixedX, fixedY; 
};

struct Element {
    int id;
    int nodeStart, nodeEnd;
    double E; 
    double A; 
};

int main() {
    // 1. GERAÇÃO AUTOMÁTICA DE 20 NÓS (Em linha no eixo X)
    int numNodes = 20;
    vector<Node> nodes;
    for (int i = 0; i < numNodes; i++) {
        // Primeiro nó (0) é fixo em X e Y
        // Último nó (19) é fixo apenas em Y
        bool fx = (i == 0);
        bool fy = (i == 0 || i == numNodes - 1);
        nodes.push_back({i, (double)i, 0.0, fx, fy});
    }

    // 2. GERAÇÃO AUTOMÁTICA DE 19 ELEMENTOS (Conectando i -> i+1)
    vector<Element> elements;
    for (int i = 0; i < numNodes - 1; i++) {
        elements.push_back({i, i, i + 1, 210e6, 0.01});
    }

    int ndofPerNode = 2; 
    int totalDofs = numNodes * ndofPerNode;

    // 3. MONTAGEM DA MATRIZ ID (Numeração Global)
    MatrixXi idMatrix(ndofPerNode, numNodes);
    int countFree = 0;
    for (auto& n : nodes) {
        if (!n.fixedX) countFree++;
        if (!n.fixedY) countFree++;
    }
    
    int nextFree = 1;
    int nextFixed = countFree + 1;
    for (int n = 0; n < numNodes; n++) {
        idMatrix(0, n) = (!nodes[n].fixedX) ? nextFree++ : nextFixed++;
        idMatrix(1, n) = (!nodes[n].fixedY) ? nextFree++ : nextFixed++;
    }

    // 4. PREPARAÇÃO DAS MATRIZES
    MatrixXd K_dense = MatrixXd::Zero(totalDofs, totalDofs);
    vector<Triplet<double>> tripletList;
    tripletList.reserve(elements.size() * 16);

    // 5. LOOP NOS ELEMENTOS E ASSEMBLY
    for (auto& e : elements) {
        Node& n1 = nodes[e.nodeStart];
        Node& n2 = nodes[e.nodeEnd];

        double L = sqrt(pow(n2.x - n1.x, 2) + pow(n2.y - n1.y, 2));
        double c = (n2.x - n1.x) / L;
        double s = (n2.y - n1.y) / L;

        // Mapeamento GLE
        vector<int> gle = { idMatrix(0, e.nodeStart), idMatrix(1, e.nodeStart), 
                            idMatrix(0, e.nodeEnd),   idMatrix(1, e.nodeEnd) };

        double k_coeff = (e.E * e.A) / L;
        Matrix4d ke;
        ke <<  c*c,  c*s, -c*c, -c*s,
               c*s,  s*s, -c*s, -s*s,
              -c*c, -c*s,  c*c,  c*s,
              -c*s, -s*s,  c*s,  s*s;
        ke *= k_coeff;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int row = gle[i] - 1;
                int col = gle[j] - 1;
                double val = ke(i, j);

                K_dense(row, col) += val;
                tripletList.push_back(Triplet<double>(row, col, val));
            }
        }
    }

    // Matriz Esparsa Final
    SparseMatrix<double> K_sparse(totalDofs, totalDofs);
    K_sparse.setFromTriplets(tripletList.begin(), tripletList.end());

    // 6. EXIBIÇÃO E PLOT DA TOPOLOGIA
    cout << fixed << setprecision(2);
    cout << "=========================================================" << endl;
    cout << "   ESTRUTURA COM 20 NOS E 19 ELEMENTOS (40 DOFs)" << endl;
    cout << "=========================================================\n" << endl;

    cout << "Estatisticas da Matriz Global:" << endl;
    cout << "-> Dimensao: " << totalDofs << "x" << totalDofs << endl;
    cout << "-> Espacos na Matriz Densa: " << K_dense.size() << " (100% de ocupacao)" << endl;
    cout << "-> Entradas Nao-Nulas na Esparsa: " << K_sparse.nonZeros() 
         << " (" << (double)K_sparse.nonZeros() / K_dense.size() * 100.0 << "% de ocupacao)" << endl;

    cout << "\n--- PLOT DA TOPOLOGIA (SPARSITY PATTERN) ---" << endl;
    cout << "Legenda: X = Rigidez armazenada, . = Posicao vazia (zero)" << endl;
    
    // Cabeçalho dos índices
    cout << "     ";
    for(int j=0; j<totalDofs; j++) {
        if (j % 5 == 0) cout << setw(2) << j << " ";
        else cout << "   ";
    }
    cout << endl;

    for (int i = 0; i < totalDofs; i++) {
        cout << setw(2) << i << " [ ";
        for (int j = 0; j < totalDofs; j++) {
            if (abs(K_sparse.coeff(i, j)) > 1e-9) cout << "X  ";
            else cout << ".  ";
        }
        cout << "]" << endl;
    }

    cout << "\nANALISE DO PLOT:" << endl;
    cout << "1. Note como os 'X' formam uma banda diagonal." << endl;
    cout << "2. A grande maioria da matriz eh preenchida por '.' (zeros)." << endl;
    cout << "3. Quanto maior a estrutura, mais pontos '.' teremos e mais eficiente sera a Esparsa." << endl;

    return 0;
}
