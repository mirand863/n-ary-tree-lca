#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

class LCA {
    private:
        unordered_map<string, vector<string>> parents;
        vector<string> euler;
        vector<int> depth;
        unordered_map<string, int> firstAppearance;
//        int **M;
        void depthFirstSearch(string current, int depth);
//        void preProcessRMQ();
//        int queryRMQ(int i, int j);
    public:
        LCA();
        virtual ~LCA();
        void addEdge(string father, string son);
        void doEulerWalk();
        void printEuler();
        void printLevel();
        void printFirstAppearance();
//        string getLCA(string u, string v);
};

void LCA::printEuler() {
    cout << "Euler:" << endl;
    for(unsigned int i = 0; i < euler.size(); i++) {
        cout << euler[i] << " ";
    }
    cout << endl;
}

void LCA::printLevel() {
    cout << "Level:" << endl;
    for(unsigned int i = 0; i < depth.size(); i++) {
        cout << depth[i] << " ";
    }
    cout << endl;
}

void LCA::printFirstAppearance() {
    cout << "First appearance:" << endl;
    for(unordered_map<string, int>::iterator it = firstAppearance.begin(); it != firstAppearance.end(); ++it) {
        cout << it->first << "->" << it->second << " ";
    }
    cout << endl;
}

LCA::LCA() {
}

LCA::~LCA() {
}

void LCA::addEdge(string father, string son) {
    if(parents.count(father) == 0) {
        vector<string> children;
        children.push_back(son);
        parents[father] = children;
        firstAppearance[father] = -1;
    } else {
        parents[father].push_back(son);
    }
    firstAppearance[son] = -1;
}

void LCA::depthFirstSearch(string current, int depth) {
    // marking first appearance for current node
    if(firstAppearance[current] == -1) {
        firstAppearance[current] = euler.size();
    }
    // pushing root to euler walk
    euler.push_back(current);
    // pushing depth of current node
    this->depth.push_back(depth);
    for(unsigned int i = 0; i < parents[current].size(); i++) {
        depthFirstSearch(parents[current][i], depth + 1);
        euler.push_back(current);
        this->depth.push_back(depth);
    }
}

void LCA::doEulerWalk() {
    depthFirstSearch("1", 0);
//    preProcessRMQ();
}

//// <O(N logN) Preprocessing time, O(1) Query time>
//void LCA::preProcessRMQ() {
//    M = new int*[depth.size()];
//    int logDepth = log2(depth.size());
//    for(unsigned int i = 0; i < depth.size(); i++) {
//        M[i] = new int[logDepth];
//    }
//
//    //initialize M for the intervals with length 1
//    for(unsigned int i = 0; i < depth.size(); i++) {
//        M[i][0] = i;
//    }
//    //compute values from smaller to bigger intervals
//    for(unsigned int j = 1; 1 << j <= depth.size(); j++) {
//        for(unsigned int i = 0; i + (1 << j) - 1 < depth.size(); i++) {
//            if(depth[M[i][j - 1]] < depth[M[i + (1 << (j - 1))][j - 1]]) {
//                M[i][j] = M[i][j - 1];
//            } else {
//                M[i][j] = M[i + (1 << (j - 1))][j - 1];
//            }
//        }
//    }
//}
//
//int LCA::queryRMQ(int i, int j) {
//    if(i > j) {
//        swap(i, j);
//    }
//
//    int k = log2(j - i + 1);
//
//    if(depth[M[i][k]] <= depth[M[j - (1 << k) + 1][k]]) {
//        return M[i][k];
//    } else {
//        return M[j - (1 << k) + 1][k];
//    }
//}

//string LCA::getLCA(string u, string v) {
////    u = encode[u];
////    v = encode[v];
//	// trivial case
//	if (u.compare(v) == 0) {
//        return u;
//	}
//
//	if(firstAppearance[u] > firstAppearance[v]) {
//        swap(u, v);
//	}
//
//	// doing RMQ in the required range
//	return euler[queryRMQ(firstAppearance[u], firstAppearance[v])];
//}

//string getLCA(LCA &tree, vector<string> &taxIds) {
//    string lca;
//    if(taxIds.size() >= 2) {
//        lca = tree.getLCA(taxIds[0], taxIds[1]);
//        for(int i = 2; i < taxIds.size(); i++) {
//            lca = tree.getLCA(lca, taxIds[i]);
//        }
//        return lca;
//    } else {
//        return taxIds[0];
//    }
//}

int main(int argc, char *argv[]) {
    if(argc == 3) {
        string father, son;
        LCA lca;
        ifstream tree(argv[1]);
        while(tree >> father >> son) {
            lca.addEdge(father, son);
        }
        lca.doEulerWalk();
        lca.printEuler();
        lca.printLevel();
        lca.printFirstAppearance();

//        string currentRead, nextRead;
//        string currentTaxId, nextTaxId;
//        int currentKmer, nextKmer;
//        vector<string> taxIds;
//        ifstream queries(argv[2]);
//        queries >> currentRead >> currentTaxId >> currentKmer;
//        taxIds.push_back(currentTaxId);
//        while(queries >> nextRead >> nextTaxId >> nextKmer) {
//            if(nextRead.compare(currentRead) == 0) {
//                taxIds.push_back(nextTaxId);
//            } else {
//                cout << currentRead << "\t";
//                cout << getLCA(lca, taxIds) << "\n";
//                taxIds.clear();
//                currentRead = nextRead;
//                taxIds.push_back(nextTaxId);
//            }
//        }
//        cout << currentRead << "\t";
//        cout << getLCA(lca, taxIds) << "\n";
    } else {
        cout << "Usage: " << argv[0] << " tree.tsv queries.tsv\n";
    }
}
