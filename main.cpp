#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <cmath>
#include <set>

using namespace std;

class LCA {
    private:
        bool **tree;
        int *parent;
        int vertices;
        vector<int> euler;
        vector<int> depth;
        int *firstAppearance;
        unordered_map<int, int> encode, decode;
        int **M;
        int getRoot();
        void depthFirstSearch(int current, int depth);
        void preProcessRMQ();
        int queryRMQ(int i, int j);
        void doEncode(int node);
    public:
        LCA(int vertices);
        virtual ~LCA();
        void addEdge(int father, int son);
        void doEulerWalk();
        int getLCA(int u, int v);
};

LCA::LCA(int vertices) {
    this->vertices = vertices;
    tree = new bool*[vertices];
    parent = new int[vertices];
    firstAppearance = new int[vertices];
    for(int i = 0; i < vertices; i++) {
        tree[i] = new bool[vertices]();
        parent[i] = -1;
        firstAppearance[i] = -1;
    }
}

LCA::~LCA() {
    for(int i = 0; i < vertices; i++) {
        delete [] tree[i];
    }
    delete [] tree;
    delete [] parent;
    delete [] firstAppearance;
    for(unsigned int i = 0; i < depth.size(); i++) {
        delete [] M[i];
    }
    delete [] M;
}

void LCA::doEncode(int node) {
    if(encode.count(node) == 0) {
        encode.insert({node, encode.size()});
        decode.insert({decode.size(), node});
    }
}

void LCA::addEdge(int father, int son) {
    doEncode(father);
    doEncode(son);
    tree[encode[father]][encode[son]] = true;
    parent[encode[son]] = encode[father];
}

void LCA::depthFirstSearch(int current, int depth) {
    // marking first appearance for current node
    if(firstAppearance[current] == -1) {
        firstAppearance[current] = euler.size();
    }
    // pushing root to euler walk
    euler.push_back(current);
    // pushing depth of current node
    this->depth.push_back(depth);
    for(int i = 0; i < vertices; i++) {
        if(tree[current][i] == true) {
            depthFirstSearch(i, depth + 1);
            euler.push_back(current);
            this->depth.push_back(depth);
        }
    }
}

void LCA::doEulerWalk() {
    int root = getRoot();
    depthFirstSearch(root, 0);
    preProcessRMQ();
}

int LCA::getRoot() {
    for(int i = 0; i < vertices; i++) {
        if(parent[i] == -1) {
            return i;
        }
    }
    return -1;
}

// <O(N logN) Preprocessing time, O(1) Query time>
void LCA::preProcessRMQ() {
    M = new int*[depth.size()];
    int logDepth = log2(depth.size());
    for(unsigned int i = 0; i < depth.size(); i++) {
        M[i] = new int[logDepth];
    }

    //initialize M for the intervals with length 1
    for(unsigned int i = 0; i < depth.size(); i++) {
        M[i][0] = i;
    }
    //compute values from smaller to bigger intervals
    for(unsigned int j = 1; 1 << j <= depth.size(); j++) {
        for(unsigned int i = 0; i + (1 << j) - 1 < depth.size(); i++) {
            if(depth[M[i][j - 1]] < depth[M[i + (1 << (j - 1))][j - 1]]) {
                M[i][j] = M[i][j - 1];
            } else {
                M[i][j] = M[i + (1 << (j - 1))][j - 1];
            }
        }
    }
}

int LCA::queryRMQ(int i, int j) {
    if(i > j) {
        swap(i, j);
    }

    int k = log2(j - i + 1);

    if(depth[M[i][k]] <= depth[M[j - (1 << k) + 1][k]]) {
        return M[i][k];
    } else {
        return M[j - (1 << k) + 1][k];
    }
}

int LCA::getLCA(int u, int v) {
    u = encode[u];
    v = encode[v];
	// trivial case
	if (u == v) {
        return u;
	}

	if(firstAppearance[u] > firstAppearance[v]) {
        swap(u, v);
	}

	// doing RMQ in the required range
	return decode[euler[queryRMQ(firstAppearance[u], firstAppearance[v])]];
}

int getLCA(LCA &tree, vector<int> &taxIds) {
    int lca;
    if(taxIds.size() >= 2) {
        lca = tree.getLCA(taxIds[0], taxIds[1]);
        for(int i = 2; i < taxIds.size(); i++) {
            lca = tree.getLCA(lca, taxIds[i]);
        }
        return lca;
    } else {
        return taxIds[0];
    }
}

int main(int argc, char *argv[]) {
    if(argc == 3) {
        int father, son;

        ifstream vertices(argv[1]);
        // counts the number of distinct vertices in the tree
        set<int> distinctVertices;
        while(vertices >> father >> son) {
            distinctVertices.insert(father);
            distinctVertices.insert(son);
        }


        LCA lca(distinctVertices.size());
        ifstream tree(argv[1]);
        while(tree >> father >> son) {
            lca.addEdge(father, son);
        }

        lca.doEulerWalk();

        string currentRead, nextRead;
        int currentTaxId, nextTaxId;
        int currentKmer, nextKmer;
        vector<int> taxIds;
        ifstream queries(argv[2]);
        queries >> currentRead >> currentTaxId >> currentKmer;
        taxIds.push_back(currentTaxId);
        while(queries >> nextRead >> nextTaxId >> nextKmer) {
            if(nextRead.compare(currentRead) == 0) {
                taxIds.push_back(nextTaxId);
            } else {
                cout << currentRead << "\t";
                cout << getLCA(lca, taxIds) << "\n";
                taxIds.clear();
                currentRead = nextRead;
                taxIds.push_back(nextTaxId);
            }
        }
        cout << currentRead << "\t";
        cout << getLCA(lca, taxIds) << "\n";
    } else {
        cout << "Usage: " << argv[0] << " tree.tsv queries.tsv\n";
    }
}
