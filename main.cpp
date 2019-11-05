#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <cmath>
#include <set>

using namespace std;

class Tree {
    private:
        bool **tree;
        int *parent;
        int vertices;
        vector<int> euler;
        vector<int> depth;
        int *firstAppearance;
        int **M;
        int getRoot();
        void depthFirstSearch(int current, int depth);
        void preProcessRMQ();
        int queryRMQ(int i, int j);
    public:
        Tree(int vertices);
        virtual ~Tree();
        void addEdge(int father, int son);
        void doEulerWalk();
        int getLCA(int u, int v);
};

Tree::Tree(int vertices) {
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

Tree::~Tree() {
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

void Tree::addEdge(int father, int son) {
    tree[father - 1][son - 1] = true;
    parent[son - 1] = father - 1;
}

void Tree::depthFirstSearch(int current, int depth) {
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

void Tree::doEulerWalk() {
    int root = getRoot();
    depthFirstSearch(root, 0);
    preProcessRMQ();
}

int Tree::getRoot() {
    for(int i = 0; i < vertices; i++) {
        if(parent[i] == -1) {
            return i;
        }
    }
    return -1;
}

// <O(N logN) Preprocessing time, O(1) Query time>
void Tree::preProcessRMQ() {
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

int Tree::queryRMQ(int i, int j) {
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

int Tree::getLCA(int u, int v) {
    u--; v--;
	// trivial case
	if (u == v) {
        return u + 1;
	}

	if(firstAppearance[u] > firstAppearance[v]) {
        swap(u, v);
	}

	// doing RMQ in the required range
	return euler[queryRMQ(firstAppearance[u], firstAppearance[v])] + 1;
}

int getLCA(Tree &tree, vector<int> &taxIds, unordered_map<int, int> &encode, unordered_map<int, int> &decode) {
    int lca;
    if(taxIds.size() >= 2) {
        lca = tree.getLCA(encode[taxIds[0]], encode[taxIds[1]]);
        lca = decode[lca];
        for(int i = 2; i < taxIds.size(); i++) {
            lca = tree.getLCA(encode[lca], encode[taxIds[i]]);
            lca = decode[lca];
        }
        return lca;
    } else {
        return taxIds[0];
    }
}

int main(int argc, char *argv[]) {
    if(argc == 3) {
        int u, v;
        ifstream tree(argv[1]);
        ifstream vertices(argv[1]);
        set<int> distinctVertices;
        while(vertices >> u >> v) {
            distinctVertices.insert(u);
            distinctVertices.insert(v);
        }
        Tree lca(distinctVertices.size());
        ifstream queries(argv[2]);
        unordered_map<int, int> encode, decode;
        int counter = 1;
        while(tree >> u >> v) {
            if(encode.count(u) == 0) {
                encode.insert({u, counter});
                decode.insert({counter, u});
                counter++;
            }
            if(encode.count(v) == 0) {
                encode.insert({v, counter});
                decode.insert({counter, v});
                counter++;
            }
            lca.addEdge(encode[u], encode[v]);
        }
        lca.doEulerWalk();
        string currentRead, nextRead;
        int currentTaxId, nextTaxId;
        int currentKmer, nextKmer;
        vector<int> taxIds;
        queries >> currentRead >> currentTaxId >> currentKmer;
        taxIds.push_back(currentTaxId);
        while(queries >> nextRead >> nextTaxId >> nextKmer) {
            if(nextRead.compare(currentRead) == 0) {
                taxIds.push_back(nextTaxId);
            } else {
                cout << currentRead << "\t";
                cout << getLCA(lca, taxIds, encode, decode) << "\n";
                taxIds.clear();
                currentRead = nextRead;
                taxIds.push_back(nextTaxId);
            }
        }
        cout << currentRead << "\t";
        cout << getLCA(lca, taxIds, encode, decode) << "\n";
    } else {
        cout << "Usage: " << argv[0] << " tree.tsv queries.tsv\n";
    }
}
