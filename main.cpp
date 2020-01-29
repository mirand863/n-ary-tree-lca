#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

class LCA {
    private:
        unordered_map<string, vector<string>> parents;
        vector<int> euler;
        vector<int> depth;
        int *firstAppearance;
        int vertices;
        unordered_map<string, int> encode;
        unordered_map<int, string> decode;
        int **M;
        void depthFirstSearch(string current, int depth);
        void preProcessRMQ();
        int queryRMQ(int i, int j);
    public:
        LCA();
        virtual ~LCA();
        void addEdge(string father, string son);
        void doEulerWalk(string root);
        int getLCA(int u, int v);
        string getLCA(vector<string> &taxIds);
};

LCA::LCA() {
    vertices = 0;
}

LCA::~LCA() {
    delete [] firstAppearance;
    for(unsigned int i = 0; i < depth.size(); i++) {
        delete [] M[i];
    }
    delete [] M;
}

void LCA::addEdge(string father, string son) {
    if(encode.count(father) == 0) {
        encode.insert({father, vertices});
        decode.insert({vertices, father});
//        firstAppearance[vertices] = -1;
        vertices++;
    }
    if(encode.count(son) == 0) {
        encode.insert({son, vertices});
        decode.insert({vertices, son});
//        firstAppearance[vertices] = -1;
        vertices++;
    }
    if(parents.count(father) == 0) {
        vector<string> children;
        children.push_back(son);
        parents[father] = children;
    } else {
        parents[father].push_back(son);
    }
}

void LCA::depthFirstSearch(string current, int depth) {
    // marking first appearance for current node
    if(firstAppearance[encode[current]] == -1) {
        firstAppearance[encode[current]] = euler.size();
    }
    // pushing root to euler walk
    euler.push_back(encode[current]);
    // pushing depth of current node
    this->depth.push_back(depth);
    for(unsigned int i = 0; i < parents[current].size(); i++) {
        depthFirstSearch(parents[current][i], depth + 1);
        euler.push_back(encode[current]);
        this->depth.push_back(depth);
    }
}

void LCA::doEulerWalk(string root) {
    firstAppearance = new int[vertices];
    for(int i = 0; i < vertices; i++) {
        firstAppearance[i] = -1;
    }
    depthFirstSearch(root, 0);
    preProcessRMQ();
}

// <O(N logN) Preprocessing time, O(1) Query time>
void LCA::preProcessRMQ() {
    M = new int*[depth.size()];
    int logDepth = ceil(log2(depth.size()));
    for(unsigned int i = 0; i < depth.size(); i++) {
        M[i] = new int[logDepth];
        M[i][0] = i; //initialize M for the intervals with length 1
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
    // trivial case
    if (u == v) {
        return u;
    }

    if(firstAppearance[u] > firstAppearance[v]) {
        swap(u, v);
    }

    // doing RMQ in the required range
    return euler[queryRMQ(firstAppearance[u], firstAppearance[v])];
}

string LCA::getLCA(vector<string> &taxIds) {
    int lca;
    if(taxIds.size() >= 2) {
        lca = getLCA(encode[taxIds[0]], encode[taxIds[1]]);
        for(int i = 2; i < taxIds.size(); i++) {
            lca = getLCA(lca, encode[taxIds[i]]);
        }
        return decode[lca];
    } else {
        return taxIds[0];
    }
}

int main(int argc, char *argv[]) {
    if(argc == 4) {
        string father, son, root;
        LCA lca;
        ifstream tree(argv[1]);
        while(tree >> father >> son) {
            lca.addEdge(father, son);
        }

        lca.doEulerWalk(argv[3]);

        string currentQuery, nextQuery;
        string currentNode, nextNode;
        vector<string> nodes;
        ifstream queries(argv[2]);
        queries >> currentQuery >> currentNode;
        nodes.push_back(currentNode);
        while(queries >> nextQuery >> nextNode) {
            if(nextQuery.compare(currentQuery) == 0) {
                nodes.push_back(nextNode);
            } else {
                cout << currentQuery << "\t";
                cout << lca.getLCA(nodes) << "\n";
                nodes.clear();
                currentQuery = nextQuery;
                nodes.push_back(nextNode);
            }
        }
        cout << currentQuery << "\t";
        cout << lca.getLCA(nodes) << "\n";
    } else {
        cout << "Usage: " << argv[0] << " tree.tsv queries.tsv root\n";
    }
}
