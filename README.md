# n-ary-tree-lca
Based on the implementation provided at https://www.geeksforgeeks.org/lca-n-ary-tree-constant-query-o1/

# Compilation
```
git clone https://github.com/mirand863/n-ary-tree-lca.git lca
cd lca
g++ -o lca main.cpp
```

# Input Files
## Tree
The input for the tree should be in the tsv format, where the first column is the father node and the second coluns contains the son. For example, supposing we have the given tree:
```
   a
 / | \
b  c  d
 / | \
e  f  g
```
The input for the tree would be the following tsv file:
```
a   b
a   c
a   d
c   e
c   f
c   g
```
## Queries
The input for the queries should also be in the tsv format. However, the first column is the query id and the second column is the node to be queried. For example, the queries below
```
LCA(a, b)
LCA(e, f, g)
LCA(g, b)
```
should be in the following format:
```
Query0	a
Query0	b
Query1	e
Query1	f
Query1	g
Query2	g
Query2	b
```
# Usage:
```
./lca tree.tsv queries.tsv root
e.g.: ./lca test_tree.tsv test_queries.tsv a
```
