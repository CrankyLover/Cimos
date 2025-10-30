#ifndef MATCHING_CIMOS
#define MATCHING_CIMOS

#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

#include "graph/graph.h"
#include "matching/matching.h"

struct TreeNode {
    std::vector<uint> forwards_;
    std::vector<uint> backwards_;
    std::vector<uint> neighbors_;
};

struct ExtendableVertex {
    uint E;
    uint matched_nbrs;
    uint u_min;
    uint index_min;
    uint label_min;

    ExtendableVertex()
    : E(NOT_EXIST), matched_nbrs(0u), u_min(NOT_EXIST), index_min(NOT_EXIST), label_min(NOT_EXIST) {}
    ExtendableVertex(uint E_arg, uint matched_nbrs_arg, uint u_min_arg, uint index_min_arg, uint label_min_arg)
    : E(E_arg), matched_nbrs(matched_nbrs_arg), u_min(u_min_arg), index_min(index_min_arg), label_min(label_min_arg) {}
};

struct Node {
    uint label;
    // below two are matched in the same position
    std::unordered_map<uint, uint> queryToVertex;
    std::unordered_set<uint> queries;
    std::unordered_map<uint, std::vector<uint>> truncatedQueries;
    std::unordered_set<uint> endQueries;

    std::vector<std::pair<uint, uint>> topology;

};

struct MatchingTree {

   std::vector<std::vector<Node*>> tree;
   std::vector<std::vector<std::vector<uint>>> children;
   // each matching tree is a vector<vector> e.g. [[2,3,4], [5], [], [], []], which represents a tree in vector
   std::unordered_map<uint, std::vector<uint>> paths;
   // store each path by the index of children for each query graph

   std::unordered_map<uint, uint> QueryToDepth;

   bool is_reverse;
   std::unordered_map<std::pair<uint, uint>, bool, PairHash> isTreeNodeStoresIR;
};


class QueryDAG {
private:
    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_;
    std::vector<uint> serialized_tree_;
public:
    void ModifyEdgeIndex(Graph& queryGraph);
    void BuildDAG(Graph& dataGraph, Graph& queryGraph);

    std::vector<std::vector<uint>>& GetEdgeIndex() {
        return eidx_;
    }

    std::vector<TreeNode>& GetTreeNodes() {
        return treeNode_;
    }

    uint& GetRootNode() {
        return q_root_;
    }

    std::vector<uint>& GetSerializedTree() {
        return serialized_tree_;
    }
};


class Cimos : public matching
{
private:

    std::vector<QueryDAG> DAGs;

    std::vector<std::unordered_map<uint, bool>> F;
    std::vector<std::unordered_map<uint, bool>> B;
    std::vector<std::unordered_map<uint, std::unordered_map<uint, uint>>> F1;
    std::vector<std::unordered_map<uint, uint>> FP;
    std::vector<std::unordered_map<uint, std::unordered_map<uint, uint>>> Num;
    std::vector<std::unordered_map<uint, uint>> BC;
    std::vector<uint> D;

    std::vector<std::unordered_map<uint, std::unordered_set<uint>>> EdgeIndex;

    std::unordered_map<std::pair<uint, uint>, std::vector<std::vector<uint>>, PairHash> intermediateResults;
    std::unordered_map<std::pair<uint, uint>, double, PairHash> globalN;

    std::queue<std::pair<uint, uint>> Q1;
    std::queue<std::pair<uint, uint>> Q2;

    std::set<std::tuple<uint, uint, uint>> all_unique_edges;
    std::vector<std::tuple<uint, uint, uint>> all_vector_unique_edges;
    // (label1, label2) to store all unique edges in query graph set

    std::vector<std::unordered_map<uint, std::vector<std::pair<uint, uint>>>> edge_in_graphs;

    std::vector<MatchingTree> matching_tree_set;
    std::vector<std::unordered_map<uint, uint>> global_hash;
    std::unordered_map<std::pair<uint, uint>, uint, PairHash> vertexMapping;

    std::vector<uint> queryResults;

public:

    Cimos(std::vector<Graph>& query_graphs, Graph& data_graph);
    ~Cimos() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    
    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2, uint label) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void LoadAllEdges(std::vector<std::string>& paths);
    void CalculateEdgeMapping();
    void ConstructMatchingTrees();
    void PrintQueryPlan();
    void PrintQueryReuslts();
    void CalculateCandidateSize();
    void SampleEdges(uint label1, uint label2, uint sampleNumber);
    void SamplingSubgraphs(uint i, uint j, uint treeIndex, MatchingTree& tree, std::vector<uint>& path);
    uint CalculateMinimalTrie(std::vector<uint>& maxLabels, std::vector<std::vector<std::pair<uint, uint>>>& maxVector, uint i, uint j);
    void CalculateOverallStatistic();
    void PrintTimeOfEachPhrase();
    std::chrono::high_resolution_clock::time_point vs;

private:
    void BuildDAG();
    void BuildIndex();

    void CalculateF();
    void CalculateB();
    void PrintFAndB();
    void PrintEdgeIndex();

    MatchingTree ConstructOneTree(std::vector<std::pair<uint, uint>> &maps, std::vector<uint> &ids, uint index, uint sampleNumber);
    uint GetMaxTreeSize();
    
    void InsertionModifyF(uint qnum, uint u, uint u_c, uint v, uint v_c);
    void InsertionModifyB(uint qnum, uint u, uint u_p, uint v, uint v_p);
    void DeletionModifyF(uint qnum, uint u, uint u_c, uint v, uint v_c);
    void DeletionModifyB(uint qnum, uint u, uint u_p, uint v, uint v_p);

    void StaticSearchMatches(uint depth, uint treeIndex, uint position, std::vector<uint>& m);
    void DynamicSearchMatches(uint depth, uint treeIndex, uint q, std::vector<uint>& m, std::vector<uint>& queryV, std::vector<ExtendableVertex>& dynamicOrders);
};

struct VectorHash {
    size_t operator () (const std::vector<std::pair<uint, uint>>& v) const {
        size_t seed = 0;
        for (const auto& p : v) {
            seed ^= pair_hash{}(p) + 0x9e3779b9 + (seed << 6) + (seed >> 2);  
        }
        return seed;
    }
};

struct VectorIntHash {
    size_t operator()(const std::vector<uint>& vec) const {
        std::hash<uint> hasher;
        size_t seed = 0;
        for (uint i : vec) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

struct tuple_hash {
    template <typename T1, typename T2, typename T3>
    size_t operator () (const std::tuple<T1, T2, T3>& p) const {
        size_t seed = 0;
        auto h1 = std::hash<T1>{}(std::get<0>(p));  
        auto h2 = std::hash<T2>{}(std::get<1>(p)); 
        auto h3 = std::hash<T3>{}(std::get<2>(p));
        seed ^= (h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        seed ^= (h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        seed ^= (h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        return seed;  
    }
};

struct VectorTupleHash {
    size_t operator () (const std::vector<std::tuple<uint, uint, uint>>& v) const {
        size_t seed = 0;
        for (const auto& p : v) {
            seed ^= tuple_hash{}(p) + 0x9e3779b9 + (seed << 6) + (seed >> 2); 
        }
        return seed;
    }
};

struct KeyHash {
    size_t operator()(const std::tuple<char, uint, uint, uint>& k) const {
        size_t seed = 0;
        seed ^= std::hash<char>()(std::get<0>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= std::hash<uint>()(std::get<1>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= std::hash<uint>()(std::get<2>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= std::hash<uint>()(std::get<3>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};


class UnionFind {
public:

    UnionFind() {}

    uint find(uint x) {
        if (parent.find(x) == parent.end()) {
            parent[x] = x; 
        }
        if (parent[x] != x) {
            parent[x] = find(parent[x]); 
        }
        return parent[x];
    }

    void unionSets(uint x, uint y) {
        uint rootX = find(x);
        uint rootY = find(y);
        if (rootX != rootY) {
            parent[rootY] = rootX; 
        }
    }

private:
    std::unordered_map<uint, uint> parent; 
};

struct CompareTree {
    bool operator()(const MatchingTree& mtr1, const MatchingTree& mtr2) {
        uint count1 = 0, count2 = 0;
        for (const auto& layer: mtr1.tree) {
            count1 += layer.size();
        }
        for (const auto& layer: mtr2.tree) {
            count2 += layer.size();
        }
        return count1 > count2;
    }
};

#endif 
