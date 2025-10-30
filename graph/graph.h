#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
#include <vector>
#include <unordered_map>
#include "utils/types.h"
#include "utils/utils.h"

struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

class Graph
{
public:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;

public:
    std::queue<InsertUnit> updates_;
    std::vector<uint> vlabels_;
    std::vector<std::tuple<uint, uint, uint, uint, uint>> edges;
    std::unordered_map<std::pair<uint, uint>, uint, PairHash> vertexPairToEdgeLabel;
public:
    Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    uint GetDegree(uint v) const;

    void LoadFromFile(const std::string &path);
    void LoadFromFile(const std::string &path, bool isQuery);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
};

#endif //GRAPH_GRAPH
