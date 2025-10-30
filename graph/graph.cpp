#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include "utils/types.h"
#include "utils/utils.h"
#include "graph/graph.h"

Graph::Graph()
: edge_count_(0)
, vlabel_count_(0)
, elabel_count_(0)
, neighbors_{}
, updates_{}
, vlabels_{}
{}

void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    
    vlabel_count_ = std::max(vlabel_count_, label + 1);
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
}

void Graph::AddEdge(uint v1, uint v2, uint label)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;
    
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].insert(lower, v2);
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);

    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
}

void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);
   
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);

    edge_count_--;
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)
    if (GetVertexLabel(i) != NOT_EXIST)
    {
        std::queue<uint> bfs_queue;
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}

void Graph::LoadFromFile(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    uint counterV = 0;
    uint counterE = 0;
    while (ifs >> type)
    {
        if (type == 't')
        {
            uint temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label, unknown;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
            counterV ++;
            if (counterV % 1000000 == 0) std::cout << counterV << " vertices have been loaded\n";
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id ;
            label = 0;
            AddEdge(from_id, to_id, label);
            counterE ++;
            if (counterE % 1000000 == 0) std::cout << counterE << " edges have been loaded\n";
        }
    }
    ifs.close();
}

void Graph::LoadFromFile(const std::string &path, bool isQuery)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label, unk;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label = 0;
            ifs >> from_id >> to_id;
            label = 0;
            AddEdge(from_id, to_id, label);
            uint label1 = GetVertexLabel(from_id);
            uint label2 = GetVertexLabel(to_id);
            if (label1 > label2) {
                std::swap(label1, label2);
                std::swap(from_id, to_id);
            }
            this->edges.emplace_back(GetVertexLabel(from_id), GetVertexLabel(to_id), label, from_id, to_id);
            this->vertexPairToEdgeLabel[std::make_pair(from_id, to_id)] = label;
            this->vertexPairToEdgeLabel[std::make_pair(to_id, from_id)] = label;
        }
    }
    ifs.close();
}

void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id ;
            updates_.emplace('e', type == "e", from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() <<
        "\n# edges = " << NumEdges() << std::endl;
    //         // print graph
    // std::cout << "labels: ";
    // for (uint i = 0; i < vlabels_.size(); i++)
    // {
    //     std::cout << i << ":" << vlabels_[i] << " (";
    //     for (uint j = 0; j < neighbors_[i].size(); j++)
    //     {
    //         std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
    //     }
    //     std::cout << ")" << std::endl;
    // }
}