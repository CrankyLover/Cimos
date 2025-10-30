#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
#include <bitset>
#include <fstream>
#include <random>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/cimos.h"

uint paddingNumber = 3e3;

Cimos::Cimos(std::vector<Graph>& query_graphs, Graph& data_graph)
: matching(
    query_graphs, data_graph)
    , DAGs(query_graphs.size())
    , F(data_.NumVertices() + paddingNumber)
    , B(data_.NumVertices() + paddingNumber)
    , F1(data_.NumVertices() + paddingNumber)
    , FP(data_.NumVertices() + paddingNumber)
    , Num(data_.NumVertices() + paddingNumber)
    , BC(data_.NumVertices() + paddingNumber)
    , D(data_.NumVertices() + paddingNumber)
{
    // for (uint i = 0; i < query_graphs.size(); i++) {
    //     DAGs[i].ModifyEdgeIndex(query_graphs[i]);
    // }
    queryResults.resize(query_graphs.size());
}

void Cimos::CalculateCandidateSize() {
    uint vertexNumber = data_.NumVertices();
    uint count = 0;
    for (uint i = 0; i < vertexNumber; i++) {
        if (D[i] > 0) {
            count ++;
        }
    }
    std::cout << count << "\n";
    std::cout << vertexNumber << "\n";
    std::cout << "Candidate ratio: " << ((1.0 * count) / (1.0 * vertexNumber)) * 100 << "\n";
}

void QueryDAG::BuildDAG(Graph& dataGraph, Graph& queryGraph) {

    q_root_ = 0;

    // build a spaning tree with the BFS order
    std::vector<bool> visited(queryGraph.NumVertices(), false);
    std::vector<bool> exist_on_tree(queryGraph.NumVertices(), false);
    serialized_tree_.resize(queryGraph.NumVertices());
    treeNode_.resize(queryGraph.NumVertices());
    std::queue<uint> bfs_queue;

    bfs_queue.push(q_root_);
    exist_on_tree[q_root_] = true;
    uint order_pos = 0u;
    serialized_tree_[order_pos++] = q_root_;

    while (!bfs_queue.empty()) {
        const uint u = bfs_queue.front();
        bfs_queue.pop();
        visited[u] = true;

        const auto& q_nbrs = queryGraph.GetNeighbors(u);

        for (uint j = 0; j < q_nbrs.size(); j++) {
            const uint u_other = q_nbrs[j];
            
            if (!exist_on_tree[u_other]) {
                bfs_queue.push(u_other);
                exist_on_tree[u_other] = true;
                serialized_tree_[order_pos++] = u_other;
            }
            if (!visited[u_other]) {
                treeNode_[u].forwards_.push_back(u_other);
            }
            else {
                treeNode_[u].backwards_.push_back(u_other);
            }
            treeNode_[u].neighbors_.push_back(u_other);
        }
    }  

        // std::cout << "DAG: " << std::endl;
        // for (uint i = 0; i < queryGraph.NumVertices(); ++i)
        // {
        //     std::cout << serialized_tree_[i] << ": (backwards: ";
        //     for (auto j: treeNode_[serialized_tree_[i]].backwards_)
        //         std::cout << j << " ";

        //     std::cout << ") (forwards: ";
        //     for (auto j: treeNode_[serialized_tree_[i]].forwards_)
        //         std::cout << j << " ";

        //     std::cout << ")" << std::endl;
        // }
        // std::cout << std::endl;
        // std::cout << std::endl;
        // std::cout << std::endl;

}

void Cimos::Preprocessing() {
    BuildDAG();
    std::cout << "Finish Building DAGs\n";
    BuildIndex();
}

void Cimos::BuildDAG() {
    std::cout << "Starting Building DAG\n";
    for (uint i = 0; i < querys_.size(); i++) {
        DAGs[i].BuildDAG(data_, querys_[i]);
    }
}

void Cimos::BuildIndex() {
    EdgeIndex.resize(data_.NumVertices() + paddingNumber);
    CalculateF();
    std::cout << "Finish Calculating F\n";
    CalculateB();
    std::cout << "Finish Calculating B\n";
    // for (uint i = 0; i < data_.NumVertices(); i++) {
    //     std::cout << "\nThe Vertex in data graph is V" << i << "\n"; 
    //     for (const auto& [q, f]: F[i]) {
    //         std::cout << "Query 's vertex " << q << " 's F is " << f << "\n";
    //     }
    // }
    // for (uint i = 0; i < data_.NumVertices(); i++) {
    //     std::cout << "\nThe Vertex in data graph is V" << i << "\n"; 
    //     for (const auto& [q, b]: B[i]) {
    //         std::cout << "Query 's vertex " << q << " 's B is " << b << "\n";
    //     }
    // }
}

void Cimos::PrintFAndB() {  

}

void Cimos::PrintEdgeIndex() {

}

void Cimos::CalculateF() {
    for (uint qnum = 0; qnum < querys_.size(); qnum++) {
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint i = 0; i < query.NumVertices(); i++) {
            uint u = qDAG.GetSerializedTree()[i];
            const uint query_label = query.GetVertexLabel(u);
            for (size_t j = 0; j < data_.NumVertices(); j++) {
                if (data_.GetVertexLabel(j) == query_label) {
                    for (uint k = 0; k < qDAG.GetTreeNodes()[u].backwards_.size(); k++) {
                        const uint u_other = qDAG.GetTreeNodes()[u].backwards_[k];

                        const auto& d_nbrs = data_.GetNeighbors(j);

                        for (uint m = 0; m < d_nbrs.size(); m++) {
                            const uint v_other = d_nbrs[m];

                            if (query.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other)) {

                                if (!EdgeIndex[j][query.GetVertexLabel(u_other)].count(v_other)) {
                                    EdgeIndex[j][query.GetVertexLabel(u_other)].insert(v_other);
                                }

                                // above 
                                if (F[v_other][vertexMapping[std::make_pair(qnum, u_other)]] > 0) {
                                    F1[j][vertexMapping[std::make_pair(qnum, u)]][u_other] += 1;
                                }
                            }
                        }
                        if (F1[j][vertexMapping[std::make_pair(qnum, u)]][u_other]) {
                            FP[j][vertexMapping[std::make_pair(qnum, u)]] += 1;
                        }
                    }

                    if (FP[j][vertexMapping[std::make_pair(qnum, u)]] == qDAG.GetTreeNodes()[u].backwards_.size()) {
                        F[j][vertexMapping[std::make_pair(qnum, u)]] = 1;
                    } else {
                        F[j][vertexMapping[std::make_pair(qnum, u)]] = 0;
                    }
                }
            }
        }
    }
}

void Cimos::CalculateB() {
    for (uint qnum = 0; qnum < querys_.size(); qnum++) {
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint i = 0; i < query.NumVertices(); i++) {
            const uint u = qDAG.GetSerializedTree()[query.NumVertices() - 1 - i];
            const uint query_label = query.GetVertexLabel(u);
            
            for (size_t j = 0; j < data_.NumVertices(); j++) {
                if (data_.GetVertexLabel(j) == query_label) {
                    for (uint k = 0; k < qDAG.GetTreeNodes()[u].forwards_.size(); k++) {
                        const uint u_other = qDAG.GetTreeNodes()[u].forwards_[k];

                        const auto& d_nbrs = data_.GetNeighbors(j);

                        for (uint m = 0; m < d_nbrs.size(); m++) {
                            const uint v_other = d_nbrs[m];
                            if (query.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other)) {
                                // edge index modify
                                if (!EdgeIndex[j][query.GetVertexLabel(u_other)].count(v_other)) {
                                    EdgeIndex[j][query.GetVertexLabel(u_other)].insert(v_other);
                                }                                
                                // finish
                                if (B[v_other][vertexMapping[std::make_pair(qnum, u_other)]] > 0) {
                                    Num[j][vertexMapping[std::make_pair(qnum, u)]][u_other] += 1;
                                }
                            }
                        }
                        if (Num[j][vertexMapping[std::make_pair(qnum, u)]][u_other]) {
                            BC[j][vertexMapping[std::make_pair(qnum, u)]] += 1;
                        }
                    }
                    if (BC[j][vertexMapping[std::make_pair(qnum, u)]] == qDAG.GetTreeNodes()[u].forwards_.size() && F[j][vertexMapping[std::make_pair(qnum, u)]] == 1) {
                        B[j][vertexMapping[std::make_pair(qnum, u)]] = 1;
                        D[j] ++;
                    } else {
                        B[j][vertexMapping[std::make_pair(qnum, u)]] = 0;
                    }
                }
            }
        }
    }

    for (uint qnum = 0; qnum < querys_.size(); qnum++) {
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint i = 0; i < query.NumVertices(); i++) {
            const uint u = qDAG.GetSerializedTree()[i];
            const uint query_label = query.GetVertexLabel(u);

            for (size_t j = 0; j < data_.NumVertices(); j++) {
                if (data_.GetVertexLabel(j) == query_label) {
                    for (uint k = 0; k < qDAG.GetTreeNodes()[u].backwards_.size(); k++) {

                        const uint u_other = qDAG.GetTreeNodes()[u].backwards_[k];

                        const auto& d_nbrs = data_.GetNeighbors(j);

                        for (uint m = 0; m < d_nbrs.size(); m++) {
                            const uint v_other = d_nbrs[m];
                            if (query.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other)) {

                                if (B[v_other][vertexMapping[std::make_pair(qnum, u_other)]] > 0) {
                                    Num[j][vertexMapping[std::make_pair(qnum, u)]][u_other] += 1;
                                }
                            }
                        }
                    }
                }
            }

        }
    }
}

void Cimos::InitialMatching()
{

}

void Cimos::InsertionModifyF(uint qnum, uint u, uint u_c, uint v, uint v_c)
{
    QueryDAG& qDAG = DAGs[qnum];

    if (F1[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] == 0) {
        FP[v_c][vertexMapping[std::make_pair(qnum, u_c)]] += 1;
        if (FP[v_c][vertexMapping[std::make_pair(qnum, u_c)]] == qDAG.GetTreeNodes()[u_c].backwards_.size()) {
            F[v_c][vertexMapping[std::make_pair(qnum, u_c)]] = 1;
            Q1.emplace(v_c, u_c);
            if (BC[v_c][vertexMapping[std::make_pair(qnum, u_c)]] == qDAG.GetTreeNodes()[u_c].forwards_.size()) {
                B[v_c][vertexMapping[std::make_pair(qnum, u_c)]] = 1;
                D[v_c] ++;
                Q2.emplace(v_c, u_c);
            }
        }
    }

    F1[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] += 1;
}

void Cimos::InsertionModifyB(uint qnum, uint u, uint u_p, uint v, uint v_p)
{
    QueryDAG& qDAG = DAGs[qnum];

    if (Num[v_p][vertexMapping[std::make_pair(qnum, u_p)]][u] == 0) {
        BC[v_p][vertexMapping[std::make_pair(qnum, u_p)]] += 1;
        if (F[v_p][vertexMapping[std::make_pair(qnum, u_p)]] && BC[v_p][vertexMapping[std::make_pair(qnum, u_p)]] == qDAG.GetTreeNodes()[u_p].forwards_.size()) {
            B[v_p][vertexMapping[std::make_pair(qnum, u_p)]] = 1;
            D[v_p] ++;
            Q2.emplace(v_p, u_p);
        }
    }

    Num[v_p][vertexMapping[std::make_pair(qnum, u_p)]][u] += 1;
}

void Cimos::DeletionModifyF(uint qnum, uint u, uint u_c, uint v, uint v_c)
{
    F1[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] -= 1;
    if (F1[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] == 0) {
        if (F[v_c][vertexMapping[std::make_pair(qnum, u_c)]] == 1) {
            if (B[v_c][vertexMapping[std::make_pair(qnum, u_c)]] == 1) {
                Q2.emplace(v_c, u_c);
                B[v_c][vertexMapping[std::make_pair(qnum, u_c)]] = 0;
                D[v_c] --;
            }
            Q1.emplace(v_c, u_c);
            F[v_c][vertexMapping[std::make_pair(qnum, u_c)]] = 0;
        }
        FP[v_c][vertexMapping[std::make_pair(qnum, u_c)]] -= 1;
    }
}

void Cimos::DeletionModifyB(uint qnum, uint u, uint u_p, uint v, uint v_p)
{
    Num[v_p][vertexMapping[std::make_pair(qnum, u_p)]][u] -= 1;
    if (Num[v_p][vertexMapping[std::make_pair(qnum, u_p)]][u] == 0) {
        if (B[v_p][vertexMapping[std::make_pair(qnum, u_p)]] == 1) {
            Q2.emplace(v_p, u_p);
            B[v_p][vertexMapping[std::make_pair(qnum, u_p)]] = 0;
            D[v_p] --;
        }
        BC[v_p][vertexMapping[std::make_pair(qnum, u_p)]] -= 1;
    }
}

void Cimos::StaticSearchMatches(uint depth, uint treeIndex, uint position, std::vector<uint>& m) {

    // get the fully matched result for some queries
    for (const auto& q: matching_tree_set[treeIndex].tree[depth][position]->endQueries) {
        queryResults[global_hash[treeIndex][q]] ++;
    }

    // handle static matching order part
    for (const auto& childrenPosition: matching_tree_set[treeIndex].children[depth][position]) {
        if (childrenPosition == 4294967295) {
            continue;
        }
        Node* node = matching_tree_set[treeIndex].tree[depth + 1][childrenPosition];
        uint targetLabel = node->label;
        uint minValue = NOT_EXIST;
        uint minIndex;
        uint exampleQ = node->queryToVertex.begin()->first;
        uint exampleU = node->queryToVertex.begin()->second;

        // find the minimal Num first (previousU)
        for (const auto& tp: node->topology) {

            uint previousV = m[tp.first];
            uint previousPosition = matching_tree_set[treeIndex].paths[exampleQ][tp.first]; 
            uint previousU = matching_tree_set[treeIndex].tree[tp.first][previousPosition]->queryToVertex[exampleQ];
            uint qu = vertexMapping[std::make_pair(global_hash[treeIndex][exampleQ], previousU)];
            uint candidatesSize = Num[previousV][qu][exampleU];
            if (candidatesSize < minValue) {
                minValue = candidatesSize;
                minIndex = tp.first;
            }
        }

        uint minV = m[minIndex];

        for (const auto& v: EdgeIndex[minV][targetLabel]) {
            // prune rule 1 & 2 : visited and index
            if (visited_[v] || !D[v]) {
                continue;
            }
            // prune rule 3：joinable
            bool isJoinable = true;
            for (const auto& tp: node->topology) {
                if (tp.first == minIndex) {
                    continue;
                } else {
                    if (!EdgeIndex[m[tp.first]][targetLabel].count(v)) {
                        isJoinable = false;
                        break;
                    }
                }
            }
            if (!isJoinable) {
                continue;
            }

            m[depth + 1] = v;
            visited_[v] = true;
            StaticSearchMatches(depth + 1, treeIndex, childrenPosition, m);
            visited_[v] = false;
            m[depth + 1] = UNMATCHED;
        }
    }

    // handle dynamic matching order part
    for (const auto& [q, qV]: matching_tree_set[treeIndex].tree[depth][position]->truncatedQueries) {
        std::vector<ExtendableVertex> dynamicOrders(querys_[global_hash[treeIndex][q]].NumVertices());
        std::vector<uint> queryV = qV;
        for (uint i = 0; i <= depth; i++) {
            uint u = matching_tree_set[treeIndex].tree[i][matching_tree_set[treeIndex].paths[q][i]]->queryToVertex[q];
            for (const auto& u_other: DAGs[global_hash[treeIndex][q]].GetTreeNodes()[u].neighbors_) {
                if (queryV[u_other] != UNMATCHED) {
                    continue;
                }
                uint qu = vertexMapping[std::make_pair(global_hash[treeIndex][q], u)];
                if (Num[m[i]][qu][u_other] < dynamicOrders[u_other].E) {
                    dynamicOrders[u_other].E = Num[m[i]][qu][u_other];
                    dynamicOrders[u_other].u_min = u;
                    dynamicOrders[u_other].index_min = i;
                }
                dynamicOrders[u_other].matched_nbrs ++;
            }
        }

        DynamicSearchMatches(depth + 1, treeIndex, q, m, queryV, dynamicOrders);
    }  
}

void Cimos::DynamicSearchMatches(uint depth, uint treeIndex, uint q, std::vector<uint>& m, std::vector<uint>& queryV, std::vector<ExtendableVertex>& dynamicOrders) {
    
    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;
    uint u, index, edgeLabel, targetLabel;

    // get fully matched result for specific query, and return
    if (depth == querys_[global_hash[treeIndex][q]].NumVertices()) {
        queryResults[global_hash[treeIndex][q]] ++;
        return;
    }

    for (uint i = 0; i < dynamicOrders.size(); i++) {
        if (queryV[i] != UNMATCHED) {
            continue;
        }
        // pick isolated vertex first
        if (dynamicOrders[i].matched_nbrs == querys_[global_hash[treeIndex][q]].GetNeighbors(i).size()) {
            if (dynamicOrders[i].E < isolate_minE) {
                isolate_minE = dynamicOrders[i].E;
                isolate_u = i;
            }
        } else { // pick normal vertex then
            if (dynamicOrders[i].E < non_isolate_minE) {
                non_isolate_minE = dynamicOrders[i].E;
                non_isolate_u = i;
            }
        }
    }

    if (non_isolate_minE == NOT_EXIST) {
        u = isolate_u;
    } else {
        u = non_isolate_u;
    }

    uint previousU = dynamicOrders[u].u_min;
    targetLabel = querys_[global_hash[treeIndex][q]].GetVertexLabel(u);

    for (const auto& v: EdgeIndex[m[queryV[previousU]]][targetLabel]) {

        uint qv = vertexMapping[std::make_pair(global_hash[treeIndex][q], u)];
        if (visited_[v] || !B[v][qv]) {
            continue;
        }

        bool isJoinable = true;
        for (auto& u_other: DAGs[global_hash[treeIndex][q]].GetTreeNodes()[u].neighbors_) {
            if (queryV[u_other] == UNMATCHED || u_other == previousU) {
                continue;
            } else {
                if (!EdgeIndex[m[queryV[u_other]]][targetLabel].count(v)) {
                    isJoinable = false;
                    break;
                }
            }
        }
        if (!isJoinable) {
            continue;
        }  
        
        m[depth] = v;
        visited_[v] = true;
        queryV[u] = depth;

        std::vector<ExtendableVertex> newDynamicOrders(dynamicOrders);
        for (const auto& u_other: DAGs[global_hash[treeIndex][q]].GetTreeNodes()[u].neighbors_) {
            if (queryV[u_other] != UNMATCHED) {
                continue;
            }
            uint qu = vertexMapping[std::make_pair(global_hash[treeIndex][q], u)];
            if (Num[v][qu][u_other] < newDynamicOrders[u_other].E) {
                newDynamicOrders[u_other].E = Num[v][qu][u_other];
                newDynamicOrders[u_other].u_min = u;   
            }
            newDynamicOrders[u_other].matched_nbrs ++;
        }

        DynamicSearchMatches(depth + 1, treeIndex, q, m, queryV, newDynamicOrders);

        queryV[u] = UNMATCHED;
        visited_[v] = false;
        m[depth] = UNMATCHED;
    }
}

uint Cimos::GetMaxTreeSize() {
    uint ans = 0;
    for (const auto& Tree: matching_tree_set) {
        ans = (Tree.tree.size() > ans) ? Tree.tree.size() : ans;
    }
    return ans;
}

void Cimos::AddEdge(uint v1, uint v2, uint label)
{
    std::chrono::high_resolution_clock::time_point s = Get_Time();
    data_.AddEdge(v1, v2, label);

    // enumerate all query edges that matches v1 --> v2
    for (uint qnum = 0; qnum < querys_.size(); qnum++) {      
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint u1 = 0; u1 < query.NumVertices(); u1++) {
            if (data_.GetVertexLabel(v1) == query.GetVertexLabel(u1)) {
                for (uint u2 = 0; u2 < query.NumVertices(); u2++) {
                    if (data_.GetVertexLabel(v2) == query.GetVertexLabel(u2)) {
                        // if (std::get<2>(query.GetEdgeLabel(u1, u2)) != label) continue;

                        bool reversed = false;
                        if (std::find(qDAG.GetTreeNodes()[u1].backwards_.begin(), qDAG.GetTreeNodes()[u1].backwards_.end(), u2) != qDAG.GetTreeNodes()[u1].backwards_.end()) {
                            std::swap(u1, u2);
                            std::swap(v1, v2);
                            reversed = true;
                        }
                        if (std::find(qDAG.GetTreeNodes()[u2].backwards_.begin(), qDAG.GetTreeNodes()[u2].backwards_.end(), u1) != qDAG.GetTreeNodes()[u2].backwards_.end()) {
                            if (!EdgeIndex[v1][data_.GetVertexLabel(v2)].count(v2)) {
                                EdgeIndex[v1][data_.GetVertexLabel(v2)].insert(v2);
                            }
                            if (!EdgeIndex[v2][data_.GetVertexLabel(v1)].count(v1)) {
                                EdgeIndex[v2][data_.GetVertexLabel(v1)].insert(v1);
                            }                            

                            bool old_p_F = F[v1][vertexMapping[std::make_pair(qnum, u1)]];
                            bool old_p_B = B[v1][vertexMapping[std::make_pair(qnum, u1)]];
                            bool old_c_B = B[v2][vertexMapping[std::make_pair(qnum, u2)]];

                            if (old_p_F) {
                                InsertionModifyF(qnum, u1, u2, v1, v2);
                            }
                            if (old_c_B) {
                                InsertionModifyB(qnum, u2, u1, v2, v1);
                            }
                            if (old_p_B) {
                                Num[v2][vertexMapping[std::make_pair(qnum, u2)]][u1] += 1;
                            }

                            while (!Q1.empty()) {
                                auto [v, u] = Q1.front();
                                Q1.pop();
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].forwards_.size(); k++) {
                                    const uint u_c = qDAG.GetTreeNodes()[u].forwards_[k];
                                    const uint u_c_label = query.GetVertexLabel(u_c);
                                    for (const auto& v_c: EdgeIndex[v][u_c_label]) {
                                        InsertionModifyF(qnum, u, u_c, v, v_c);
                                    }
                                }
                            }

                            while (!Q2.empty()) {
                                auto [v, u] = Q2.front();
                                Q2.pop();
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].backwards_.size(); k++) {
                                    const uint u_p = qDAG.GetTreeNodes()[u].backwards_[k];
                                    const uint u_p_label = query.GetVertexLabel(u_p);
                                    for (const auto& v_p: EdgeIndex[v][u_p_label]) {
                                        InsertionModifyB(qnum, u, u_p, v, v_p);
                                    }
                                }
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].forwards_.size(); k++) {
                                    const uint u_c = qDAG.GetTreeNodes()[u].forwards_[k];
                                    const uint u_c_label = query.GetVertexLabel(u_c);
                                    for (const auto& v_c: EdgeIndex[v][u_c_label]) {
                                        Num[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] += 1;
                                    }
                                }
                            }
                        }

                        if (reversed)
                        {
                            std::swap(u1, u2);
                            std::swap(v1, v2);
                        }
                    }
                }
            }
        }
    }

    if (!(D[v1] && D[v2])) {
        return;
    }

    uint maxVertexNumber = GetMaxTreeSize();
    std::vector<uint> m(maxVertexNumber, UNMATCHED);

    // enumerate all query edges that matches v1 --> v2
    // size_t num_results = 0ul;
    uint label1 = data_.GetVertexLabel(v1);
    uint label2 = data_.GetVertexLabel(v2);
    if (label1 > label2) {
        std::swap(label1, label2);
        std::swap(v1, v2);
    }
    std::vector<ExtendableVertex> dynamicOrders = {};

    for (uint i = 0; i < all_vector_unique_edges.size(); i++) {
        if (label1 == std::get<0>(all_vector_unique_edges[i]) && label2 == std::get<1>(all_vector_unique_edges[i]) && label == std::get<2>(all_vector_unique_edges[i])) {
            if (matching_tree_set[i].is_reverse) {
                m[0] = v2;
                m[1] = v1;
            } else {
                m[0] = v1;
                m[1] = v2;
            }
            visited_[v1] = true;
            visited_[v2] = true;

            StaticSearchMatches(1, i, 0, m);

            m[0] = UNMATCHED;
            m[1] = UNMATCHED;
            visited_[v1] = false;
            visited_[v2] = false;

            if (label1 == label2) {
                if (matching_tree_set[i].is_reverse) {
                    m[0] = v1;
                    m[1] = v2;
                } else {
                    m[0] = v2;
                    m[1] = v1;
                }
                visited_[v1] = true;
                visited_[v2] = true;

                StaticSearchMatches(1, i, 0, m);

                m[0] = UNMATCHED;
                m[1] = UNMATCHED;
                visited_[v1] = false;
                visited_[v2] = false;
            }
        }

    }
}

void Cimos::RemoveEdge(uint v1, uint v2, uint label)
{
    uint maxVertexNumber = GetMaxTreeSize();
    std::vector<uint> m(maxVertexNumber, UNMATCHED);

    uint label1 = data_.GetVertexLabel(v1);
    uint label2 = data_.GetVertexLabel(v2);
    if (label1 > label2) {
        std::swap(label1, label2);
        std::swap(v1, v2);
    }
    std::vector<ExtendableVertex> dynamicOrders = {};

    for (uint i = 0; i < all_vector_unique_edges.size(); i++) {
        if (label1 == std::get<0>(all_vector_unique_edges[i]) && label2 == std::get<1>(all_vector_unique_edges[i]) && label == std::get<2>(all_vector_unique_edges[i])) {
            if (matching_tree_set[i].is_reverse) {
                m[0] = v2;
                m[1] = v1;
            } else {
                m[0] = v1;
                m[1] = v2;
            }
            visited_[v1] = true;
            visited_[v2] = true;

            StaticSearchMatches(1, i, 0, m);

            m[0] = UNMATCHED;
            m[1] = UNMATCHED;
            visited_[v1] = false;
            visited_[v2] = false;

            if (label1 == label2) {
                if (matching_tree_set[i].is_reverse) {
                    m[0] = v1;
                    m[1] = v2;
                } else {
                    m[0] = v2;
                    m[1] = v1;
                }
                visited_[v1] = true;
                visited_[v2] = true;

                StaticSearchMatches(1, i, 0, m);

                m[0] = UNMATCHED;
                m[1] = UNMATCHED;
                visited_[v1] = false;
                visited_[v2] = false;
            }
        }

    }

    for (uint qnum = 0; qnum < querys_.size(); qnum++) {
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint u1 = 0; u1 < query.NumVertices(); u1++) {
            if (data_.GetVertexLabel(v1) == query.GetVertexLabel(u1)) {
                for (uint u2 = 0; u2 < query.NumVertices(); u2++) {
                    if (data_.GetVertexLabel(v2) == query.GetVertexLabel(u2)) {

                        auto it = std::lower_bound(query.GetNeighbors(u1).begin(), query.GetNeighbors(u1).end(), u2);
                        if (it == query.GetNeighbors(u1).end() || *it != u2) {
                            continue;
                        }
                        
                        bool reversed = false;
                        if (std::find(qDAG.GetTreeNodes()[u1].backwards_.begin(), qDAG.GetTreeNodes()[u1].backwards_.end(), u2) != qDAG.GetTreeNodes()[u1].backwards_.end()) {
                            std::swap(u1, u2);
                            std::swap(v1, v2);
                            reversed = true;
                        }
                        if (std::find(qDAG.GetTreeNodes()[u2].backwards_.begin(), qDAG.GetTreeNodes()[u2].backwards_.end(), u1) != qDAG.GetTreeNodes()[u2].backwards_.end()) {

                            EdgeIndex[v1][data_.GetVertexLabel(v2)].erase(v2);
                            EdgeIndex[v2][data_.GetVertexLabel(v1)].erase(v1);

                            bool old_p_F = F[v1][vertexMapping[std::make_pair(qnum, u1)]];
                            bool old_p_B = B[v1][vertexMapping[std::make_pair(qnum, u1)]];
                            bool old_c_B = B[v2][vertexMapping[std::make_pair(qnum, u2)]];

                            if (old_c_B) {
                                DeletionModifyB(qnum, u2, u1, v2, v1);
                            }
                            if (old_p_B) {
                                Num[v2][vertexMapping[std::make_pair(qnum, u2)]][u1] -= 1;
                            }
                            if (old_p_F) {
                                DeletionModifyF(qnum, u1, u2, v1, v2);
                            }

                            while (!Q1.empty())
                            {
                                auto [v, u] = Q1.front();
                                Q1.pop();
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].forwards_.size(); k++) {
                                    const uint u_c = qDAG.GetTreeNodes()[u].forwards_[k];
                                    const uint u_c_label = query.GetVertexLabel(u_c);
                                    for (const auto& v_c: EdgeIndex[v][u_c_label]) {
                                        DeletionModifyF(qnum, u, u_c, v, v_c);
                                    }
                                }
                            }
                            while (!Q2.empty())
                            {
                                auto [v, u] = Q2.front();
                                Q2.pop();
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].backwards_.size(); k++) {
                                    const uint u_p = qDAG.GetTreeNodes()[u].backwards_[k];
                                    const uint u_p_label = query.GetVertexLabel(u_p);
                                    for (const auto& v_p: EdgeIndex[v][u_p_label]) {
                                        DeletionModifyB(qnum, u, u_p, v, v_p);
                                    }
                                }
                                for (uint k = 0; k < qDAG.GetTreeNodes()[u].forwards_.size(); k++) {
                                    const uint u_c = qDAG.GetTreeNodes()[u].forwards_[k];
                                    const uint u_c_label = query.GetVertexLabel(u_c);
                                    for (const auto& v_c: EdgeIndex[v][u_c_label]) {
                                        Num[v_c][vertexMapping[std::make_pair(qnum, u_c)]][u] -= 1;
                                    }
                                }
                            }
                        }
                        if (reversed)
                        {
                            std::swap(u1, u2);
                            std::swap(v1, v2);
                        }

                    }
                }
            }
        }
    }
    data_.RemoveEdge(v1, v2);
}

void Cimos::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
    
    visited_.resize(id + 1, false);

    for (uint qnum = 0; qnum < querys_.size(); qnum++) {
        QueryDAG& qDAG = DAGs[qnum];
        Graph& query = querys_[qnum];
        for (uint i = 0; i < query.NumVertices(); i++) {
            const uint u = qDAG.GetSerializedTree()[i];
            if (data_.GetVertexLabel(id) == query.GetVertexLabel(u)) {  
                if (qDAG.GetTreeNodes()[u].backwards_.empty()) {
                    F[id][vertexMapping[std::make_pair(qnum, u)]] = 1;
                }  
                if (F[id][vertexMapping[std::make_pair(qnum, u)]] && qDAG.GetTreeNodes()[u].forwards_.empty()) {
                    B[id][vertexMapping[std::make_pair(qnum, u)]] = 1;
                }  
            }
        }
    }
}

void Cimos::RemoveVertex(uint id)
{
    // data_.RemoveVertex(id);
}

void Cimos::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{

}

void Cimos::LoadAllEdges(std::vector<std::string>& paths) {
    for (uint i = 0; i < paths.size(); i++) {
        std::ifstream ifs(paths[i]);
        char type;
        while (ifs >> type) {
            if (type == 't') {
                uint temp1, temp2;
                ifs >> temp1 >> temp2;
            }
            else if (type == 'v') {
                uint vid, vlabel, unk;
                ifs >> vid >> vlabel;
                continue;
            } else {
                uint from_id, to_id, edge_label = 0;
                ifs >> from_id >> to_id ;
                edge_label = 0;
                uint label1 = querys_[i].GetVertexLabel(from_id);
                uint label2 = querys_[i].GetVertexLabel(to_id);
                if (label1 > label2) {
                    std::swap(label1, label2);
                    std::swap(from_id, to_id);
                }
                all_unique_edges.emplace(label1, label2, edge_label);
            }
        }
        ifs.close();
    }

    all_vector_unique_edges.assign(all_unique_edges.begin(), all_unique_edges.end());
        // use a vector to store instead of set
    edge_in_graphs.resize(all_vector_unique_edges.size());  
}

void Cimos::CalculateEdgeMapping() {

    uint counter = 0;
    for (uint i = 0; i < querys_.size(); i++) {
        for (uint j = 0; j < querys_[i].NumVertices(); j++) {
            vertexMapping[std::make_pair(i, j)] = counter ++;
        }
    }

    for (uint i = 0; i < all_vector_unique_edges.size(); i++) {

        uint label1 = std::get<0>(all_vector_unique_edges[i]);
        uint label2 = std::get<1>(all_vector_unique_edges[i]);
        uint edge_label = std::get<2>(all_vector_unique_edges[i]);

        for (uint j = 0; j < querys_.size(); j++) {
            for (const auto &k: querys_[j].edges) {
                // find matched edge instances
                if ((label1 == std::get<0>(k)) && (label2 == std::get<1>(k)) && (edge_label == std::get<2>(k))) {
                    // a match occurs, from query graph [j]
                    if (edge_in_graphs[i].find(j) != edge_in_graphs[i].end()) {
                        // query graph j already exist, just emplace back
                        edge_in_graphs[i][j].emplace_back(std::get<3>(k), std::get<4>(k));
                    } else {
                        std::vector<std::pair<uint, uint>> tmp = {{std::get<3>(k), std::get<4>(k)}};
                        edge_in_graphs[i][j] = tmp;
                    }
                }
            }
        }
    }
}

void Cimos::ConstructMatchingTrees() {
    global_hash.resize(edge_in_graphs.size() + 1);
    uint sampleNumber = 10000;

    for (uint i = 0; i < edge_in_graphs.size(); i++) {

        std::vector<std::pair<uint, uint>> edge_instances;
        std::vector<uint> q_ids;
        uint query_number = 0;
        for (const auto& pair: edge_in_graphs[i]) {
            uint query_id = pair.first;
            for (const auto& edges: pair.second) {
                edge_instances.push_back(edges);
                global_hash[i][query_number] = query_id;
                q_ids.push_back(query_number);
                query_number++;
            }
        }
        matching_tree_set.push_back(ConstructOneTree(edge_instances, q_ids, i, sampleNumber));
    }
}

void Cimos::CalculateOverallStatistic() {

    std::unordered_map<std::pair<uint, uint>, uint, PairHash> globalall;
    for (uint i = 0; i < data_.NumVertices(); i++) {
        for (uint j = 0; j < querys_.size(); j++) {
            for (uint v = 0; v < querys_[j].NumVertices(); v++) {
                for (const auto& v_nei: querys_[j].GetNeighbors(v)) {
                    std::pair<uint, uint> pr = std::make_pair(vertexMapping[std::make_pair(j, v)], vertexMapping[std::make_pair(j, v_nei)]);
                    globalall[pr] += Num[i][vertexMapping[std::make_pair(j, v)]][v_nei];
                }
            }
        }
    }

    for (const auto& ele: globalall) {
        globalN[ele.first] = (1.0 * ele.second) / (1.0 * data_.NumVertices());
    }
}

void Cimos::SampleEdges(uint label1, uint label2, uint sampleNumber) {

    for (uint u = 0; u < data_.neighbors_.size(); ++u) {
        for (uint v : data_.neighbors_[u]) {
            if (u < v) {
                if (( data_.vlabels_[u] == label1 &&  data_.vlabels_[v] == label2) || 
                    ( data_.vlabels_[u] == label2 &&  data_.vlabels_[v] == label1)) {
                    std::vector<uint> tp;
                    tp.push_back(u);
                    tp.push_back(v);
                    intermediateResults[std::make_pair(1, 0)].push_back(tp);
                }
            }
        }
    }

    if (intermediateResults[std::make_pair(1, 0)].size() > sampleNumber) {
        // randomly sample K edges
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(intermediateResults[std::make_pair(1, 0)].begin(), intermediateResults[std::make_pair(1, 0)].end(), g);
        intermediateResults[std::make_pair(1, 0)].resize(sampleNumber); 
    }



}

void Cimos::SamplingSubgraphs(uint i, uint j, uint treeIndex, MatchingTree& Mtree, std::vector<uint>& path) {

    std::pair<uint, uint> startPair;
    for (uint depth = path.size() - 1; depth >= 0; depth--) {
        std::pair<uint, uint> pair = std::make_pair(depth, path[depth]);
        if (Mtree.isTreeNodeStoresIR[pair] == true) {
            startPair = pair;
            break;
        }
    }

    if (startPair.first == i && startPair.second == j) {
        // denote that this node has IRs already
        return;
    }

    for (auto IR: intermediateResults[startPair]) {

        std::queue<std::vector<uint>> newIR;
        newIR.push(IR);
        uint number = 1;

        for (uint s = startPair.first + 1; s <= i; s++) {

            number = newIR.size();

            for (uint k = 0; k < number; k++) {

                std::vector<uint> currentIR = newIR.front();
                Node* node = Mtree.tree[s][path[s]];
                uint targetLabel = node->label;
                uint minValue = NOT_EXIST;
                uint minIndex;
                uint exampleQ = node->queryToVertex.begin()->first;
                uint exampleU = node->queryToVertex.begin()->second;

                for (const auto& tp: node->topology) {
                    uint previousV = currentIR[tp.first];
                    uint previousPosition = path[tp.first]; 
                    uint previousU = Mtree.tree[tp.first][previousPosition]->queryToVertex[exampleQ];
                    uint qu = vertexMapping[std::make_pair(global_hash[treeIndex][exampleQ], previousU)];
                    uint candidatesSize = Num[previousV][qu][exampleU];

                    if (candidatesSize < minValue) {
                        minValue = candidatesSize;
                        minIndex = tp.first;
                    }
                }

                if (minValue == 0) {
                    newIR.pop();
                    continue;
                }

                uint minV = currentIR[minIndex];
                for (const auto& v: EdgeIndex[minV][targetLabel]) {

                    if (!D[v]) {
                        continue;
                    }

                    bool flag = false;
                    for (const auto& ele: currentIR) {
                        if (ele == v) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        continue;
                    }

                    bool isJoinable = true;
                    for (const auto& tp: node->topology) {
                        if (tp.first == minIndex) {
                            continue;
                        } else {
                            if (!EdgeIndex[currentIR[tp.first]][targetLabel].count(v)) {
                                isJoinable = false;
                                break;
                            }
                        }
                    }
                    if (!isJoinable) {
                        continue;
                    }
                    std::vector<uint> newOne = currentIR;
                    newOne.push_back(v);
                    newIR.push(newOne);

                }

                newIR.pop();

            }
        }

        while (!newIR.empty()) {
            intermediateResults[std::make_pair(i, j)].push_back(newIR.front());
            newIR.pop(); 
        }
        Mtree.isTreeNodeStoresIR[std::make_pair(i, j)] = true;
    }
}

uint Cimos::CalculateMinimalTrie(std::vector<uint>& maxLabels, std::vector<std::vector<std::pair<uint, uint>>>& maxVector, uint i, uint j) {

}


MatchingTree Cimos::ConstructOneTree(std::vector<std::pair<uint, uint>> &maps, std::vector<uint> &ids, uint index, uint sampleNumber) {

    MatchingTree matchingTree;

    std::unordered_map<uint, std::bitset<15>> visit;
    // we assume that the vertex number of each query graph will not exceed 15.
    for (const auto it : ids) {
        for (uint i = 0; i < 15; i++) {
            visit[it][i] = false;
        }
    }

    uint degree1 = 0, degree2 = 0;
    std::vector<uint> list1, list2;
    for (uint i = 0; i < maps.size(); i++) {
        degree1 += querys_[global_hash[index][ids[i]]].GetDegree(std::get<0>(maps[i]));
        degree2 += querys_[global_hash[index][ids[i]]].GetDegree(std::get<1>(maps[i]));
        list1.push_back(std::get<0>(maps[i]));
        list2.push_back(std::get<1>(maps[i]));
        visit[ids[i]][std::get<0>(maps[i])] = true;  // when meet with visit array, we must use ids[i] instead of gloabl hash[i]
        visit[ids[i]][std::get<1>(maps[i])] = true;
    }
    float d1 = degree1 / maps.size();
    float d2 = degree2 / maps.size();

    // construct root node and its only children node
    Node node1, node2;

    node1.label = std::get<0>(all_vector_unique_edges[index]);
    for (uint i = 0; i < ids.size(); i++) {
        node1.queryToVertex[ids[i]] = list1[i]; 
    }

    node2.label = std::get<1>(all_vector_unique_edges[index]);
    for (uint i = 0; i < ids.size(); i++) {
        node2.queryToVertex[ids[i]] = list2[i]; 
    }    

    uint label1, label2;
    if (d1 > d2) {
        matchingTree.tree.push_back({new Node(node2)});
        matchingTree.tree.push_back({new Node(node1)});
        matchingTree.is_reverse = true;
        label1 = node2.label;
        label2 = node1.label;
    } else {
        matchingTree.tree.push_back({new Node(node1)});
        matchingTree.tree.push_back({new Node(node2)});
        matchingTree.is_reverse = false;
        label1 = node2.label;
        label2 = node1.label;
    }
    matchingTree.children.push_back({{0}});
    for (const auto& i : ids) {
        matchingTree.paths[i] = {};
        matchingTree.paths[i].push_back(0);
        matchingTree.paths[i].push_back(0);
    }

    SampleEdges(label1, label2, sampleNumber);
    matchingTree.isTreeNodeStoresIR[std::make_pair(1, 0)] = true; 

    // get max vertex number in query graphs
    uint max_vertex = 0;
    for (const auto& i : ids) {
        max_vertex = (querys_[global_hash[index][i]].NumVertices() > max_vertex) ? querys_[global_hash[index][i]].NumVertices() : max_vertex;
    }

    for (uint i = 2; i <= max_vertex; i++) {
        uint parent_number = matchingTree.tree[i - 1].size();
        uint counter = 0;
        matchingTree.tree.emplace_back();
        // for each parent, we make the next vertex label to match, and decide whether split it into multiple children or not
        std::vector<std::vector<uint>> layer_children_vector = {};
        for (uint j = 0; j < parent_number; j++) {

            std::vector<uint> children_vector = {};
            Node node = *(matchingTree.tree[i - 1][j]);
            std::vector<Node*> layer_nodes;

            // first round to select which vertex to expand
            // choose random one, such as 0, because they are all the same for path
            std::vector<uint> path = matchingTree.paths[node.queryToVertex.begin()->first];
            uint last_label = 99999;
            std::vector<uint> classified_queries_remain;

            // delete queries that all verices have been visited
            std::vector<uint> delete_query;
            
            for (const auto& [q, v]: node.queryToVertex) {
                if (node.truncatedQueries.count(q) > 0) {
                    continue;
                }
                uint cc = 0;
                for (uint v = 0; v < querys_[global_hash[index][q]].NumVertices(); v++) {
                    if (visit[q][v]) cc++;
                }
                if (cc == querys_[global_hash[index][q]].NumVertices()) {
                    delete_query.push_back(q);
                } else {
                    classified_queries_remain.push_back(q);
                }           
            }

            for (const auto& dq: delete_query) {
                matchingTree.tree[i - 1][j]->endQueries.insert(dq);
            }

            // continue if there is no left queries
            if (classified_queries_remain.empty())  {
                children_vector.push_back(-1);
                layer_children_vector.push_back(children_vector);
                continue;
            }

            // below are the strategy to select next matching vertex(label),
            std::unordered_set<uint> allLabels;

            // 统计所有的邻居label
            for (const auto& query_id : classified_queries_remain) {
                for (uint u = 0; u < path.size(); u++) {
                    uint this_vertex = matchingTree.tree[u][path[u]]->queryToVertex[query_id];
                    for (const auto& ver: querys_[global_hash[index][query_id]].GetNeighbors(this_vertex)) {
                        if (visit[query_id][ver]) continue;
                        uint _label = querys_[global_hash[index][query_id]].GetVertexLabel(ver);
                        allLabels.insert(_label);
                    }
                }
            }


            // (label, (queryid, (vertexid, (hash_vector))))
            std::unordered_map<uint, std::unordered_map<uint, std::unordered_map<uint, std::vector<std::pair<uint, uint>>>>> hashElements;

            for (const auto& current_label : allLabels) {
                std::unordered_map<uint, std::unordered_map<uint, std::vector<std::pair<uint, uint>>>> hash_value;
                std::unordered_map<uint, std::vector<uint>> hash_vector;

                for (const auto& query_id: classified_queries_remain) {
                    uint note = 0;

                    for (uint u = 0; u < path.size(); u++) {
                        uint this_vertex = matchingTree.tree[u][path[u]]->queryToVertex[query_id];
                        for (const auto& ver: querys_[global_hash[index][query_id]].GetNeighbors(this_vertex)) {
                            if (visit[query_id][ver]) continue;
                            if (querys_[global_hash[index][query_id]].GetVertexLabel(ver) == current_label) {
                                hash_value[query_id][ver].resize(15);
                                note = 1;
                            }
                        }
                    }
                }  
                
                // compute hash vector
                for (auto& item: hash_value) {
                    uint query_id = item.first;
                    for (auto& pair: item.second) {
                        uint vertex_id = pair.first;

                        for (uint u = 0; u < path.size(); u++) {
                            uint this_vertex = matchingTree.tree[u][path[u]]->queryToVertex[query_id];
                            bool isconnect = false;
                            for (const auto& ver: querys_[global_hash[index][query_id]].GetNeighbors(this_vertex)) {
                                if (visit[query_id][ver]) continue;
                                if (ver == vertex_id) {
                                    uint edgeLabel = querys_[global_hash[index][query_id]].vertexPairToEdgeLabel[std::make_pair(this_vertex, ver)];
                                    hash_value[query_id][vertex_id][u] = {1, edgeLabel};
                                    isconnect = true;
                                    break;
                                }
                            }
                            if (!isconnect) {
                                // 0 represents no edge connected
                                hash_value[query_id][vertex_id][u] = {0, 0}; 
                            }
                        }
                    }
                }

                // Drop None Edges, fullfill hashElement data structure
                for (auto& item: hash_value) {
                    uint query_id = item.first;
                    for (auto& pair: item.second) {
                        uint vertex_id = pair.first;
                        std::vector<std::pair<uint, uint>> resultVector;
                        for (uint index = 0; index < pair.second.size(); index++) {
                            if (pair.second[index].first == 1) {
                                resultVector.emplace_back(index, pair.second[index].second);
                            }
                        }
                        hashElements[current_label][query_id][vertex_id] = resultVector;
                    }
                }
            }

            // (label, (hashVector, (queryid, {vertexid})))
            std::unordered_map<uint, std::unordered_map<std::vector<std::pair<uint, uint>>, std::unordered_map<uint, std::vector<uint>>, VectorHash>> sortedHashValues;

            for (const auto& [current_id, hashs]: hashElements) {
                for (const auto& [query_id, vectors]: hashs) {
                    for (const auto& [vertex_id, hashVector]: vectors) {
                        sortedHashValues[current_id][hashVector][query_id].push_back(vertex_id);
                    }
                }
            }

            std::unordered_map<uint, bool> isSelected;
            std::unordered_map<uint, std::unordered_map<std::vector<std::pair<uint, uint>>, bool, VectorHash>> isVectorSelected;
            for (const auto& queryid: classified_queries_remain) {
                isSelected[queryid] = false;
            }

            for (const auto& [current_label, hashs]: sortedHashValues) {
                for (const auto& [hashVector, queries]: hashs) {
                    isVectorSelected[current_label][hashVector] = false;
                }
            }

            bool isTruncate = false;


            while (true) {
                uint maxSize = 0;
                std::vector<uint> maxLabels;
                std::vector<std::vector<std::pair<uint, uint>>> maxVectors;

                uint maxLabel;
                std::vector<std::pair<uint, uint>> maxVector;
                for (const auto& [current_label, hashs]: sortedHashValues) {
                    for (const auto& [hashVector, queries]: hashs) {
                        if (isVectorSelected[current_label][hashVector]) {
                            // if this label's hashVector has been selected，iterate to the next directly
                            continue;
                        }
                        bool flag = false;
                        for (const auto& [query, vertexVector]: queries) {
                            if (isSelected[query]) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            continue;
                        }
                        uint querySize = queries.size();
                        if (querySize > maxSize) {
                            maxSize = querySize;
                        }
                    }
                }

                uint mN = 0;

                for (const auto& [current_label, hashs]: sortedHashValues) {
                    for (const auto& [hashVector, queries]: hashs) {
                        if (isVectorSelected[current_label][hashVector]) {
                            continue;
                        }
                        bool flag = false;
                        for (const auto& [query, vertexVector]: queries) {
                            if (isSelected[query]) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            continue;
                        }
                        uint querySize = queries.size();
                        if (querySize == maxSize) {
                            mN ++;
                            maxLabels.push_back(current_label);
                            maxVectors.push_back(hashVector);
                        }
                    }
                }
                if (maxSize > 1 && mN != 1) {
                    // trie occurs, and it needs sampling approach to decide which vertex to choose
                    SamplingSubgraphs(i - 1, j, index, matchingTree, path);
                    float Min = UNMATCHED;
                    uint MinK;
                    for (uint k = 0; k < maxLabels.size(); k++) {
                        std::vector<std::pair<uint, uint>> topoVec = maxVectors[k];
                        uint exampleQ = sortedHashValues[maxLabels[k]][topoVec].begin()->first;
                        uint exampleU = sortedHashValues[maxLabels[k]][topoVec].begin()->second[0];
                        uint tot = 0;

                        for (const auto& IR: intermediateResults[std::make_pair(i - 1, j)]) {
                            uint minValue = UNMATCHED;
                            for (const auto& tp: topoVec) {
                                uint previousV = IR[tp.first];
                                uint previousPosition = path[tp.first]; 
                                uint previousU = matchingTree.tree[tp.first][previousPosition]->queryToVertex[exampleQ];
                                uint qu = vertexMapping[std::make_pair(global_hash[index][exampleQ], previousU)];
                                uint candidatesSize = Num[previousV][qu][exampleU];
                                if (candidatesSize < minValue) {
                                    minValue = candidatesSize;
                                }
                            }
                            tot += minValue;
                        }

                        float finalResult;
                        if (intermediateResults[std::make_pair(i - 1, j)].size() == 0) {
                            finalResult = 0;
                        } else {
                            finalResult = (1.0 * tot) / (1.0 * intermediateResults[std::make_pair(i - 1, j)].size());
                        } 

                        if (finalResult < Min) {
                            Min = finalResult;
                            MinK = k;
                        }
                    }

                    maxLabel = maxLabels[MinK];
                    maxVector = maxVectors[MinK];

                    goto Merge;

                } else if (maxSize == 0) {
                    break;
                } else if (maxSize == 1) {
                    std::vector<uint> trunc;
                    for (const auto& [q, isVisit]: isSelected) {
                        if (!isVisit) {
                            trunc.push_back(q);
                        }
                    }

                    for (const auto& q: trunc) {
                        matchingTree.tree[i - 1][j]->truncatedQueries[q].resize(querys_[global_hash[index][q]].NumVertices(), UNMATCHED);
                        for (uint u = 0; u < path.size(); u++) {
                            uint pathV = matchingTree.tree[u][path[u]]->queryToVertex[q];
                            matchingTree.tree[i - 1][j]->truncatedQueries[q][pathV] = u;
                        }
                        for (uint u = 0; u < querys_[global_hash[index][q]].NumVertices(); u++) {
                            visit[q][u] = true;
                        }
                    }
                    break;
                } else {
                    
                    maxLabel = maxLabels[0];
                    maxVector = maxVectors[0];
                Merge: 
                    std::vector<uint> queries;
                    std::vector<uint> vertices;
                    for (const auto& [query, vertexVector]: sortedHashValues[maxLabel][maxVector]) {
                        queries.push_back(query);
                        vertices.push_back(vertexVector[0]);
                        isSelected[query] = true;
                    }
                    std::unordered_set<uint> queriesSet;
                    for (const auto& q: queries) {
                        queriesSet.insert(q);
                    }
                    isVectorSelected[maxLabel][maxVector] = true;
                    Node new_node;
                    new_node.label = maxLabel;
                    for (uint k = 0; k < queries.size(); k++) {
                        new_node.queryToVertex[queries[k]] = vertices[k];
                        visit[queries[k]][vertices[k]] = true;
                    }
                    new_node.topology = maxVector;
                    layer_nodes.push_back(new Node(new_node));


                    for (auto& [current_label, hashs]: sortedHashValues) {
                        for (auto& [hashVector, queries]: hashs) {
                            for (auto it = queries.begin(); it != queries.end(); ) {
                                if (isSelected[it->first]) {
                                    it = queries.erase(it); 
                                } else {
                                    ++it;
                                }
                            }
                        }
                    }                    
                }



            }

            // maintain tree information
            if (layer_nodes.empty()) {
                children_vector.push_back(-1);
            }
            for (const auto& item: layer_nodes) {
                for (const auto& [each_query, v] : (*item).queryToVertex) {
                    matchingTree.paths[each_query].push_back(counter);
                }
                matchingTree.tree[i].push_back(item);
                children_vector.push_back(counter);
                counter++;
            }
            layer_children_vector.push_back(children_vector);

        }
        matchingTree.children.push_back(layer_children_vector);
    }

    // clear all subgraph matching IRs for each tree construction
    intermediateResults.clear();
    return matchingTree;
}

void Cimos::PrintQueryPlan() {
    for (uint i = 0; i < matching_tree_set.size(); i++) {
        std::cout <<  "EDGE: " << i << " NEXT::::::::::::::::::::::::::::::::::::::::::::::\n";
        for (auto& iii : matching_tree_set[i].tree) {
            std::cout << "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLAYER: " << "\n";
            for (auto & jjj : iii) {
                std::cout << "LABEL: " << jjj->label << "\nSURVIVE: ";
                for (const auto& [q, v] : jjj->queryToVertex) {
                    std::cout << q << "(" << global_hash[i][q] << ")" << " -> " << v << " ";
                }
                std::cout << "\n";
                for (const auto& pp : jjj->topology) {
                    std::cout << "(" << std::get<0>(pp) << ", " << std::get<1>(pp) << ")" << " ";
                    std::cout << "-----------\n";
                }
                if (jjj->truncatedQueries.size() != 0) {
                    std::cout << "TRUNCATED: ";
                    for (const auto& qqq: jjj->truncatedQueries) {
                        std::cout << qqq.first << " ";
                    }
                    std::cout << "\n";
                }
                 if (jjj->endQueries.size() != 0) {
                    std::cout << "ENDED: ";
                    for (const auto& qqq: jjj->endQueries) {
                        std::cout << qqq << " ";
                    }
                    std::cout << "\n";
                }               
                std::cout << "\n";
            }
        }
    }

}

void Cimos::PrintQueryReuslts() {
    for (uint i = 0; i < queryResults.size(); i++) {
        std::cout << "Query " << i << " , the number of matching is " << queryResults[i] << "\n";
    }
}
