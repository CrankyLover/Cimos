#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <type_traits>
#include <memory>

#include "utils/CLI11.hpp"
#include "utils/globals.h"
#include "utils/types.h"

#include "graph/graph.h"
#include "matching/matching.h"
#include "matching/cimos.h"

int main(int argc, char *argv[])
{

    CLI::App app{"App description"};

    std::string query_path = "", initial_path = "", stream_path = "";
    std::string query_folder_path = "";
    std::string data_path = "";
    std::string update_path = "";
    std::vector<std::string> query_paths;

    app.add_option("-q,--query", query_path, "query graph path")->required();
    app.add_option("-d,--data", initial_path, "initial data graph path")->required();
    app.add_option("-u,--update", stream_path, "data graph update stream path")->required();
    
    CLI11_PARSE(app, argc, argv);
    
    std::chrono::high_resolution_clock::time_point start, lstart;
    query_paths = GetQueryList(query_path);
    
    start = Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;

    std::vector<Graph> query_graphs;
    for (const auto& path: query_paths) {
        Graph query{};
        query.LoadFromFile(path, true);
        query.PrintMetaData();
        query_graphs.push_back(query);
    }

    Graph data_graph {};
    data_graph.LoadFromFile(initial_path);
    data_graph.PrintMetaData();
    Print_Time("Load Graphs: ", start);

    std::cout << "------------ Preprocessing ------------" << std::endl;
    Cimos *cimos = nullptr;

    start = Get_Time();

    cimos = new Cimos(query_graphs, data_graph);

    cimos->LoadAllEdges(query_paths);
    cimos->CalculateEdgeMapping();
    start = Get_Time();
    cimos->Preprocessing();
    Print_Time("Preprocessing: ", start);

    start = Get_Time();
    cimos->ConstructMatchingTrees();
    Print_Time("Generate Static Query plan: ", start);
    // cimos->PrintQueryPlan();
    
    std::cout << "--------- Incremental Matching --------" << std::endl;
    data_graph.LoadUpdateStream(stream_path);
    
    size_t num_v_updates = 0ul, num_e_updates = 0ul;

    std::vector<double> timeList;

    auto IncrementalFun = [&data_graph, &cimos, &num_v_updates, &num_e_updates]()
    {
        uint count = 0;
        std::vector<double> timeList;
        std::chrono::high_resolution_clock::time_point s = Get_Time();
        cimos->vs = Get_Time();
        while (!data_graph.updates_.empty())
        {
            InsertUnit insert = data_graph.updates_.front();
            data_graph.updates_.pop();
            if (insert.type == 'v' && insert.is_add)
            {
                cimos->AddVertex(insert.id1, insert.label);
                num_v_updates ++;
           }
            else if (insert.type == 'v' && !insert.is_add)
            {
                cimos->RemoveVertex(insert.id1);
                num_v_updates ++;
            }
            else if (insert.type == 'e' && insert.is_add)
            {
                cimos->AddEdge(insert.id1, insert.id2, insert.label);
                num_e_updates ++;
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                cimos->RemoveEdge(insert.id1, insert.id2, insert.label);
                num_e_updates ++;
            }
            if (reach_time_limit) break;
            count ++;
        }
    };

    start = Get_Time();
    uint time_limit = UINT_MAX;
    execute_with_time_limit(IncrementalFun, time_limit, reach_time_limit);
    Print_Time("Incremental Matching: ", start);
    cimos->PrintQueryReuslts();
    std::cout << "\n\n----------------- End -----------------" << std::endl;

    delete cimos;
}
