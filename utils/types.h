#ifndef UTILS_TYPES
#define UTILS_TYPES

#include <chrono>
#include <iostream>
#include <vector>
#include <climits>
#include <functional>
#include <stdlib.h>
#include <sys/stat.h> /* For stat() */
#include <dirent.h>

#define NOT_EXIST UINT_MAX
#define UNMATCHED UINT_MAX

// Time counting
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast<\
    std::chrono::microseconds>(Get_Time() - start).count()/(float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << \
    "ms" << std::endl


struct InsertUnit {
    char type;  // 'v' or 'e' 
    bool is_add;// addition or deletion
    uint id1;   // vertex id or edge source id
    uint id2;   // edge target id
    uint label; // vertex or edge label
    InsertUnit(char type_arg, bool is_add_arg, uint id1_arg, uint id2_arg, uint label_arg)
    : type(type_arg), is_add(is_add_arg), id1(id1_arg), id2(id2_arg), label(label_arg) {}
};

// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
// auxiliary generic functions to create a hash value using a seed
template <typename T> inline void hash_val(std::size_t &seed, const T &val) {
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return hash_val(p.first, p.second);
    }
};

inline bool isFile(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    return S_ISREG(info.st_mode);
}

inline std::vector<std::string> GetQueryList(std::string query_path) {
    std::vector<std::string> query_path_list;
    DIR* dir = opendir(query_path.c_str());
    if (dir == nullptr) {
        std::cerr << "Could not open directorty: " << query_path << "\n";
    }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string file_name = entry->d_name;
        if (file_name == "." || file_name == "..") {
            continue;
        }
        std::string file_path = query_path + "/" + file_name;
        if (isFile(file_path)) {
            query_path_list.push_back(file_path);
        }
    }
    closedir(dir);

    return query_path_list;
}


#endif //UTILS_TYPES
