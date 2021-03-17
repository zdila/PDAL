#if !defined(LIDARTRAJECTORY_PARAMS_HPP)
#define LIDARTRAJECTORY_PARAMS_HPP 1

#include <map>
#include <string>
#include <iostream>
#include <stdexcept>

namespace LidarTrajectory {
  class Params {
    std::map<std::string, std::string> _params;
  public:
    Params() {}
    void Clear() { _params.clear(); }
    void Add(std::istream& str);
    void Add(const std::string& string);
    void Init(std::istream& str) { Clear(); Add(str); }
    void Init(const std::string& file);
    Params(std::istream& str) { Add(str); }
    Params(std::string& file) { Init(file); }
    // Lookup a parameter with a default
    template <typename T>
    T lookup(const std::string& key, const T& defval) const;
    // Lookup a parameter
    template <typename T>
    T lookup(const std::string& key) const;
  };

}

#endif
