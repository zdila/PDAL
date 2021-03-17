#include <fstream>
#include <LidarTrajectory/Params.hpp>
#include <GeographicLib/Utility.hpp>

namespace LidarTrajectory {

  void Params::Add(std::istream& str) {
    std::string line;
    while (std::getline(str, line))
      Add(line);
  }

  void Params::Add(const std::string& string) {
    std::string key, val, temp = string;
    std::string::size_type n = temp.find('=');
    if (n != std::string::npos)
      temp[n] = ' ';        // This will be unnecessary with GeographicLib 1.51
    if (GeographicLib::Utility::ParseLine(temp, key, val))
        _params[key] = val;
  }

  void Params::Init(const std::string& file) {
    Clear();
    std::ifstream fstr(file.c_str());
    Add(fstr);
  }

  template <typename T>
  T Params::lookup(const std::string& key, const T& defval) const {
    auto r = _params.find(key);
    return r == _params.end() ? defval :
      // This will throw is value can't be read as a T
      GeographicLib::Utility::val<T>(r->second);
  }
  // Lookup a parameter
  template <typename T>
  T Params::lookup(const std::string& key) const {
    auto r = _params.find(key);
    if (r == _params.end())
      throw std::runtime_error("Unknown key " + key);
    // This will throw is value can't be read as a T
    return GeographicLib::Utility::val<T>(r->second);
  }

#define PARAMS_LOOKUP_INSTANTIATE(T) \
  template T Params::lookup<T>(const std::string&, const T&) const; \
  template T Params::lookup<T>(const std::string&) const;

  PARAMS_LOOKUP_INSTANTIATE(std::string)
  PARAMS_LOOKUP_INSTANTIATE(bool)
  PARAMS_LOOKUP_INSTANTIATE(int)
  PARAMS_LOOKUP_INSTANTIATE(double)

}
