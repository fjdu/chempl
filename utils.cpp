#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <algorithm>
#include <vector>

namespace UTILS {

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}

double interpol(const std::vector<double>& ts,
                const std::vector<double>& vs,  double t) {
  if (t <= ts[0]) {
    return vs[0];
  }
  if (t >= ts.back()) {
    return vs.back();
  }

  auto up = std::upper_bound(ts.begin(), ts.end(), t);
  int i = up - ts.begin() - 1;
  double k = (vs[i+1] - vs[i]) / (ts[i+1] - ts[i]);
  return vs[i] + k * (t - ts[i]);
}


}

#endif // UTILS_H
