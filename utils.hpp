#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <algorithm>
#include <vector>

namespace UTILS {

extern std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
 
extern std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
 
extern std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

double interpol(const std::vector<double>& ts,
                const std::vector<double>& vs, double t);

}

#endif // UTILS_H
