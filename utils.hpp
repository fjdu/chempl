#ifndef UTILS_H
#define UTILS_H

#include <string>

namespace UTILS {

extern std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
 
extern std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
 
extern std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

}

#endif // UTILS_H
