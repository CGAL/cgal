/**
 * @file   util/StringFuncs.h
 * @author Gernot Walzl
 * @date   2013-12-16
 */

#ifndef UTIL_STRINGFUNCS_H
#define UTIL_STRINGFUNCS_H

#include <string>
#include <vector>

namespace util {

class StringFuncs {
public:
    virtual ~StringFuncs();
    static bool startsWith(const std::string& str, const std::string& prefix);
    static bool endsWith(const std::string& str, const std::string& suffix);
    static std::string trim(const std::string& str);
    static std::vector<std::string> split(const std::string& str, const std::string& delimiter, bool keep_empty);
};

}

#endif /* UTIL_STRINGFUNCS_H */
