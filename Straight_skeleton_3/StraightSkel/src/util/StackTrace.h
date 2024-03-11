/**
 * @file   util/StackTrace.h
 * @author Gernot Walzl
 * @date   2012-02-28
 */

#ifndef UTIL_STACKTRACE_H
#define UTIL_STACKTRACE_H

#ifdef __linux__
#include <execinfo.h>
#include <cxxabi.h>
#endif
#include <iostream>
#include <string>

namespace util {

class StackTrace {
public:
    virtual ~StackTrace();
    static std::string demangle(const std::string& symbol_mangled);
    static void print(std::ostream& os);
};

}

#endif /* UTIL_STACKTRACE_H */
