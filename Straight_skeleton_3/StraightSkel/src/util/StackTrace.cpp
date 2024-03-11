/**
 * @file   util/StackTrace.cpp
 * @author Gernot Walzl
 * @date   2012-02-28
 */

#include "util/StackTrace.h"

namespace util {

StackTrace::~StackTrace() {
    // intentionally does nothing
}

std::string StackTrace::demangle(const std::string& symbol_mangled) {
    std::string result;
#ifdef __linux__
    size_t pos_begin = symbol_mangled.find('(');
    size_t pos_end = symbol_mangled.find('+');
    if (pos_begin != std::string::npos && pos_end != std::string::npos &&
            pos_begin < pos_end) {
        std::string symbol = symbol_mangled.substr(pos_begin+1, pos_end-pos_begin-1);
        int status;
        char * realname = abi::__cxa_demangle(symbol.c_str(), 0, 0, &status);
        if (realname) {
            result.append(symbol_mangled.substr(0, pos_begin+1));
            result.append(realname);
            result.append(symbol_mangled.substr(pos_end));
            free(realname);
        } else {
            result.append(symbol_mangled);
        }
    } else {
        result.append(symbol_mangled);
    }
#endif
    return result;
}

void StackTrace::print(std::ostream& os) {
#ifdef __linux__
    const int size = 100;
    void* buffer[size];
    int nptrs = backtrace(buffer, size);
    char** strings = backtrace_symbols(buffer, nptrs);
    for (int i = 1; i < nptrs-2; i++) {
        os << demangle(strings[i]) << std::endl;
    }
    free(strings);
#endif
}

}
