/**
 * @file   util/StringFuncs.cpp
 * @author Gernot Walzl
 * @date   2013-12-16
 */

#include "util/StringFuncs.h"

#include <cstddef>

namespace util {

bool StringFuncs::startsWith(const std::string& str, const std::string& prefix) {
    size_t len_str = str.length();
    size_t len_prefix = prefix.length();
    if (len_prefix > len_str) {
        return false;
    }
    return (0 == str.compare(0, len_prefix, prefix));
}

bool StringFuncs::endsWith(const std::string& str, const std::string& suffix) {
    size_t len_str = str.length();
    size_t len_suffix = suffix.length();
    if (len_suffix > len_str) {
        return false;
    }
    return (0 == str.compare(len_str-len_suffix, len_suffix, suffix));
}

std::string StringFuncs::trim(const std::string& str) {
    std::string result;
    unsigned int length = str.length();
    for (int i = str.length()-1; i >= 0; i--) {
        const char chr = str.at(i);
        if (chr == ' ' || chr == '\t') {
            length = i;
        } else {
            break;
        }
    }
    result = str.substr(0, length);
    unsigned int pos = 0;
    for (unsigned int i = 0; i < length; i++) {
        const char chr = str.at(i);
        if (chr == ' ' || chr == '\t') {
            pos = i+1;
        } else {
            break;
        }
    }
    result = result.substr(pos, length-pos);
    return result;
}

std::vector<std::string> StringFuncs::split(const std::string& str, const std::string& delim, bool keep_empty) {
    std::vector<std::string> result;
    std::string element;
    std::size_t pos = 0;
    std::size_t pos_delim = str.find_first_of(delim);
    while (pos_delim != std::string::npos) {
        element = str.substr(pos, pos_delim - pos);
        if (element.length() > 0 || keep_empty) {
            result.push_back(element);
        }
        pos = pos_delim+1;
        pos_delim = str.find_first_of(delim, pos);
    }
    element = str.substr(pos);
    if (element.length() > 0 || keep_empty) {
        result.push_back(element);
    }
    return result;
}

}
