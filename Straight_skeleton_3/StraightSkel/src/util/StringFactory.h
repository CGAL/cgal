/**
 * @file   util/StringFactory.h
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#ifndef UTIL_STRINGFACTORY_H
#define UTIL_STRINGFACTORY_H

#include <string>

namespace util {

class StringFactory {
public:
    virtual ~StringFactory();
    static std::string fromBoolean(bool value);
    static std::string fromInteger(int value);
    static std::string fromFloat(float value);
    static std::string fromFloatArr(int length, float value[]);
    static std::string fromDouble(double value);
    static std::string fromDoubleArr(int length, double value[]);
    static std::string fromPointer(const void* value);
    static std::string replaceAll(const std::string& str, const std::string& search, const std::string& replace);

    static const std::string DATE_FORMAT;

    /**
     * @param format  %Y-%m-%d_%H%M%S
     */
    static std::string now(const std::string& format);
};

}

#endif /* UTIL_STRINGFACTORY_H */
