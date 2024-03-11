/**
 * @file   util/StringFactory.cpp
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#include "util/StringFactory.h"

#include <boost/date_time/local_time/local_time.hpp>
#include <sstream>

namespace util {

const std::string StringFactory::DATE_FORMAT = "%Y-%m-%d_%H%M%S";

StringFactory::~StringFactory() {
    // intentionally does nothing
}

std::string StringFactory::fromBoolean(bool value) {
    std::string result;
    if (value) {
        result = "true";
    } else {
        result = "false";
    }
    return result;
}

std::string StringFactory::fromInteger(int value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

std::string StringFactory::fromFloat(float value) {
    std::stringstream sstr;
    sstr.precision(4);
    sstr << value;
    return sstr.str();
}

std::string StringFactory::fromFloatArr(int length, float value[]) {
    std::string result;
    result += "<";
    for (int i = 0; i < length; i++) {
        if (i > 0) {
            result += ", ";
        }
        result += fromFloat(value[i]);
    }
    result += ">";
    return result;
}

std::string StringFactory::fromDouble(double value) {
    std::stringstream sstr;
    sstr.precision(4);
    sstr << value;
    return sstr.str();
}

std::string StringFactory::fromDoubleArr(int length, double value[]) {
    std::string result;
    result += "<";
    for (int i = 0; i < length; i++) {
        if (i > 0) {
            result += ", ";
        }
        result += fromDouble(value[i]);
    }
    result += ">";
    return result;
}

std::string StringFactory::fromPointer(const void* value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

std::string StringFactory::replaceAll(const std::string& str, const std::string& search, const std::string& replace) {
    std::string result = str;
    size_t pos = result.find(search);
    while (pos != std::string::npos) {
        result.replace(pos, search.size(), replace);
        pos = result.find(search, pos+replace.size());
    }
    return result;
}

std::string StringFactory::now(const std::string& format) {
    boost::local_time::local_time_facet* facet =
            new boost::local_time::local_time_facet(format.c_str());
    std::stringstream date_stream;
    date_stream.imbue(std::locale(date_stream.getloc(), facet));
    date_stream << boost::local_time::local_microsec_clock::local_time(
            boost::local_time::time_zone_ptr());
    return date_stream.str();
//    time_t rawtime;
//    time(&rawtime);
//    struct tm * timeinfo = localtime(&rawtime);
//    char result[256];
//    strftime(result, sizeof(result), format.c_str(), timeinfo);
//    return string(result);
}

}
