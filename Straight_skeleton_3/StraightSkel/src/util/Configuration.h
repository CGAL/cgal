/**
 * @file   util/Configuration.h
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#ifndef UTIL_CONFIGURATION_H
#define UTIL_CONFIGURATION_H

#include "util/ptrs.h"
#include <string>
#include <iostream>
#include <map>

namespace util {

class Configuration {
public:
    virtual ~Configuration();
    static ConfigurationSPtr getInstance();

    std::string findDefaultFilename();

    void parse(std::istream& input);
    bool load(const std::string& filename);
    bool isLoaded() const;

    bool contains(const std::string& section, const std::string& key);

    std::string getString(const std::string& section, const std::string& key);
    int getInt(const std::string& section, const std::string& key);
    double getDouble(const std::string& section, const std::string& key);
    bool getBool(const std::string& section, const std::string& key);

protected:
    Configuration();
    static ConfigurationSPtr instance_;

    std::map<std::string, std::string> properties_;
};

}

#endif /* UTIL_CONFIGURATION_H */
