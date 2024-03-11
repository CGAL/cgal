/**
 * @file   util/Configuration.cpp
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#include "util/Configuration.h"

#include "debug.h"
#include "util/StringFuncs.h"
#include <cstdlib>
#include <fstream>
#include <set>

namespace util {

ConfigurationSPtr Configuration::instance_ = ConfigurationSPtr();

Configuration::Configuration() {
    // intentionally does nothing
}

Configuration::~Configuration() {
    // intentionally does nothing
}

ConfigurationSPtr Configuration::getInstance() {
    if (!instance_) {
        instance_ = ConfigurationSPtr(new Configuration());
    }
    return instance_;
}

std::string Configuration::findDefaultFilename() {
    std::string name("StraightSkel");
    std::string result = name + ".ini";
    std::string home(getenv("HOME"));
    std::string sysconfdir("/etc");
    std::string filenames[2];
    filenames[0] = home+"/."+name+"/"+name+".ini";
    filenames[1] = sysconfdir+"/"+name+"/"+name+".ini";
    for (unsigned int i = 0; i < 2; i++) {
        std::ifstream input(filenames[i].c_str());
        if (input.is_open()) {
            input.close();
            result = filenames[i];
            break;
        }
    }
    return result;
}

void Configuration::parse(std::istream& input) {
    properties_.clear();
    std::string section;
    std::string line;
    while (std::getline(input, line)) {
        line = StringFuncs::trim(line);
        if (line.empty()) {
            continue;
        }
        if (StringFuncs::startsWith(line, "[") &&
                StringFuncs::endsWith(line, "]")) {
            section = line.substr(1, line.length()-2);
            continue;
        }
        if (!section.empty()) {
            std::size_t pos = line.find("=");
            if (pos != std::string::npos) {
                std::string key = StringFuncs::trim(line.substr(0, pos));
                std::string value = StringFuncs::trim(
                        line.substr(pos+1, line.length()-pos-1));
                std::string mapkey = section + "." + key;
                properties_[mapkey] = value;
            }
        }
    }
//    namespace pod = boost::program_options::detail;
//    std::set<std::string> options;
//    options.insert("*");
//    for (pod::config_file_iterator i(input, options), e; i != e; i++) {
//        properties_[i->string_key] = i->value[0];
//    }
}

bool Configuration::load(const std::string& filename) {
    DEBUG_VAR(filename);
    bool result = false;
    std::ifstream input(filename.c_str());
    if (input.is_open()) {
        parse(input);
        result = true;
        input.close();
    } else {
        DEBUG_VAL("Error: Config file not found.");
    }
    return result;
}

bool Configuration::isLoaded() const {
    bool result = (properties_.size() > 0);
    return result;
}

bool Configuration::contains(const std::string& section, const std::string& key) {
    std::string mapkey = section + "." + key;
    bool result = (properties_.find(mapkey) != properties_.end());
    return result;
}

std::string Configuration::getString(const std::string& section, const std::string& key) {
    std::string result;
    std::string mapkey = section + "." + key;
    if (properties_.find(mapkey) == properties_.end()) {
        // map does not contain this key
        DEBUG_VAL("key=" << mapkey << " not found.");
    } else {
        result = properties_[mapkey];
    }
    return result;
}

int Configuration::getInt(const std::string& section, const std::string& key) {
    int result = 0;
    std::string value = getString(section, key);
    if (value.length() != 0) {
        result = atoi(value.c_str());
    }
    return result;
}

double Configuration::getDouble(const std::string& section, const std::string& key) {
    double result = 0.0;
    std::string value = getString(section, key);
    if (value.length() != 0) {
        result = atof(value.c_str());
    }
    return result;
}

bool Configuration::getBool(const std::string& section, const std::string& key) {
    bool result = false;
    std::string value = getString(section, key);
    if (value.length() != 0) {
        if (value.compare("1") == 0 ||
                value.compare("t") == 0 ||
                value.compare("T") == 0 ||
                value.compare("true") == 0 ||
                value.compare("TRUE") == 0 ||
                value.compare("True") == 0) {
            result = true;
        }
    }
    return result;
}

}
