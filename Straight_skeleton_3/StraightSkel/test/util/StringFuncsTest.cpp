#include <boost/test/unit_test.hpp>

#include "util/StringFuncs.h"

using std::string;
using std::vector;
using util::StringFuncs;

BOOST_AUTO_TEST_SUITE(StringFuncsTest)

BOOST_AUTO_TEST_CASE(testStartsWith) {
    string str("abcdefg");
    string prefix("abc");
    bool result = StringFuncs::startsWith(str, prefix);
    BOOST_CHECK_EQUAL(true, result);
}

BOOST_AUTO_TEST_CASE(testEndsWith) {
    string str("README.txt");
    string suffix("txt");
    bool result = StringFuncs::endsWith(str, suffix);
    BOOST_CHECK_EQUAL(true, result);

    suffix = "pdf";
    result = StringFuncs::endsWith(str, suffix);
    BOOST_CHECK_EQUAL(false, result);
}

BOOST_AUTO_TEST_CASE(testTrim) {
    string str("  as df ");
    string result = StringFuncs::trim(str);
    string expected("as df");
    BOOST_CHECK_EQUAL(0, expected.compare(result));

    str = "   ";
    result = StringFuncs::trim(str);
    BOOST_CHECK_EQUAL(true, result.empty());
}

BOOST_AUTO_TEST_CASE(testSplit) {
    string str("ab cd   ef");
    vector<string> result = StringFuncs::split(str, " \t", false);
    BOOST_CHECK_EQUAL(3, result.size());
    BOOST_CHECK_EQUAL(0, result[0].compare("ab"));
    BOOST_CHECK_EQUAL(0, result[1].compare("cd"));
    BOOST_CHECK_EQUAL(0, result[2].compare("ef"));
    str = "asdf;;jkl;s";
    result = StringFuncs::split(str, ";", true);
    BOOST_CHECK_EQUAL(4, result.size());
    BOOST_CHECK_EQUAL(0, result[0].compare("asdf"));
    BOOST_CHECK_EQUAL(0, result[1].compare(""));
    BOOST_CHECK_EQUAL(0, result[2].compare("jkl"));
    BOOST_CHECK_EQUAL(0, result[3].compare("s"));
}

BOOST_AUTO_TEST_SUITE_END()
