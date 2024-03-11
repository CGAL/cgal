#include <boost/test/unit_test.hpp>

#include "algo/3d/CombiVertexSplitter.h"

BOOST_AUTO_TEST_SUITE(CombiVertexSplitterTest)

using std::list;
using std::vector;
using algo::_3d::vec2i;
using algo::_3d::combi;
using algo::_3d::CombiVertexSplitter;

BOOST_AUTO_TEST_CASE(testCompareSplits) {
    vec2i split1 = CombiVertexSplitter::createSplit(1, 5);
    vec2i split2 = CombiVertexSplitter::createSplit(2, 4);
    BOOST_CHECK_EQUAL(1, CombiVertexSplitter::compareSplits(split1, split2));
    BOOST_CHECK_EQUAL(-1, CombiVertexSplitter::compareSplits(split2, split1));
    BOOST_CHECK_EQUAL(0, CombiVertexSplitter::compareSplits(split1, split1));
}

BOOST_AUTO_TEST_CASE(testSplitLabels) {
    vector<int> labels1 = CombiVertexSplitter::initLabels(5);
    vec2i split = CombiVertexSplitter::createSplit(2, 4);
    vector<int> labels2 = CombiVertexSplitter::splitLabels(labels1, split);
    BOOST_CHECK_EQUAL(labels1.size(), 4);
    if (labels1.size() == 4) {
        BOOST_CHECK_EQUAL(0, labels1[0]);
        BOOST_CHECK_EQUAL(1, labels1[1]);
        BOOST_CHECK_EQUAL(2, labels1[2]);
        BOOST_CHECK_EQUAL(4, labels1[3]);
    }
    BOOST_CHECK_EQUAL(labels2.size(), 3);
    if (labels2.size() == 3) {
        BOOST_CHECK_EQUAL(2, labels2[0]);
        BOOST_CHECK_EQUAL(3, labels2[1]);
        BOOST_CHECK_EQUAL(4, labels2[2]);
    }
}

BOOST_AUTO_TEST_CASE(testCreateSingleSplitCombinations) {
    vector<int> labels = CombiVertexSplitter::initLabels(6);
    list<vec2i> combinations = CombiVertexSplitter::createSingleSplitCombinations(labels);
    BOOST_CHECK_EQUAL(combinations.size(), 9);
    if (combinations.size() == 9) {
        unsigned int i = 0;
        list<vec2i>::iterator it_combis = combinations.begin();
        while (it_combis != combinations.end()) {
            vec2i split = *it_combis++;
            if (i == 0) {
                BOOST_CHECK_EQUAL(0, split[0]);
                BOOST_CHECK_EQUAL(2, split[1]);
            } else if (i == 1) {
                BOOST_CHECK_EQUAL(0, split[0]);
                BOOST_CHECK_EQUAL(3, split[1]);
            } else if (i == 2) {
                BOOST_CHECK_EQUAL(0, split[0]);
                BOOST_CHECK_EQUAL(4, split[1]);
            } else if (i == 3) {
                BOOST_CHECK_EQUAL(1, split[0]);
                BOOST_CHECK_EQUAL(3, split[1]);
            } else if (i == 4) {
                BOOST_CHECK_EQUAL(1, split[0]);
                BOOST_CHECK_EQUAL(4, split[1]);
            } else if (i == 5) {
                BOOST_CHECK_EQUAL(1, split[0]);
                BOOST_CHECK_EQUAL(5, split[1]);
            } else if (i == 6) {
                BOOST_CHECK_EQUAL(2, split[0]);
                BOOST_CHECK_EQUAL(4, split[1]);
            } else if (i == 7) {
                BOOST_CHECK_EQUAL(2, split[0]);
                BOOST_CHECK_EQUAL(5, split[1]);
            } else if (i == 8) {
                BOOST_CHECK_EQUAL(3, split[0]);
                BOOST_CHECK_EQUAL(5, split[1]);
            }
            i++;
        }
    }
}

BOOST_AUTO_TEST_CASE(testGenerateAllCombinations) {
    // number of results equals catalan number
    // each label has to occur as often as each other label
    for (unsigned int degree = 4; degree < 9; degree++) {
        list<combi> result = CombiVertexSplitter::generateAllCombinations(degree);
        int num_labels = (result.size() * (2*(degree-3)));
        BOOST_CHECK_EQUAL(0, num_labels%degree);
    }
}

BOOST_AUTO_TEST_SUITE_END()
