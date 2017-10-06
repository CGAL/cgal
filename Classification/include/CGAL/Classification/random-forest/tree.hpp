#ifndef LIBLEARNING_RANDOMFOREST_TREE_H
#define LIBLEARNING_RANDOMFOREST_TREE_H 
#include "../dataview.h"
#include "common-libraries.hpp"
#include <boost/serialization/scoped_ptr.hpp>

namespace liblearning {
namespace RandomForest {

template <typename NodeT>
class Tree {
public:
    typedef NodeT NodeType;
    typedef typename NodeType::ParamType ParamType;
    typedef typename NodeType::FeatureType FeatureType;
    boost::scoped_ptr<NodeT> root_node;
    ParamType const* params;

    Tree() : params(0) {}
    Tree(ParamType const* params) :
        params(params)
    { }

    template<typename SplitGenerator>
    void train(DataView2D<FeatureType> samples, 
               DataView2D<int> labels, 
               int* sample_idxes, 
               size_t n_samples, 
               SplitGenerator split_generator,
               RandomGen const& gen
               ) 
    {
        // copy generator
        RandomGen my_gen = gen;
        // initialize root node
        root_node.reset(new NodeT(0, params));
        // train root node (other notes get trained recursively)
        root_node->train(samples, labels, sample_idxes, n_samples, split_generator, my_gen);
    }
    float const* evaluate(FeatureType const* sample) {
        // start with root
        NodeT const* node = root_node.get();
        // split until leaf
        while (node && !node->is_leaf) {
            node = node->split(sample);
        }
        if (!node) {
            // bug if this is reached
            return 0;
        } else {
            return node->votes();
        }
    }
    std::list<NodeT const*> all_nodes() {
        return root_node->get_all_childs();
    }

#if 0
    float similarity_endnode(float const* sample_1, float const* sample_2) {
        NodeT const* node_1 = root_node.get();
        NodeT const* node_2 = root_node.get();
        while (!node_1->is_leaf) {
            node_1 = node_1->split(sample_1);
        }
        while (!node_2->is_leaf) {
            node_2 = node_2->split(sample_2);
        }
        if (node_1 == node_2) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    float similarity_path(float const* sample_1, float const* sample_2) {
        NodeT const* node_1 = root_node.get();
        NodeT const* node_2 = root_node.get();
        int n_common = 0;
        int n_1 = 1;
        int n_2 = 1;

        while (node_1 == node_2) {
            if (node_1->is_leaf) {
                return 1;
            }
            n_common++;
            node_1 = node_1->split(sample_1);
            node_2 = node_2->split(sample_2);
        }
        while (!node_1->is_leaf) {
            node_1 = node_1->split(sample_1);
            n_1++;
        }
        while (!node_2->is_leaf) {
            node_2 = node_2->split(sample_2);
            n_2++;
        }
        return n_common/sqrt((n_common + n_1)*(n_common + n_2));
    }
#endif
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(params);
        ar & BOOST_SERIALIZATION_NVP(root_node);
    }
};

}
}
#endif
