// Copyright (c) 2014 Stefan Walk
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: MIT
//
// Author(s)     : Stefan Walk

// Modifications from original library:
//  * changed inclusion protection tag
//  * moved to namespace CGAL::internal::
//  * add parameter "reset_trees" to train() to be able to construct
//    forest with several iterations

#ifndef CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_FOREST_H
#define CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_FOREST_H
#include "common-libraries.hpp"
#include "tree.hpp"
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#if VERBOSE_TREE_PROGRESS
#include <cstdio>
#endif

namespace CGAL { namespace internal {

namespace liblearning {
namespace RandomForest {

template <typename NodeT>
class RandomForest {
public:
    typedef typename NodeT::ParamType ParamType;
    typedef typename NodeT::FeatureType FeatureType;
    typedef Tree<NodeT> TreeType;
    ParamType params;

    std::vector<uint8_t> was_oob_data;
    DataView2D<uint8_t> was_oob;

    boost::ptr_vector< Tree<NodeT> > trees;

    RandomForest() {}
    RandomForest(ParamType const& params) : params(params) {}

    template<typename SplitGenerator>
    void train(DataView2D<FeatureType> samples, 
               DataView2D<int> labels, 
               DataView2D<int> train_sample_idxes, 
               SplitGenerator const& split_generator,
               size_t seed_start = 1,
               bool register_oob = true,
               bool reset_trees = true
               ) 
    {
        if (reset_trees)
          trees.clear();

        params.n_classes  = *std::max_element(&labels(0,0), &labels(0,0)+labels.num_elements()) + 1;
        params.n_features = samples.cols;
        params.n_samples  = samples.rows;

        std::vector<int> sample_idxes;

        if (train_sample_idxes.empty()) {
            // no indexes were passed, generate vector with all indexes
            sample_idxes.resize(params.n_samples);
            for (size_t i_sample = 0; i_sample < params.n_samples; ++i_sample) {
                sample_idxes[i_sample] = i_sample;
            }
        } else {
            // copy indexes
            sample_idxes.assign(&train_sample_idxes(0,0), &train_sample_idxes(0,0)+train_sample_idxes.num_elements());
        }

        size_t n_idxes = sample_idxes.size();
        params.n_in_bag_samples = n_idxes * (1 - params.sample_reduction);

        // Random distribution over indexes
        UniformIntDist dist(0, n_idxes - 1);

        // Store for each sample and each tree if sample was used for tree
        if (register_oob) {
            was_oob_data.assign(n_idxes*params.n_trees, 1);
            was_oob = DataView2D<uint8_t>(&was_oob_data[0], n_idxes, params.n_trees);
        }

        std::size_t nb_trees = trees.size();
        for (size_t i_tree = nb_trees; i_tree < nb_trees + params.n_trees; ++i_tree) {
#if VERBOSE_TREE_PROGRESS
            std::printf("Training tree %zu/%zu, max depth %zu\n", i_tree+1, nb_trees + params.n_trees, params.max_depth);
#endif
            // new tree
            trees.push_back(new TreeType(&params));
            // initialize random generator with sequential seeds (one for each
            // tree)
            RandomGen gen(seed_start + i_tree);
            // Bagging: draw random sample indexes used for this tree
            std::vector<int> in_bag_samples(params.n_in_bag_samples);
            for (size_t i_sample = 0; i_sample < in_bag_samples.size(); ++i_sample) {
                int random_idx = dist(gen);
                in_bag_samples[i_sample] = sample_idxes[random_idx];
                if (register_oob && was_oob(random_idx, i_tree)) {
                    was_oob(random_idx, i_tree) = 0;
                }
            }
#ifdef TREE_GRAPHVIZ_STREAM
            TREE_GRAPHVIZ_STREAM << "digraph Tree {" << std::endl;
#endif
            // Train the tree
            trees.back().train(samples, labels, &in_bag_samples[0], in_bag_samples.size(), split_generator, gen);
#ifdef TREE_GRAPHVIZ_STREAM
            TREE_GRAPHVIZ_STREAM << "}" << std::endl << std::endl;
#endif
        }
    }
    int evaluate(FeatureType const* sample, float* results) {
        // initialize output probabilities to 0
        std::fill_n(results, params.n_classes, 0.0f);
        // accumulate votes of the trees
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            float const* tree_result = trees[i_tree].evaluate(sample);
            for (size_t i_cls = 0; i_cls < params.n_classes; ++i_cls) {
                results[i_cls] += tree_result[i_cls];
            }
        }
        float best_val   = 0.0;
        int   best_class = 0;
        float scale      = 1.0 / trees.size();
        for (size_t i_cls = 0; i_cls < params.n_classes; ++i_cls) {
            // divide by number of trees to normalize probability
            results[i_cls] *= scale;
            // determine best class
            if (results[i_cls] > best_val) {
                best_val = results[i_cls];
                best_class = i_cls;
            }
        }
        return best_class;
    }
#if 0
    float similarity_endnode(float const* sample_1, float const* sample_2) {
        double sum = 0.0;
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            sum += trees[i_tree].similarity_endnode(sample_1, sample_2);
        }
        return sum/trees.size();
    }
    float similarity_path(float const* sample_1, float const* sample_2) {
        double sum = 0.0;
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            sum += trees[i_tree].similarity_path(sample_1, sample_2);
        }
        return sum/trees.size();
    }
#endif
    template <typename Archive>
    void serialize(Archive& ar, unsigned /* version */)
    {
        ar & BOOST_SERIALIZATION_NVP(params);
        ar & BOOST_SERIALIZATION_NVP(trees);
    }
};

}
}

}} // namespace CGAL::internal::

#endif
