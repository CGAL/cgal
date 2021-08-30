// Copyright (c) 2014 Stefan Walk
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LicenseRef-RFL
// License notice in Installation/LICENSE.RFL
//
// Author(s)     : Stefan Walk

// Modifications from original library:
//  * changed inclusion protection tag
//  * moved to namespace CGAL::internal::
//  * fix computation of node_dist[label] so that results are always <= 1.0
//  * change serialization functions to avoid a bug with boost and some
//    compilers (that leads to dereferencing a null pointer)
//  * add a method to get feature usage

#ifndef CGAL_INTERNAL_LIBLEARNING_RANDOMFORESTS_NODE_H
#define CGAL_INTERNAL_LIBLEARNING_RANDOMFORESTS_NODE_H
#include "../dataview.h"
#include "common-libraries.hpp"

#include <CGAL/IO/binary_file_io.h>

#if defined(CGAL_LINKED_WITH_BOOST_IOSTREAMS) && defined(CGAL_LINKED_WITH_BOOST_SERIALIZATION)
#include <boost/serialization/scoped_ptr.hpp>
#include <boost/serialization/vector.hpp>
#else
#include <boost/scoped_ptr.hpp>
#include <vector>
#endif

#if VERBOSE_NODE_LEARNING
#include <cstdio>
#endif

namespace CGAL { namespace internal {

namespace liblearning {
namespace RandomForest {

template <typename Derived, typename ParamT, typename Splitter>
class Node {
public:
    typedef typename Splitter::FeatureType FeatureType;
    bool is_leaf;
    size_t n_samples;
    size_t depth;
    typedef ParamT ParamType;
    ParamType const* params;
    Splitter splitter;

    boost::scoped_ptr<Derived> left;
    boost::scoped_ptr<Derived> right;
    std::vector<float> node_dist;

    Node() : is_leaf(true), n_samples(0), depth(-1), params(0) {}
    Node(size_t depth, ParamType const* params) :
        is_leaf(true), n_samples(0), depth(depth), params(params)
    {}

    bool pure(DataView2D<int> labels, int* sample_idxes) const {
        if (n_samples < 2)
            return true; // an empty node is by definition pure
        int first_sample_idx = sample_idxes[0];
        int seen_class = labels(first_sample_idx, 0);
        // check if all classes are equal to the first class
        for (size_t i_sample = 1; i_sample < n_samples; ++i_sample) {
            int sample_idx = sample_idxes[i_sample];
            if (labels(sample_idx, 0) != seen_class)
                return false;
        }
        return true;
    }

    float const* votes() const {
        return (float const*)&node_dist[0];
    }

    int partition_samples(DataView2D<FeatureType> samples, int* sample_idxes) {
        // sort samples in bag so that left-samples precede right-samples
        // works like std::partition
        int low  = 0;
        int high = n_samples;

        while (true) {
            while (true) {
                if (low == high) {
                    return low;
                } else if (!splitter.classify_sample(samples.row_pointer(sample_idxes[low]))) {
                    ++low;
                } else {
                    break;
                }
            }
            --high;
            while (true) {
                if (low == high) {
                    return low;
                } else if (splitter.classify_sample(samples.row_pointer(sample_idxes[high]))) {
                    --high;
                } else {
                    break;
                }
            }
            std::swap(sample_idxes[low], sample_idxes[high]);
            ++low;
        }
    }

    Derived const* split (FeatureType const* sample) const {
        if (splitter.classify_sample(sample)) {
            return right.get();
        } else {
            return left.get();
        }
    }

    typedef std::list<Derived const*> NodeList;

    NodeList get_all_childs() {
        NodeList ret;
        ret.push_back(this);
        if (!is_leaf) {
            NodeList left_childs  = left->get_all_childs();
            ret.splice(ret.end(), left_childs);
            NodeList right_childs = right->get_all_childs();
            ret.splice(ret.end(), right_childs);
        }
        return ret;
    }

    template<typename SplitGenerator>
    void determine_best_split(DataView2D<FeatureType> samples,
                              DataView2D<int>         labels,
                              int*                  sample_idxes,
                              SplitGenerator        split_generator,
                              RandomGen&            gen
                              )
    {
        typename Splitter::FeatureClassData data_points;
        init_feature_class_data(data_points, params->n_classes, n_samples);
        float best_loss = std::numeric_limits<float>::infinity();

        std::vector<uint64_t> classes_l;
        std::vector<uint64_t> classes_r;

        // pass information about data to split generator
        split_generator.init(samples,
                             labels,
                             sample_idxes,
                             n_samples,
                             params->n_classes,
                             gen);

        size_t n_proposals = split_generator.num_proposals();

        std::pair<FeatureType, float> results; // (threshold, loss)

        for (size_t i_proposal = 0; i_proposal < n_proposals; ++i_proposal) {
            // generate proposal
            Splitter split = split_generator.gen_proposal(gen);
            // map samples to numbers using proposal
            split.map_points(samples, labels, sample_idxes, n_samples, data_points);
            // check best loss using this proposal
            results = static_cast<Derived*>(this)->determine_best_threshold(data_points, classes_l, classes_r, gen);
            if (results.second < best_loss) {
                // Proposal resulted into new optimum
                best_loss = results.second;
                split.set_threshold(results.first);
                splitter = split;
            }
        }
    }

    template<typename SplitGenerator>
    void train(DataView2D<FeatureType> samples,
               DataView2D<int> labels,
               int* sample_idxes,
               size_t n_samples_,
               SplitGenerator const& split_generator,
               RandomGen& gen
               )
    {
        n_samples = n_samples_;
        node_dist.resize(params->n_classes, 0.0f);
        for (size_t i_sample = 0; i_sample < n_samples; ++i_sample) {
            int label = labels(sample_idxes[i_sample], 0);
            node_dist[label] += 1.0f;
        }

        if (n_samples != 0)
          for (std::size_t i = 0; i < node_dist.size(); ++ i)
            node_dist[i] /= n_samples;

        bool do_split = // Only split if ...
            (n_samples >= params->min_samples_per_node) && // enough samples are available
            !pure(labels, sample_idxes) && // this node is not already pure
            (depth < params->max_depth); // we did not reach max depth
        if (!do_split) {
            splitter.threshold = 0.0;
            return;
        }

        is_leaf = false;

#if VERBOSE_NODE_LEARNING
        std::printf("Determining the best split at depth %zu/%zu\n", depth, params->max_depth);
#endif
        determine_best_split(samples, labels, sample_idxes, split_generator, gen);

        left.reset(new Derived(depth + 1, params));
        right.reset(new Derived(depth + 1, params));

        // sort samples in bag so that left-samples precede right-samples
        int low = partition_samples(samples, sample_idxes);
        int n_samples_left  = low;
        int n_samples_right = n_samples - low;
        int offset_left     = 0;
        int offset_right    = low;
#ifdef TREE_GRAPHVIZ_STREAM
        if (depth <= TREE_GRAPHVIZ_MAX_DEPTH) {
            TREE_GRAPHVIZ_STREAM << "p" << std::hex << (unsigned long)this
                << " -> "
                << "p" << std::hex << (unsigned long)left.get()
                << std::dec << " [label=\"" << n_samples_left << "\"];" <<  std::endl;
            TREE_GRAPHVIZ_STREAM << "p" << std::hex << (unsigned long)this
                << " -> "
                << "p" << std::hex << (unsigned long)right.get()
                << std::dec << " [label=\"" << n_samples_right << "\"];" << std::endl;
        }
#endif
        // train left and right side of split
        left->train (samples, labels, sample_idxes + offset_left,  n_samples_left,  split_generator, gen);
        right->train(samples, labels, sample_idxes + offset_right, n_samples_right, split_generator, gen);
    }

#if defined(CGAL_LINKED_WITH_BOOST_IOSTREAMS) && defined(CGAL_LINKED_WITH_BOOST_SERIALIZATION)
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(is_leaf);
        ar & BOOST_SERIALIZATION_NVP(n_samples);
        ar & BOOST_SERIALIZATION_NVP(depth);
        ar & BOOST_SERIALIZATION_NVP(params);
        ar & BOOST_SERIALIZATION_NVP(splitter);
        ar & BOOST_SERIALIZATION_NVP(node_dist);
        if (!is_leaf)
        {
          ar & BOOST_SERIALIZATION_NVP(left);
          ar & BOOST_SERIALIZATION_NVP(right);
        }
    }
#endif

    void write (std::ostream& os)
    {
      I_Binary_write_bool (os, is_leaf);
      I_Binary_write_size_t_into_uinteger32 (os, n_samples);
      I_Binary_write_size_t_into_uinteger32 (os, depth);
      splitter.write(os);

      for (const float& f : node_dist)
        I_Binary_write_float32 (os, f);

      if (!is_leaf)
      {
        left->write(os);
        right->write(os);
      }
    }

    void read (std::istream& is)
    {
      I_Binary_read_bool (is, is_leaf);
      I_Binary_read_size_t_from_uinteger32 (is, n_samples);
      I_Binary_read_size_t_from_uinteger32 (is, depth);
      splitter.read(is);

      node_dist.resize(params->n_classes, 0.0f);
      for (std::size_t i = 0; i < node_dist.size(); ++ i)
        I_Binary_read_float32 (is, node_dist[i]);

      if (!is_leaf)
      {
        left.reset(new Derived(depth + 1, params));
        right.reset(new Derived(depth + 1, params));
        left->read(is);
        right->read(is);
      }
    }

    void get_feature_usage (std::vector<std::size_t>& count) const
    {
      if (!is_leaf && splitter.feature != -1)
      {
        count[std::size_t(splitter.feature)] ++;
        left->get_feature_usage(count);
        right->get_feature_usage(count);
      }
    }
};

}
}

}} // namespace CGAL::internal::

#endif
