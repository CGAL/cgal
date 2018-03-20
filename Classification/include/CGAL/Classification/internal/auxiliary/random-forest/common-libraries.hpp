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

#ifndef CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_COMMON_LIBRARIES_H
#define CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_COMMON_LIBRARIES_H
#include <algorithm>
#include <numeric>
#include <limits>
#include <list>
#include <boost/version.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/random/mersenne_twister.hpp>
#if BOOST_VERSION >= 104700
#  include <boost/random/uniform_int_distribution.hpp>
#else 
#  include <boost/random/uniform_int.hpp>
#endif
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/unordered_set.hpp>
#include <iostream>
#include <cstdio>

#include "../dataview.h"

namespace CGAL { namespace internal {

namespace liblearning {
namespace RandomForest {

typedef std::vector< std::pair<float, int> > FeatureClassDataFloat;
inline void init_feature_class_data(FeatureClassDataFloat& data, int /*n_classes*/, int n_samples)
{
    data.resize(n_samples);
}
typedef boost::unordered_set<int> FeatureSet;

#if BOOST_VERSION >= 104700
typedef boost::random::uniform_int_distribution<> UniformIntDist;
typedef boost::random::normal_distribution<> NormalDist;
typedef boost::random::mt19937 RandomGen;
typedef boost::random::uniform_01<> UnitDist;
#else 
typedef boost::uniform_int<> UniformIntDist;
typedef boost::normal_distribution<> NormalDist;
typedef boost::uniform_01<> UnitDist;
typedef boost::mt19937 RandomGen;
#endif

struct ForestParams { 
    size_t n_classes;
    size_t n_features;
    size_t n_samples;
    size_t n_in_bag_samples;
    size_t max_depth;
    size_t n_trees;
    size_t min_samples_per_node;
    float  sample_reduction;
    ForestParams() :
        n_classes(0),
        n_features(0),
        n_samples(0),
        n_in_bag_samples(0),
        max_depth(42),
        n_trees(100),
        min_samples_per_node(5),
        sample_reduction(0)
    {}
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(n_classes);
        ar & BOOST_SERIALIZATION_NVP(n_features);
        ar & BOOST_SERIALIZATION_NVP(n_samples);
        ar & BOOST_SERIALIZATION_NVP(n_in_bag_samples);
        ar & BOOST_SERIALIZATION_NVP(max_depth);
        ar & BOOST_SERIALIZATION_NVP(n_trees);
        ar & BOOST_SERIALIZATION_NVP(min_samples_per_node);
        ar & BOOST_SERIALIZATION_NVP(sample_reduction);
    }
};

struct QuadraticSplitter {
    typedef float FeatureType;
    typedef FeatureClassDataFloat FeatureClassData;
    int n_features;
    std::vector<FeatureType> w;
    FeatureType threshold;
    QuadraticSplitter() : n_features(0) {}
    QuadraticSplitter(int n_features, std::vector<FeatureType> const& w) :
        n_features(n_features), w(w)
    {}
    void set_threshold(FeatureType new_threshold) {
        threshold = new_threshold;
    }
    FeatureType map_sample(FeatureType const* v) const {
        double result = 0.0;
        int i_feature = 0;
        for (; i_feature < n_features; ++i_feature) {
            result += w[i_feature] * v[i_feature];
        }
        for (int i_1 = 0; i_1 < n_features; ++i_1) {
            for (int i_2 = 0; i_2 < n_features; ++i_2) {
                result += w[i_feature++] * v[i_1] * v[i_2];
            }
        }
        return result;
    }
    bool classify_sample(FeatureType const* v) const {
        return map_sample(v) > threshold;
    }
    void map_points(DataView2D<FeatureType> samples,
                    DataView2D<int>   labels,
                    int const*      sample_idxes,
                    int             n_samples,
                    FeatureClassData&        data_points) const
    {
        for (int i_sample = 0; i_sample < n_samples; ++i_sample) {
            int sample_idx    = sample_idxes[i_sample];
            int sample_class  = labels(sample_idx, 0);
            FeatureType sample_fval = map_sample(samples.row_pointer(sample_idx));
            data_points[i_sample] = std::make_pair(sample_fval, sample_class);
        }
    }
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(n_features);
        ar & BOOST_SERIALIZATION_NVP(w);
        ar & BOOST_SERIALIZATION_NVP(threshold);
    }
};

struct LinearSplitter {
    typedef float FeatureType;
    typedef FeatureClassDataFloat FeatureClassData;
    std::vector<FeatureType> w;
    FeatureType threshold;
    LinearSplitter() {}
    LinearSplitter(std::vector<FeatureType> const& w) :
        w(w)
    {}
    void set_threshold(FeatureType new_threshold) {
        threshold = new_threshold;
    }
    bool classify_sample(FeatureType const* v) const {
        return std::inner_product(w.begin(), w.end(), v, 0.0f) > threshold;
    }
    void map_points(DataView2D<FeatureType> samples,
                    DataView2D<int>   labels,
                    int const*      sample_idxes,
                    int             n_samples,
                    FeatureClassData&        data_points) const
    {
        for (int i_sample = 0; i_sample < n_samples; ++i_sample) {
            int sample_idx    = sample_idxes[i_sample];
            int sample_class  = labels(sample_idx, 0);
            FeatureType sample_fval = std::inner_product(w.begin(), w.end(), 
                                                   samples.row_pointer(sample_idx), 0.0f);
            data_points[i_sample] = std::make_pair(sample_fval, sample_class);
        }
    }
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(w);
        ar & BOOST_SERIALIZATION_NVP(threshold);
    }
};

struct AxisAlignedSplitter {
    typedef float FeatureType;
    typedef FeatureClassDataFloat FeatureClassData;
    int feature;
    FeatureType threshold;
    AxisAlignedSplitter() : feature(-1) {}
    AxisAlignedSplitter(int feature) : 
        feature(feature) 
    {}
    void set_threshold(FeatureType new_threshold) {
        threshold = new_threshold;
    }
    bool classify_sample(FeatureType const* v) const {
        return v[feature] > threshold;
    }
    void map_points(DataView2D<FeatureType> samples,
                    DataView2D<int>   labels,
                    int const*      sample_idxes,
                    int             n_samples,
                    FeatureClassData&        data_points) const
    {
        for (int i_sample = 0; i_sample < n_samples; ++i_sample) {
            // determine index of this sample ...
            int sample_idx    = sample_idxes[i_sample];
            // determine class ...
            int sample_class  = labels(sample_idx, 0);
            // determine value of the selected feature for this sample
            FeatureType sample_fval = samples(sample_idx, feature);
            data_points[i_sample] = std::make_pair(sample_fval, sample_class);
        }
    }
    template <typename Archive>
    void serialize(Archive& ar, unsigned /*version*/)
    {
        ar & BOOST_SERIALIZATION_NVP(feature);
        ar & BOOST_SERIALIZATION_NVP(threshold);
    }
};

struct AxisAlignedRandomSplitGenerator {
    typedef float FeatureType;
    FeatureSet features;
    FeatureSet::const_iterator it;

    void init(DataView2D<FeatureType> samples,
              DataView2D<int>   /*labels*/,
              int*            /*sample_idxes*/,
              int             /*n_samples*/,
              size_t          /*n_classes*/,
              RandomGen&      gen)
    {
        features.clear();
        int n_features = samples.cols;
        size_t n_used_features = std::sqrt(n_features);
        UniformIntDist dist(0, n_features - 1);
        // insert into set until required number of unique features is found
        while (features.size() < n_used_features) {
            features.insert(dist(gen));
        }
        it = features.begin();
    }
    AxisAlignedSplitter gen_proposal(RandomGen& /*gen*/)
    {
        if (it == features.end()) {
            it = features.begin();
        }
        return AxisAlignedSplitter(*it++);
    }
    size_t num_proposals() const {
        return features.size();
    }
};

struct LinearSplitGenerator {
    typedef float FeatureType;
    size_t n_features;
    size_t n_proposals;
    LinearSplitGenerator(size_t n_proposals = 5) : 
        n_proposals(n_proposals)
    {}
    void init(DataView2D<FeatureType> samples,
              DataView2D<int>   /*labels*/,
              int*            /*sample_idxes*/,
              int             /*n_samples*/,
              size_t          /*n_classes*/,
              RandomGen&      /*gen*/)
    {
        n_features = samples.cols;
    }
    size_t num_proposals() const {
        return n_proposals;
    }
    LinearSplitter gen_proposal(RandomGen& gen) {
        NormalDist dist;
        std::vector<FeatureType> weights(n_features);
        for (size_t i_feature = 0; i_feature < n_features; ++i_feature) {
            weights[i_feature] = dist(gen);
        }
        return LinearSplitter(weights);
    }
};

struct QuadraticSplitGenerator {
    typedef float FeatureType;
    size_t n_features;
    size_t n_proposals;
    QuadraticSplitGenerator(size_t n_proposals = 5) : 
        n_proposals(n_proposals)
    {}
    void init(DataView2D<FeatureType> samples,
              DataView2D<int>   /*labels*/,
              int*            /*sample_idxes*/,
              int             /*n_samples*/,
              size_t          /*n_classes*/,
              RandomGen&      /*gen*/)
    {
        n_features = samples.cols;
    }
    size_t num_proposals() const {
        return n_proposals;
    }
    QuadraticSplitter gen_proposal(RandomGen& gen) {
        NormalDist dist;
        std::vector<FeatureType> weights(n_features + n_features*n_features);
        for (size_t i_feature = 0; i_feature < weights.size(); ++i_feature) {
            weights[i_feature] = dist(gen);
        }
        return QuadraticSplitter(n_features, weights);
    }
};

}
}

}} // namespace CGAL::internal::
    
#endif
