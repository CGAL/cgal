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

//  * improve sorting algorithm by only comparing the first of pair
//    (second is useless)



#ifndef CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_NODE_GINI_H
#define CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_NODE_GINI_H
#include "node.hpp"
#include "common-libraries.hpp"

namespace CGAL { namespace internal {

namespace liblearning {
namespace RandomForest {

/*
template <typename T>
class X : Y<X> {}
-> https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern#Static_polymorphism
*/

template <typename Splitter>
class NodeGini : public Node< NodeGini<Splitter>, ForestParams, Splitter > {
public:
    typedef typename Node< NodeGini<Splitter>, ForestParams, Splitter>::ParamType ParamType;
    typedef typename Splitter::FeatureType FeatureType;
    typedef typename Splitter::FeatureClassData FeatureClassData;
    using Node< NodeGini<Splitter>, ForestParams, Splitter>::params;
    NodeGini() {}
    NodeGini(int depth, ParamType const* params) :
        Node< NodeGini<Splitter>, ForestParams, Splitter>(depth, params)
    {
    }

    uint64_t gini_square_term(std::vector<uint64_t> const& frequencies) const
    {
        return std::inner_product( frequencies.begin(), frequencies.end(), frequencies.begin(), uint64_t(0));
    }
    std::pair<FeatureType, float> determine_best_threshold(FeatureClassData& data_points,
                                                     std::vector<uint64_t>& classes_l,
                                                     std::vector<uint64_t>& classes_r,
                                                     RandomGen&        gen)
    {
        double best_loss = std::numeric_limits<double>::infinity();
        float best_thresh = 0;

        UnitDist fraction_dist;
        classes_l.assign(params->n_classes, 0);
        classes_r.assign(params->n_classes, 0);
        double n_l = 0;
        double n_r = 0;
        for (size_t i_sample = 0; i_sample < data_points.size(); ++i_sample) {
            classes_r[data_points[i_sample].second]++;
            n_r += 1;
        }
        // sort data so thresholding is easy based on position in array
        std::sort(data_points.begin(), data_points.end(),
                  [&](const std::pair<float, int>& a,
                      const std::pair<float, int>& b) -> bool
                  {
                    return a.first < b.first;
                  });

        // loop over data, update class distributions left&right
        for (size_t i_point = 1; i_point < data_points.size(); ++i_point) {
            int cls = data_points[i_point-1].second;
            classes_l[cls]++; // sample with class cls gets moved to left ...
            classes_r[cls]--; // from right
            n_l += 1;
            n_r -= 1;
            // don't split here if values are the same
            if (data_points[i_point-1].first == data_points[i_point].first)
                continue;
            // weighted average
            double gini = n_l - gini_square_term(classes_l) * 1.0 / n_l + n_r - gini_square_term(classes_r) * 1.0 / n_r;
            if (gini < best_loss) {
                best_loss = gini;
                double fraction = fraction_dist(gen);
                best_thresh = fraction * data_points[i_point-1].first + (1-fraction) * data_points[i_point].first;
            }
        }
        return std::make_pair(best_thresh, float(best_loss));
    }

#if defined(CGAL_LINKED_WITH_BOOST_IOSTREAMS) && defined(CGAL_LINKED_WITH_BOOST_SERIALIZATION)
    template <typename Archive>
    void serialize(Archive& ar, unsigned /* version */)
    {
        ar & boost::serialization::make_nvp("base",  boost::serialization::base_object< Node< NodeGini<Splitter>, ForestParams, Splitter > >(*this));
    }
#endif
};

}
}

}} // namespace CGAL::internal::

#endif
