// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_CLASSIFY_H
#define CGAL_CLASSIFICATION_CLASSIFY_H

#include <CGAL/license/Classification.h>

#include <CGAL/boost/graph/alpha_expansion_graphcut.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/for_each.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/property_map.h>
#include <CGAL/iterator.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <mutex>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationMain

    \brief runs the classification algorithm without any regularization.

    There is no relationship between items, the classification energy
    is only minimized itemwise. This method is quick but produces
    suboptimal results.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_if_available_tag`, `Parallel_tag` or `Sequential_tag`.

    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`. Its value type depends on the data that is
    classified (for example, `CGAL::Point_3` or `CGAL::Triangle_3`).

    \tparam Classifier model of `Classifier`.

    \tparam Model of `Range` with random access iterators whose value
    type is an integer type.

    \param input input range.
    \param labels set of input labels.
    \param classifier input classifier.
    \param output where to store the result. It is stored as a sequence,
    ordered like the input range, containing for each point the index
    (in the `Label_set`) of the assigned label.

  */
  template <typename ConcurrencyTag,
            typename ItemRange,
            typename Classifier,
            typename LabelIndexRange>
  void classify (const ItemRange& input,
                 const Label_set& labels,
                 const Classifier& classifier,
                 LabelIndexRange& output)
  {
    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_counting_range<std::size_t> (0, input.size()),
       [&](const std::size_t& s) -> bool
       {
         std::size_t nb_class_best=0;
         std::vector<float> values;
         classifier (s, values);

         float val_class_best = 0.f;
         for(std::size_t k = 0; k < labels.size(); ++ k)
         {
           if(val_class_best < values[k])
           {
             val_class_best = values[k];
             nb_class_best = k;
           }
         }
         output[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);

         return true;
       });
  }

  /// \cond SKIP_IN_MANUAL
  // variant to get a detailed output (not documented yet)
  template <typename ConcurrencyTag,
            typename ItemRange,
            typename Classifier,
            typename LabelIndexRange,
            typename ProbabilitiesRanges>
  void classify (const ItemRange& input,
                 const Label_set& labels,
                 const Classifier& classifier,
                 LabelIndexRange& output,
                 ProbabilitiesRanges& probabilities)
  {
    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_counting_range<std::size_t> (0, input.size()),
       [&](const std::size_t& s) -> bool
       {
         std::size_t nb_class_best=0;
         std::vector<float> values;
         classifier (s, values);

         float val_class_best = 0.f;
         for(std::size_t k = 0; k < labels.size(); ++ k)
         {
           probabilities[k][s] = values[k];
           if(val_class_best < values[k])
           {
             val_class_best = values[k];
             nb_class_best = k;
           }
         }

         output[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);

         return true;
       });
  }
  /// \endcond

  /*!
    \ingroup PkgClassificationMain

    \brief runs the classification algorithm with a local smoothing.

    The computed classification energy is smoothed on a user defined
    local neighborhood of items. This method is a compromise between
    efficiency and better quality results.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_if_available_tag`, `Parallel_tag` or `Sequential_tag`.
    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`.
    \tparam ItemMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `ItemRange` and value type
    is the type of item to classify (for example, `CGAL::Point_3`).
    \tparam NeighborQuery model of `NeighborQuery`.
    \tparam Classifier model of `Classifier`.
    \tparam Model of `Range` with random access iterators whose value
    type is an integer type.

    \param input input range.
    \param item_map property map to access the input items.
    \param labels set of input labels.
    \param classifier input classifier.
    \param neighbor_query used to access neighborhoods of items.
    \param output where to store the result. It is stored as a sequence,
    ordered like the input range, containing for each point the index
    (in the `Label_set`) of the assigned label.

  */
  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ItemMap,
            typename NeighborQuery,
            typename Classifier,
            typename LabelIndexRange>
  void classify_with_local_smoothing (const ItemRange& input,
                                      const ItemMap item_map,
                                      const Label_set& labels,
                                      const Classifier& classifier,
                                      const NeighborQuery& neighbor_query,
                                      LabelIndexRange& output)
  {
    std::vector<std::vector<float> > values
      (labels.size(), std::vector<float> (input.size(), -1.));

    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_counting_range<std::size_t> (0, input.size()),
       [&](const std::size_t& s) -> bool
       {
         std::vector<float> v;
         classifier(s, v);
         for(std::size_t k = 0; k < labels.size(); ++ k)
           values[k][s] = v[k];

         return true;
       });

    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_counting_range<std::size_t> (0, input.size()),
       [&](const std::size_t& s) -> bool
       {
         std::vector<std::size_t> neighbors;
         neighbor_query (get (item_map, *(input.begin()+s)), std::back_inserter (neighbors));

         std::vector<float> mean (values.size(), 0.);
         for (std::size_t n = 0; n < neighbors.size(); ++ n)
           for (std::size_t j = 0; j < values.size(); ++ j)
             mean[j] += values[j][neighbors[n]];

         std::size_t nb_class_best=0;
         float val_class_best = 0.f;
         for(std::size_t k = 0; k < mean.size(); ++ k)
         {
           mean[k] /= neighbors.size();
           if(val_class_best < mean[k])
           {
             val_class_best = mean[k];
             nb_class_best = k;
           }
         }

         output[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);

         return true;
       });
  }

  /*!
    \ingroup PkgClassificationMain

    \brief runs the classification algorithm with a global
    regularization based on a graph cut.

    The computed classification energy is globally regularized through
    an alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    To speed up computation, the input domain can be subdivided into
    smaller subsets such that several smaller graph cuts are applied
    instead of a big one. The computation of these smaller graph cuts can
    be done in parallel. Increasing the number of subsets allows for
    faster computation times but can also reduce the quality of the
    results.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_if_available_tag`, `Parallel_tag` or `Sequential_tag`.
    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`.
    \tparam ItemMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `ItemRange` and value type
    is the type of item to classify (for example, `CGAL::Point_3`).
    \tparam NeighborQuery model of `NeighborQuery`.
    \tparam Classifier model of `Classifier`.
    \tparam Model of `Range` with random access iterators whose value
    type is an integer type.

    \param input input range.
    \param item_map property map to access the input items.
    \param labels set of input labels.
    \param classifier input classifier.
    \param neighbor_query used to access neighborhoods of items.
    \param strength strength of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.
    \param min_number_of_subdivisions minimum number of subdivisions
    (for parallel processing to be efficient, this should be at least
    the number of cores of the processor).
    \param output where to store the result. It is stored as a sequence,
    ordered like the input range, containing for each point the index
    (in the `Label_set`) of the assigned label.

  */
  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ItemMap,
            typename NeighborQuery,
            typename Classifier,
            typename LabelIndexRange>
  void classify_with_graphcut (const ItemRange& input,
                               const ItemMap item_map,
                               const Label_set& labels,
                               const Classifier& classifier,
                               const NeighborQuery& neighbor_query,
                               const float strength,
                               const std::size_t min_number_of_subdivisions,
                               LabelIndexRange& output)
  {
    CGAL::Bbox_3 bbox = CGAL::bbox_3
      (CGAL::make_transform_iterator_from_property_map (input.begin(), item_map),
       CGAL::make_transform_iterator_from_property_map (input.end(), item_map));

    double Dx = double(bbox.xmax() - bbox.xmin());
    double Dy = double(bbox.ymax() - bbox.ymin());
    double A = Dx * Dy;
    double a = A / min_number_of_subdivisions;
    double l = std::sqrt(a);
    std::size_t nb_x = std::size_t(Dx / l) + 1;
    std::size_t nb_y = std::size_t((A / nb_x) / a) + 1;
    std::size_t nb = nb_x * nb_y;

    std::vector<CGAL::Bbox_3> bboxes;
    bboxes.reserve(nb);
    for (std::size_t x = 0; x < nb_x; ++ x)
      for (std::size_t y = 0; y < nb_y; ++ y)
      {
        bboxes.push_back
          (CGAL::Bbox_3 (bbox.xmin() + Dx * (x / double(nb_x)),
                         bbox.ymin() + Dy * (y / double(nb_y)),
                         bbox.zmin(),
                         (x == nb_x - 1 ? bbox.xmax() : bbox.xmin() + Dx * ((x+1) / double(nb_x))),
                         (y == nb_y - 1 ? bbox.ymax() : bbox.ymin() + Dy * ((y+1) / double(nb_y))),
                         bbox.zmax()));
      }

#ifdef CGAL_CLASSIFICATION_VERBOSE
    std::cerr << "Number of divisions = " << nb_x * nb_y << std::endl;
    std::cerr << " -> Size of division: " << Dx / nb_x << " " << Dy / nb_y << std::endl;
#endif

    std::vector<std::vector<std::size_t> > indices (nb);
    std::vector<std::pair<std::size_t, std::size_t> > input_to_indices(input.size());

    for (std::size_t s = 0; s < input.size(); ++ s)
    {
      CGAL::Bbox_3 b = get(item_map, *(input.begin() + s)).bbox();

      std::size_t i = 0;
      for (; i < bboxes.size(); ++ i)
        if (CGAL::do_overlap (b, bboxes[i]))
        {
          input_to_indices[s] = std::make_pair (i, indices[i].size());
          indices[i].push_back (s);
          break;
        }
      CGAL_assertion_msg (i != bboxes.size(), "Point was not assigned to any subdivision.");
    }

    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_counting_range<std::size_t> (0, indices.size()),
       [&](const std::size_t& sub) -> bool
       {
         if (indices[sub].empty())
           return true;

         std::vector<std::pair<std::size_t, std::size_t> > edges;
         std::vector<double> edge_weights;
         std::vector<std::vector<double> > probability_matrix
           (labels.size(), std::vector<double>(indices[sub].size(), 0.));
         std::vector<std::size_t> assigned_label (indices[sub].size());

         for (std::size_t j = 0; j < indices[sub].size(); ++ j)
         {
           std::size_t s = indices[sub][j];

           std::vector<std::size_t> neighbors;

           neighbor_query (get(item_map, *(input.begin()+s)), std::back_inserter (neighbors));

           for (std::size_t i = 0; i < neighbors.size(); ++ i)
             if (sub == input_to_indices[neighbors[i]].first
                 && j != input_to_indices[neighbors[i]].second)
             {
               edges.push_back (std::make_pair (j, input_to_indices[neighbors[i]].second));
               edge_weights.push_back (strength);
             }

           std::vector<float> values;
           classifier(s, values);
           std::size_t nb_class_best = 0;
           float val_class_best = 0.f;
           for(std::size_t k = 0; k < labels.size(); ++ k)
           {
             float value = values[k];
             probability_matrix[k][j] = -std::log(value);

             if(val_class_best < value)
             {
               val_class_best = value;
               nb_class_best = k;
             }
           }
           assigned_label[j] = nb_class_best;
         }

         CGAL::alpha_expansion_graphcut (edges, edge_weights, probability_matrix, assigned_label);

         for (std::size_t i = 0; i < assigned_label.size(); ++ i)
           output[indices[sub][i]] = static_cast<typename LabelIndexRange::iterator::value_type>(assigned_label[i]);
         return true;
       });
  }


}

}

#endif // CGAL_CLASSIFICATION_CLASSIFY_H
