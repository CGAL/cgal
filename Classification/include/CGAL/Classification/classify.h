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

#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/property_map.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <mutex>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {


/// \cond SKIP_IN_MANUAL
namespace internal {

  template <typename Classifier, typename LabelIndexRange>
  class Classify_functor
  {
    const Label_set& m_labels;
    const Classifier& m_classifier;
    LabelIndexRange& m_out;
    
  public:

    Classify_functor (const Label_set& labels,
                      const Classifier& classifier,
                      LabelIndexRange& out)
      : m_labels (labels), m_classifier (classifier), m_out (out)
    { }
    
#ifdef CGAL_LINKED_WITH_TBB
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply(s);
    }
#endif // CGAL_LINKED_WITH_TBB
    
    inline void apply (std::size_t s) const
    {
      std::size_t nb_class_best=0; 
      std::vector<float> values;
      m_classifier (s, values);
        
      float val_class_best = 0.f;
      for(std::size_t k = 0; k < m_labels.size(); ++ k)
      {
        if(val_class_best < values[k])
        {
          val_class_best = values[k];
          nb_class_best = k;
        }
      }

      m_out[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);
    }

  };

  template <typename Classifier, typename LabelIndexRange, typename ProbabilitiesRanges>
  class Classify_detailed_output_functor
  {
    const Label_set& m_labels;
    const Classifier& m_classifier;
    LabelIndexRange& m_out;
    ProbabilitiesRanges& m_prob;
    
  public:

    Classify_detailed_output_functor (const Label_set& labels,
                                      const Classifier& classifier,
                                      LabelIndexRange& out,
                                      ProbabilitiesRanges& prob)
      : m_labels (labels), m_classifier (classifier), m_out (out), m_prob (prob)
    { }
    
#ifdef CGAL_LINKED_WITH_TBB
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply(s);
    }
#endif // CGAL_LINKED_WITH_TBB
    
    inline void apply (std::size_t s) const
    {
      std::size_t nb_class_best=0; 
      std::vector<float> values;
      m_classifier (s, values);
        
      float val_class_best = 0.f;
      for(std::size_t k = 0; k < m_labels.size(); ++ k)
      {
        m_prob[k][s] = values[k];
        if(val_class_best < values[k])
        {
          val_class_best = values[k];
          nb_class_best = k;
        }
      }

      m_out[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);
    }

  };

  template <typename Classifier>
  class Classify_functor_local_smoothing_preprocessing
  {
    const Label_set& m_labels;
    const Classifier& m_classifier;
    std::vector<std::vector<float> >& m_values;
    
  public:

    Classify_functor_local_smoothing_preprocessing
    (const Label_set& labels,
     const Classifier& classifier,
     std::vector<std::vector<float> >& values)
      : m_labels (labels), m_classifier (classifier), m_values (values)
    { }

#ifdef CGAL_LINKED_WITH_TBB    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply (s);
    }
#endif

    inline void apply (std::size_t s) const
    {
      std::vector<float> values;
      m_classifier(s, values);
      for(std::size_t k = 0; k < m_labels.size(); ++ k)
        m_values[k][s] = values[k];
    }
  };
  
  template <typename ItemRange, typename ItemMap, typename NeighborQuery, typename LabelIndexRange>
  class Classify_functor_local_smoothing
  {
    const ItemRange& m_input;
    const ItemMap m_item_map;
    const Label_set& m_labels;
    const std::vector<std::vector<float> >& m_values;
    const NeighborQuery& m_neighbor_query;
    LabelIndexRange& m_out;
    
  public:

    Classify_functor_local_smoothing (const ItemRange& input,
                                      ItemMap item_map,
                                      const Label_set& labels,
                                      const std::vector<std::vector<float> >& values,
                                      const NeighborQuery& neighbor_query,
                                      LabelIndexRange& out)
      : m_input (input), m_item_map (item_map), m_labels (labels),
        m_values(values),
        m_neighbor_query (neighbor_query),
        m_out (out)
    { }

#ifdef CGAL_LINKED_WITH_TBB    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply (s);
    }
#endif

    inline void apply (std::size_t s) const
    {
      std::vector<std::size_t> neighbors;
      m_neighbor_query (get (m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

      std::vector<float> mean (m_values.size(), 0.);
      for (std::size_t n = 0; n < neighbors.size(); ++ n)
        for (std::size_t j = 0; j < m_values.size(); ++ j)
          mean[j] += m_values[j][neighbors[n]];

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

      m_out[s] = static_cast<typename LabelIndexRange::iterator::value_type>(nb_class_best);
    }


  };

  template <typename ItemRange, typename ItemMap,
            typename Classifier, typename NeighborQuery,
            typename LabelIndexRange>
  class Classify_functor_graphcut
  {
    const ItemRange& m_input;
    ItemMap m_item_map;
    const Label_set& m_labels;
    const Classifier& m_classifier;
    const NeighborQuery& m_neighbor_query;
    float m_strength;
    const std::vector<std::vector<std::size_t> >& m_indices;
    const std::vector<std::pair<std::size_t, std::size_t> >& m_input_to_indices;
    LabelIndexRange& m_out;

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
    typedef CGAL::internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
    typedef CGAL::internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
    
  public:

    Classify_functor_graphcut (const ItemRange& input,
                               ItemMap item_map,
                               const Label_set& labels,
                               const Classifier& classifier,
                               const NeighborQuery& neighbor_query,
                               float strength,
                               const std::vector<std::vector<std::size_t> >& indices,
                               const std::vector<std::pair<std::size_t, std::size_t> >& input_to_indices,
                               LabelIndexRange& out)
    : m_input (input), m_item_map (item_map), m_labels (labels),
      m_classifier (classifier), m_neighbor_query (neighbor_query),
      m_strength (strength), m_indices (indices), m_input_to_indices (input_to_indices), m_out (out)
    { }

#ifdef CGAL_LINKED_WITH_TBB
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        apply(s);
    }
#endif // CGAL_LINKED_WITH_TBB

    
    inline void apply (std::size_t sub) const
    {
      if (m_indices[sub].empty())
        return;
        
      std::vector<std::pair<std::size_t, std::size_t> > edges;
      std::vector<double> edge_weights;
      std::vector<std::vector<double> > probability_matrix
        (m_labels.size(), std::vector<double>(m_indices[sub].size(), 0.));
      std::vector<std::size_t> assigned_label (m_indices[sub].size());

      for (std::size_t j = 0; j < m_indices[sub].size(); ++ j)
      {
        std::size_t s = m_indices[sub][j];
            
        std::vector<std::size_t> neighbors;

        m_neighbor_query (get(m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

        for (std::size_t i = 0; i < neighbors.size(); ++ i)
          if (sub == m_input_to_indices[neighbors[i]].first
              && j != m_input_to_indices[neighbors[i]].second)
          {
            edges.push_back (std::make_pair (j, m_input_to_indices[neighbors[i]].second));
            edge_weights.push_back (m_strength);
          }

        std::vector<float> values;
        m_classifier(s, values);
        std::size_t nb_class_best = 0;
        float val_class_best = 0.f;
        for(std::size_t k = 0; k < m_labels.size(); ++ k)
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
    
      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, probability_matrix, assigned_label);

      for (std::size_t i = 0; i < assigned_label.size(); ++ i)
        m_out[m_indices[sub][i]] = static_cast<typename LabelIndexRange::iterator::value_type>(assigned_label[i]);
    }

  };

} // namespace internal
  
/// \endcond
  

  /*! 
    \ingroup PkgClassificationMain

    \brief Runs the classification algorithm without any regularization.

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
    internal::Classify_functor<Classifier, LabelIndexRange>
      f (labels, classifier, output);

#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); ++ i)
        f.apply(i);
    }
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
    internal::Classify_detailed_output_functor<Classifier, LabelIndexRange, ProbabilitiesRanges>
      f (labels, classifier, output, probabilities);

#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); ++ i)
        f.apply(i);
    }
  }
  /// \endcond

  /*! 
    \ingroup PkgClassificationMain

    \brief Runs the classification algorithm with a local smoothing.

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
    internal::Classify_functor_local_smoothing_preprocessing<Classifier>
      f1 (labels, classifier, values);
    internal::Classify_functor_local_smoothing<ItemRange, ItemMap, NeighborQuery, LabelIndexRange>
      f2 (input, item_map, labels, values, neighbor_query, output);
    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f1);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f2);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); ++ i)
        f1.apply(i);
      for (std::size_t i = 0; i < input.size(); ++ i)
        f2.apply(i);
    }
  }

  /*! 
    \ingroup PkgClassificationMain

    \brief Runs the classification algorithm with a global
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
      (boost::make_transform_iterator (input.begin(), CGAL::Property_map_to_unary_function<ItemMap>(item_map)),
       boost::make_transform_iterator (input.end(), CGAL::Property_map_to_unary_function<ItemMap>(item_map)));

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

    internal::Classify_functor_graphcut<ItemRange, ItemMap, Classifier, NeighborQuery, LabelIndexRange>
      f (input, item_map, labels, classifier, neighbor_query, strength, indices, input_to_indices, output);
    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, indices.size ()), f);
    }
    else
#endif
    {
      for (std::size_t sub = 0; sub < indices.size(); ++ sub)
        f.apply (sub);
    }
  }


}

}

#endif // CGAL_CLASSIFICATION_CLASSIFY_H
