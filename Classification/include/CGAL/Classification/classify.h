// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_CLASSIFY_H
#define CGAL_CLASSIFICATION_CLASSIFY_H

#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

#include <CGAL/Classification/Label_set.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {

namespace internal {

  template <typename ClassificationPredicate>
  class Classifier
  {
    const Label_set& m_labels;
    const ClassificationPredicate& m_predicate;
    std::vector<std::size_t>& m_out;
    
  public:

    Classifier (const Label_set& labels,
                const ClassificationPredicate& predicate,
                std::vector<std::size_t>& out)
      : m_labels (labels), m_predicate (predicate), m_out (out)
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
      m_predicate.probabilities (s, values);
        
      float val_class_best = (std::numeric_limits<float>::max)();      
      for(std::size_t k = 0; k < m_labels.size(); ++ k)
        {
          if(val_class_best > values[k])
            {
              val_class_best = values[k];
              nb_class_best = k;
            }
        }

      m_out[s] = nb_class_best;
    }

  };

  template <typename ClassificationPredicate>
  class Classifier_local_smoothing_preprocessing
  {
    const Label_set& m_labels;
    const ClassificationPredicate& m_predicate;
    std::vector<std::vector<float> >& m_values;
    
  public:

    Classifier_local_smoothing_preprocessing
    (const Label_set& labels,
     const ClassificationPredicate& predicate,
     std::vector<std::vector<float> >& values)
      : m_labels (labels), m_predicate (predicate), m_values (values)
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
      m_predicate.probabilities(s, values);
      for(std::size_t k = 0; k < m_labels.size(); ++ k)
        m_values[k][s] = values[k];
    }
  };
  
  template <typename ItemRange, typename ItemMap, typename NeighborQuery>
  class Classifier_local_smoothing
  {
    const ItemRange& m_input;
    const ItemMap m_item_map;
    const Label_set& m_labels;
    const std::vector<std::vector<float> >& m_values;
    const NeighborQuery& m_neighbor_query;
    std::vector<std::size_t>& m_out;
    
  public:

    Classifier_local_smoothing (const ItemRange& input,
                                ItemMap item_map,
                                const Label_set& labels,
                                const std::vector<std::vector<float> >& values,
                                const NeighborQuery& neighbor_query,
                                std::vector<std::size_t>& out)
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
      float val_class_best = (std::numeric_limits<float>::max)();
      for(std::size_t k = 0; k < mean.size(); ++ k)
        {
          mean[k] /= neighbors.size();
          if(val_class_best > mean[k])
            {
              val_class_best = mean[k];
              nb_class_best = k;
            }
        }

      m_out[s] = nb_class_best;
    }


  };

  template <typename ItemRange, typename ItemMap,
            typename ClassificationPredicate, typename NeighborQuery>
  class Classifier_graphcut
  {
    const ItemRange& m_input;
    ItemMap m_item_map;
    const Label_set& m_labels;
    const ClassificationPredicate& m_predicate;
    const NeighborQuery& m_neighbor_query;
    float m_weight;
    const std::vector<std::vector<std::size_t> >& m_indices;
    const std::vector<std::pair<std::size_t, std::size_t> >& m_input_to_indices;
    std::vector<std::size_t>& m_out;

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
    typedef CGAL::internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
    typedef CGAL::internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
    
  public:

    Classifier_graphcut (const ItemRange& input,
                         ItemMap item_map,
                         const Label_set& labels,
                         const ClassificationPredicate& predicate,
                         const NeighborQuery& neighbor_query,
                         float weight,
                         const std::vector<std::vector<std::size_t> >& indices,
                         const std::vector<std::pair<std::size_t, std::size_t> >& input_to_indices,
                         std::vector<std::size_t>& out)
      : m_input (input), m_item_map (item_map), m_labels (labels),
        m_predicate (predicate), m_neighbor_query (neighbor_query),
        m_weight (weight), m_indices (indices), m_input_to_indices (input_to_indices), m_out (out)
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
                edge_weights.push_back (m_weight);
              }

          std::vector<float> values;
          m_predicate.probabilities(s, values);
          std::size_t nb_class_best = 0;
          float val_class_best = (std::numeric_limits<float>::max)();
          for(std::size_t k = 0; k < m_labels.size(); ++ k)
            {
              float value = values[k];
              probability_matrix[k][j] = value;
            
              if(val_class_best > value)
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
        m_out[m_indices[sub][i]] = assigned_label[i];
    }

  };

} // namespace internal

  

  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ClassificationPredicate>
  void classify (const ItemRange& input,
                 const Label_set& labels,
                 const ClassificationPredicate& predicate,
                 std::vector<std::size_t>& output)
  {
    output.resize(input.size());
    
    internal::Classifier<ClassificationPredicate>
      f (labels, predicate, output);

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

  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ItemMap,
            typename NeighborQuery,
            typename ClassificationPredicate>
  void classify_with_local_smoothing (const ItemRange& input,
                                      const ItemMap item_map,
                                      const Label_set& labels,
                                      const ClassificationPredicate& predicate,
                                      const NeighborQuery& neighbor_query,
                                      std::vector<std::size_t>& output)
  {
    output.resize(input.size());
    
    std::vector<std::vector<float> > values
      (labels.size(), std::vector<float> (input.size(), -1.));
    internal::Classifier_local_smoothing_preprocessing<ClassificationPredicate>
      f1 (labels, predicate, values);
    internal::Classifier_local_smoothing<ItemRange, ItemMap, NeighborQuery>
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

  template <typename ConcurrencyTag,
            typename ItemRange,
            typename ItemMap,
            typename ItemWithBboxMap,
            typename NeighborQuery,
            typename ClassificationPredicate>
  void classify_with_graphcut (const ItemRange& input,
                               const ItemMap item_map,
                               const ItemWithBboxMap bbox_map,
                               const Label_set& labels,
                               const ClassificationPredicate& predicate,
                               const NeighborQuery& neighbor_query,
                               const float weight,
                               const std::size_t min_number_of_subdivisions,
                               std::vector<std::size_t>& output)
  {
    CGAL::Bbox_3 bbox = CGAL::bbox_3
      (boost::make_transform_iterator (input.begin(), CGAL::Property_map_to_unary_function<ItemWithBboxMap>(bbox_map)),
       boost::make_transform_iterator (input.end(), CGAL::Property_map_to_unary_function<ItemWithBboxMap>(bbox_map)));

    float Dx = bbox.xmax() - bbox.xmin();
    float Dy = bbox.ymax() - bbox.ymin();
    float A = Dx * Dy;
    float a = A / min_number_of_subdivisions;
    float l = std::sqrt(a);
    std::size_t nb_x = std::size_t(Dx / l) + 1;
    std::size_t nb_y = std::size_t((A / nb_x) / a) + 1;
    std::size_t nb = nb_x * nb_y;

    std::vector<CGAL::Bbox_3> bboxes;
    bboxes.reserve(nb);
    for (std::size_t x = 0; x < nb_x; ++ x)
      for (std::size_t y = 0; y < nb_y; ++ y)
        {
          bboxes.push_back
            (CGAL::Bbox_3 (bbox.xmin() + Dx * (x / float(nb_x)),
                           bbox.ymin() + Dy * (y / float(nb_y)),
                           bbox.zmin(),
                           bbox.xmin() + Dx * ((x+1) / float(nb_x)),
                           bbox.ymin() + Dy * ((y+1) / float(nb_y)),
                           bbox.zmax()));
        }
    
    std::cerr << "Number of divisions = " << nb_x * nb_y << std::endl;
    std::cerr << " -> Size of division: " << Dx / nb_x << " " << Dy / nb_y << std::endl;

    std::vector<std::vector<std::size_t> > indices (nb);
    std::vector<std::pair<std::size_t, std::size_t> > input_to_indices(input.size());
    
    for (std::size_t s = 0; s < input.size(); ++ s)
      {
        CGAL::Bbox_3 b = get(bbox_map, *(input.begin() + s)).bbox();
        
        for (std::size_t i = 0; i < bboxes.size(); ++ i)
          if (CGAL::do_overlap (b, bboxes[i]))
            {
              input_to_indices[s] = std::make_pair (i, indices[i].size());
              indices[i].push_back (s);
              break;
            }
      }

    output.resize (input.size());

    internal::Classifier_graphcut<ItemRange, ItemMap, ClassificationPredicate, NeighborQuery>
      f (input, item_map, labels, predicate, neighbor_query, weight, indices, input_to_indices, output);
    
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
