// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_3_SHAPE_BASE_H
#define CGAL_SHAPE_DETECTION_3_SHAPE_BASE_H

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

/*!
 \file Shape_base.h
 */

// CODE REVIEW
// make code more modular: connected_component()
// use const where relevant, eg wrapU
// initialize all variables including max


namespace CGAL {
  namespace Shape_detection_3 {
  namespace internal {
    template<class PointAccessor>
    class Octree;
  }
    
    /*!
     \ingroup PkgPointSetShapeDetection3Shapes
     \brief Base class for shape types defining an interface to construct a
            shape from a set of points and to compute the point distance and normal
            deviation from the surface normal. It is used during detection to
            identify the inliers from the input data and to extract the largest
            connected component in inlier points.
     */
  template <class Traits>
  class Shape_base {
    /// \cond SKIP_IN_MANUAL
    template <class T>
    friend class Efficient_RANSAC;
    template<class PointAccessor>
    friend class internal::Octree;
    /// \endcond

  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Input_range::iterator Input_iterator;
      ///< random access iterator for input data.
    typedef typename Traits::Point_map Point_map;
      ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
      ///< property map to access the unoriented normal of an input point.
    typedef Shape_base<Traits> Shape;
      ///< own type.
    /// \endcond

    typedef typename Traits::FT FT; ///< number type.
    typedef typename Traits::Point_3 Point_3; ///< point type.
    typedef typename Traits::Vector_3 Vector_3; ///< vector type.

    Shape_base() :
    m_is_valid(false),
      m_lower_bound((std::numeric_limits<FT>::min)()),
      m_upper_bound((std::numeric_limits<FT>::min)()),
      m_score(0),
      m_sum_expected_value(0),
      m_nb_subset_used(0),
      m_has_connected_component(false) {
    }

    virtual ~Shape_base() {}
      
    /*!
      returns the indices of the points in the input range assigned to this shape.
     */
    const std::vector<std::size_t> & indices_of_assigned_points() const {
      return m_indices;
    }
      
    /*!
      returns a string containing the shape type
      and the numerical parameters.
     */
    virtual std::string info() const {
      return std::string();
    }

    /*!
      Computes the squared Euclidean distance from the query point `p` to the shape.
     */
    virtual FT squared_distance(const Point_3 &p) const = 0;

  protected:
      
    /*!
      Constructs the shape based on a minimal set of samples from the
      input data.
     */
    virtual void create_shape(const std::vector<std::size_t> &indices) = 0;
    
    /*!
      Determines the largest cluster of inlier points. A point belongs to a cluster
      if there is a point in the cluster closer than `cluster_epsilon` distance.
     */
    std::size_t connected_component(std::vector<std::size_t> &indices, FT cluster_epsilon) {
      if (indices.size() == 0)
        return 0;

      if (m_has_connected_component)
        return m_score;

      m_has_connected_component = true;
      if (!this->supports_connected_component())
        return connected_component_kdTree(indices, cluster_epsilon);

      FT min[] = {0,0}, max[] = {0,0};

      std::vector<std::pair<FT, FT> > parameter_space;
      parameter_space.resize(indices.size());

      parameters(m_indices, parameter_space, min, max);
      int i_min[2], i_max[2];
      i_min[0] = (int) (min[0] / cluster_epsilon);
      i_min[1] = (int) (min[1] / cluster_epsilon);
      i_max[0] = (int) (max[0] / cluster_epsilon);
      i_max[1] = (int) (max[1] / cluster_epsilon);

      std::size_t u_extent = CGAL::abs(i_max[0] - i_min[0]) + 1;
      std::size_t v_extent = CGAL::abs(i_max[1] - i_min[1]) + 1;

      std::vector<std::vector<std::size_t> > bitmap;
      std::vector<bool> visited;
      bitmap.resize(u_extent * v_extent);
      visited.resize(u_extent * v_extent, false);

      bool wrap_u = wraps_u();
      bool wrap_v = wraps_v();

      for (std::size_t i = 0;i<parameter_space.size();i++) {
        int u = (int)((parameter_space[i].first - min[0]) / cluster_epsilon);
        int v = (int)((parameter_space[i].second - min[1]) / cluster_epsilon);

        if (u < 0 || (std::size_t)u >= u_extent) {
          if (wrap_u) {
            while (u < 0) u += (int) u_extent;
            while (u >= (int) u_extent) u-= (int)u_extent;
          }
          else {
            u = (u < 0) ? 0 : (u >= (int) u_extent) ? (int)u_extent - 1 : u;
          }
        }
        if (v < 0 || v >= (int) v_extent) {
          if (wrap_v) {
            while (v < 0) v += (int) v_extent;
            while (v >= (int) v_extent) v-= (int) v_extent;
          }
          else {
            v = (v < 0) ? 0 : (v >= (int) v_extent) ? (int) v_extent - 1 : v;
          }
        }

        bitmap[v * int(u_extent) + u].push_back(m_indices[i]);
      }

      std::vector<std::vector<std::size_t> > cluster;
      for (std::size_t i = 0;i<(u_extent * v_extent);i++) {
        cluster.push_back(std::vector<std::size_t>());
        if (bitmap[i].empty())
          continue;
        if (visited[i])
          continue;

        std::stack<std::size_t> fields;
        fields.push(i);
        while (!fields.empty()) {
          std::size_t f = fields.top();
          fields.pop();
          if (visited[f])
            continue;
          visited[f] = true;
          if (bitmap[f].empty())
            continue;

          // copy indices
          std::copy(bitmap[f].begin(), bitmap[f].end(),
            std::back_inserter(cluster.back()));

          // grow 8-neighborhood
          int v_index = int(f / u_extent);
          int u_index = int(f % u_extent);
          bool upper_border = v_index == 0;
          bool lower_border = v_index == ((int)v_extent - 1);
          bool left_border = u_index == 0;
          bool right_border = u_index == ((int)u_extent - 1);

          int n;
          if (!upper_border) {
            n = int(f - u_extent);
            if (!visited[n])
              fields.push(n);
          }
          else if (wrap_v) {
            n = int((f + v_extent - 1) * u_extent);
            if (!visited[n]) fields.push(n);
          }

          if (!left_border) {
            n = int(f - 1);
            if (!visited[n]) fields.push(n);
          }
          else if (wrap_u) {
            n = int(f + u_extent - 1);
            if (!visited[n]) fields.push(n);
          }

          if (!lower_border) {
            n = int(f + u_extent);
            if (!visited[n]) fields.push(n);
          }
          else if (wrap_v) {
            n = int((f - (v_extent - 1)) * u_extent);
            if (!visited[n]) fields.push(n);
          }

          if (!right_border) {
            n = int(f) + 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrap_u) {
            n = int(f - u_extent + 1);
            if (!visited[n]) fields.push(n);
          }
        }
      }

      std::size_t max_cluster = 0;
      for (std::size_t i = 1;i<cluster.size();i++) {
        if (cluster[i].size() > cluster[max_cluster].size()) {
          max_cluster = i;
        }
      }

      indices = cluster[max_cluster];

      return m_score = indices.size();
    }
    
    /*!
      Determines the largest cluster with a point-to-point
      distance not larger than `cluster_epsilon`. This general version performs
      a region growing within the inliers using a kd-tree.
     */
    std::size_t connected_component_kdTree(std::vector<std::size_t> &indices,
                                           FT cluster_epsilon) {
      typedef boost::tuple<Point_3, std::size_t> Point_and_size_t;
      typedef CGAL::Search_traits_adapter<Point_and_size_t,
        CGAL::Nth_of_tuple_property_map<0, Point_and_size_t>,
        typename Traits::Search_traits> Search_traits_adapter;

      typedef CGAL::Kd_tree<Search_traits_adapter> Kd_Tree;
      typedef CGAL::Fuzzy_sphere<Search_traits_adapter> Fuzzy_sphere;

      m_has_connected_component = true;
      
      std::vector<Point_and_size_t> pts;
      std::vector<std::size_t> label_map;
      pts.resize(indices.size());
      label_map.resize(indices.size(), 0);

      for (std::size_t i = 0;i < indices.size();i++) {
        pts[i] = Point_and_size_t(point(indices[i]), i);
      }

      // construct kd tree
      Kd_Tree tree(pts.begin(), pts.end());
      
      std::stack<std::size_t> stack;
      std::size_t unlabeled = pts.size();
      std::size_t label = 1;
      std::size_t best = 0;
      std::size_t best_size = 0;

      for (std::size_t i = 0;i<pts.size();i++) {
        if (label_map[i] != 0)
          continue;

        std::size_t assigned = 0;

        stack.push(i);
        while(!stack.empty()) {
          std::vector<Point_and_size_t> near_points;

          std::size_t p = stack.top();
          stack.pop();

          Fuzzy_sphere fs(pts[p], cluster_epsilon, 0);
          tree.search(std::back_inserter(near_points), fs);

          for (std::size_t j = 0;j<near_points.size();j++) {
            std::size_t index = boost::get<1>(near_points[j]);
            if (index == p)
              continue;

            if (label_map[index] != label) {
              label_map[index] = label;
              assigned++;
              stack.push(index);
            }
          }
        }
        
        // Track most prominent label and remaining points
        unlabeled -= assigned;
        if (assigned > best_size) {
          best = label;
          best_size = assigned;
        }

        label++;

        // Can we stop already?
        if (unlabeled <= best_size)
          break;
      }

      std::vector<std::size_t> tmp_indices;
      tmp_indices.reserve(best_size);
      for (std::size_t i = 0;i<pts.size();i++) {
        if (label_map[i] == best)
          tmp_indices.push_back(indices[i]);
      }

      indices = tmp_indices;
      
      return indices.size();
    }

    /*!
      Computes the squared Euclidean distance from a set of points to the shape.
      The distances will be stored in the so called parameter.
     */
    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &distances) = 0;

    /*!
      Computes the deviation of the point normal from the surface normal at the
      projected point in form of the dot product and writes the result into the
      provided `angles` vector.
     */
    virtual void cos_to_normal(const std::vector<std::size_t> &indices,
                               std::vector<FT> &angles) const = 0;

    /*!
      Returns minimal number of sample points required for construction.
     */
    virtual std::size_t minimum_sample_size() const = 0;

    /*!
      Retrieves the point location from its index.
     */
    typename boost::property_traits< typename Traits::Point_map >::reference
    point(std::size_t i) const {
      return get(this->m_point_pmap, *(this->m_first + i));
    }
    
    /*!
      Retrieves the normal vector from its index.
     */
    typename boost::property_traits< typename Traits::Normal_map >::reference
    normal(std::size_t i) const {
      return get(this->m_normal_pmap, *(this->m_first + i));
    }

    /*!
      Retrieves the traits class.
     */
    const Traits&
    traits() const
    {
      return m_traits;
    }

    /// \cond SKIP_IN_MANUAL
    struct Compare_by_max_bound {
        bool operator() (const Shape *a, const Shape *b) {
            return a->max_bound() < b->max_bound();
        }
    };
      
    FT expected_value() const {
      return (m_lower_bound + m_upper_bound) / 2.f;
    }

    FT inline min_bound() const {
      return  m_lower_bound;
    }

    FT inline max_bound() const {
      return  m_upper_bound;
    }

    // return last computed score, or -1 if no score yet
    FT inline score() const {
      return m_score;
    } 

    int inline subsets() const {
      return m_nb_subset_used;
    }

    // sorting is performed by expected value
    operator FT() const {
      return expected_value();
    }

    void update_points(const std::vector<int> &shape_index) {
      if (!m_indices.size())
        return;
      std::size_t start = 0, end = m_indices.size() - 1;
      while (start < end) {
        while (shape_index[m_indices[start]] == -1
          && start < end) start++;

        while (shape_index[m_indices[end]] != -1
          && start < end) end--;

        if (shape_index[m_indices[start]] != -1
          && shape_index[m_indices[end]] == -1
          && start < end) {
          std::size_t tmp = m_indices[start];
          m_indices[start] = m_indices[end];
          m_indices[end] = tmp;
        }
      }
      m_indices.resize(end);
      m_score = m_indices.size();
    }

    bool is_valid() const {
      return m_is_valid;
    }

    virtual void parameters(const std::vector<std::size_t>& indices,
                            std::vector<std::pair<FT, FT> >& parameter_space,
                            FT min[2],
                            FT max[2]) const {
      // Avoid compiler warnings about unused parameters.
      (void)indices;
      (void)parameter_space;
      (void)min;
      (void)max;
    }

    void compute(const std::set<std::size_t>& indices,
                 Input_iterator first,
                 Traits traits,
                 Point_map point_pmap,
                 Normal_map normal_pmap,
                 FT epsilon,
                 FT normal_threshold) {
      if (indices.size() < minimum_sample_size())
        return;

      m_first = first;
      m_traits = traits;
      m_point_pmap = point_pmap;
      m_normal_pmap = normal_pmap;
      m_epsilon = epsilon;
      m_normal_threshold = normal_threshold;

      std::vector<std::size_t> output(indices.begin(), indices.end());

      create_shape(output);
    }

    inline bool operator<(const Shape &c) const {
      return expected_value() < c.expected_value();
    }

    std::size_t cost_function(FT epsilon,
                         FT normal_threshold,
                         const std::vector<std::size_t> &indices) {
      std::vector<FT> dists, angles;
      dists.resize(indices.size());
      squared_distance(indices, dists);
      angles.resize(indices.size());
      cos_to_normal(indices, angles);

      std::size_t score_before = m_indices.size();

      FT eps = epsilon * epsilon;
      for (std::size_t i = 0;i<indices.size();i++) {
          if (dists[i] <= eps && angles[i] > normal_threshold)
            m_indices.push_back(indices[i]);
        }

      return m_indices.size() - score_before;
    }

    template<typename T> bool is_finite(T arg) {
      return arg == arg && 
        arg != std::numeric_limits<T>::infinity() &&
        arg != -std::numeric_limits<T>::infinity();
    }

    void compute_bound(const std::size_t num_evaluated_points, const std::size_t num_available_points) {
      hypergeometrical_dist(-2 - num_evaluated_points,
                            -2 - num_available_points,
                            -1 - signed(m_indices.size()),
                            m_lower_bound, m_upper_bound);

      m_lower_bound = -1 - m_lower_bound;
      m_upper_bound = -1 - m_upper_bound;
    }

    void hypergeometrical_dist(const std::ptrdiff_t UN, 
                               const std::ptrdiff_t x,
                               const std::ptrdiff_t n, 
                               FT &low,
                               FT &high) {                           
      const FT q = FT(x * n * double(UN - x) * (UN - n) / (UN - 1));
      const FT sq = CGAL::sqrt(q);
      low  = (x * n - sq) / UN;
      high = (x * n + sq)/UN;

      if (!is_finite<FT>(low) || !is_finite<FT>(high)) {
        low = high = 0;
      }
    }


    virtual bool supports_connected_component() const {
      return false;
    };

    virtual bool wraps_u() const {
      return false;
    };

    virtual bool wraps_v() const {
      return false;
    };

  protected:
    /// \endcond
    // 
    /// \cond SKIP_IN_MANUAL
    /*!
      Contains indices of the inliers of the candidate, access to the point and normal data is provided via property maps.
     */
    std::vector<std::size_t> m_indices;

    FT m_epsilon;

    //deviation of normal, used during first check of the 3 normal
    FT m_normal_threshold;

    bool m_is_valid;
    FT m_lower_bound;
    FT m_upper_bound;

    std::size_t m_score;

    FT m_sum_expected_value;

    //count the number of subset used so far for the score,
    //and thus indicate the next one to use
    std::size_t m_nb_subset_used;
    bool m_has_connected_component;

    Input_iterator m_first;

    Traits m_traits;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;
    /// \endcond
  };  
}
}
#endif
