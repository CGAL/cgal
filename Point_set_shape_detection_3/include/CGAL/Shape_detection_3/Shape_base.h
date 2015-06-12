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

      std::vector<std::pair<FT, FT> > parameterSpace;
      parameterSpace.resize(indices.size());

      parameters(m_indices, parameterSpace, min, max);
      int iMin[2], iMax[2];
      iMin[0] = (int) (min[0] / cluster_epsilon);
      iMin[1] = (int) (min[1] / cluster_epsilon);
      iMax[0] = (int) (max[0] / cluster_epsilon);
      iMax[1] = (int) (max[1] / cluster_epsilon);

      std::size_t uExtent = abs(iMax[0] - iMin[0]) + 2;
      std::size_t vExtent = abs(iMax[1] - iMin[1]) + 2;

      std::vector<std::vector<int> > bitmap;
      std::vector<bool> visited;
      bitmap.resize(uExtent * vExtent);
      visited.resize(uExtent * vExtent, false);

      bool wrapU = wraps_u();
      bool wrapV = wraps_v();

      for (std::size_t i = 0;i<parameterSpace.size();i++) {
        int u = (int)((parameterSpace[i].first - min[0]) / cluster_epsilon);
        int v = (int)((parameterSpace[i].second - min[1]) / cluster_epsilon);
        if (u < 0 || (std::size_t)u >= uExtent) {
          if (wrapU) {
            while (u < 0) u += uExtent;
            while ((std::size_t)u >= uExtent) u-= uExtent;
          }
          else {
            u = (u < 0) ? 0 : ((std::size_t)u >= uExtent) ? (int)uExtent - 1 : u;
          }
        }
        if (v < 0 || (std::size_t)v >= vExtent) {
          if (wrapV) {
            while (v < 0) v += vExtent;
            while ((std::size_t)v >= vExtent) v-= vExtent;
          }
          else {
            v = (v < 0) ? 0 : ((std::size_t)v >= vExtent) ? (int)vExtent - 1 : v;
          }
        }
        bitmap[v * uExtent + u].push_back(m_indices[i]);
      }

      std::vector<std::vector<std::size_t> > cluster;
      for (std::size_t i = 0;i<(uExtent * vExtent);i++) {
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
          int vIndex = f / uExtent;
          int uIndex = f % uExtent;
          bool upperBorder = vIndex == 0;
          bool lowerBorder = vIndex == ((int)vExtent - 1);
          bool leftBorder = uIndex == 0;
          bool rightBorder = uIndex == ((int)uExtent - 1);

          int n;
          if (!upperBorder) {
            n = f - uExtent;
            if (!visited[n])
              fields.push(n);
          }
          else if (wrapV) {
            n = f + (vExtent - 1) * uExtent;
            if (!visited[n]) fields.push(n);
          }

          if (!leftBorder) {
            n = f - 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f + uExtent - 1;
            if (!visited[n]) fields.push(n);
          }

          if (!lowerBorder) {
            n = f + uExtent;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapV) {
            n = f - (vExtent - 1) * uExtent;
            if (!visited[n]) fields.push(n);
          }

          if (!rightBorder) {
            n = f + 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f - uExtent + 1;
            if (!visited[n]) fields.push(n);
          }
        }
      }

      int maxCluster = 0;
      for (std::size_t i = 1;i<cluster.size();i++) {
        if (cluster[i].size() > cluster[maxCluster].size()) {
          maxCluster = i;
        }
      }

      indices = cluster[maxCluster];

      return m_score = indices.size();
    }
    
    /*!
      Determines the largest cluster with a point-to-point
      distance not larger than `cluster_epsilon`. This general version performs
      a region growing within the inliers using a kd-tree.
     */
    std::size_t connected_component_kdTree(std::vector<std::size_t> &indices,
                                           FT cluster_epsilon) {
      typedef boost::tuple<Point_3,int> Point_and_int;
      typedef CGAL::Search_traits_adapter<Point_and_int,
        CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
        typename Traits::Search_traits> Search_traits_adapter;

      typedef CGAL::Kd_tree<Search_traits_adapter> Kd_Tree;
      typedef CGAL::Fuzzy_sphere<Search_traits_adapter> Fuzzy_sphere;

      m_has_connected_component = true;
      
      std::vector<Point_and_int> pts;
      std::vector<std::size_t> labelMap;
      pts.resize(indices.size());
      labelMap.resize(indices.size(), 0);

      for (std::size_t i = 0;i < indices.size();i++) {
        pts[i] = Point_and_int(point(indices[i]), i);
      }

      // construct kd tree
      Kd_Tree tree(pts.begin(), pts.end());
      
      std::stack<int> stack;
      std::size_t unlabeled = pts.size();
      std::size_t label = 1;
      std::size_t best = 0;
      std::size_t bestSize = 0;

      for (std::size_t i = 0;i<pts.size();i++) {
        if (labelMap[i] != 0)
          continue;

        std::size_t assigned = 0;

        stack.push(i);
        while(!stack.empty()) {
          std::vector<Point_and_int> nearPoints;

          std::size_t p = stack.top();
          stack.pop();

          Fuzzy_sphere fs(pts[p], cluster_epsilon, 0);
          tree.search(std::back_inserter(nearPoints), fs);

          for (std::size_t j = 0;j<nearPoints.size();j++) {
            std::size_t index = boost::get<1>(nearPoints[j]);
            if (index == p)
              continue;

            if (labelMap[index] != label) {
              labelMap[index] = label;
              assigned++;
              stack.push(index);
            }
          }
        }
        
        // Track most prominent label and remaining points
        unlabeled -= assigned;
        if (assigned > bestSize) {
          best = label;
          bestSize = assigned;
        }

        label++;

        // Can we stop already?
        if (unlabeled <= bestSize)
          break;
      }

      std::vector<std::size_t> tmpIndices;
      tmpIndices.reserve(bestSize);
      for (std::size_t i = 0;i<pts.size();i++) {
        if (labelMap[i] == best)
          tmpIndices.push_back(indices[i]);
      }

      indices = tmpIndices;
      
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

    void update_points(const std::vector<int> &shapeIndex) {
      if (!m_indices.size())
        return;
      std::size_t start = 0, end = m_indices.size() - 1;
      while (start < end) {
        while (shapeIndex[m_indices[start]] == -1
          && start < end) start++;

        while (shapeIndex[m_indices[end]] != -1
          && start < end) end--;

        if (shapeIndex[m_indices[start]] != -1
          && shapeIndex[m_indices[end]] == -1
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
                            std::vector<std::pair<FT, FT> >& parameterSpace,
                            FT min[2],
                            FT max[2]) const {
      // Avoid compiler warnings about unused parameters.
      (void)indices;
      (void)parameterSpace;
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

      std::size_t scoreBefore = m_indices.size();

      FT eps = epsilon * epsilon;
      for (std::size_t i = 0;i<indices.size();i++) {
          if (dists[i] <= eps && angles[i] > normal_threshold)
            m_indices.push_back(indices[i]);
        }

      return m_indices.size() - scoreBefore;
    }

    template<typename T> bool is_finite(T arg) {
      return arg == arg && 
        arg != std::numeric_limits<T>::infinity() &&
        arg != -std::numeric_limits<T>::infinity();
    }

    void compute_bound(const int sizeS1, const int sizeP) {
      hypergeometrical_dist(-2 - sizeS1,
                            -2 - sizeP,
                            -1 - signed(m_indices.size()),
                            m_lower_bound, m_upper_bound);

      m_lower_bound = -1 - m_lower_bound;
      m_upper_bound = -1 - m_upper_bound;
    }

    void hypergeometrical_dist(const int UN, 
                               const int x,
                               const FT n, 
                               FT &low,
                               FT &high) {                           
      const FT sq = sqrt(x * n * (UN- x) * (UN - n) / (UN - 1));
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
