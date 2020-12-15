// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_SHAPE_BASE_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_SHAPE_BASE_H

#include <CGAL/license/Shape_detection.h>

#include <vector>
#include <set>
#include <stack>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/Random.h>
#include <CGAL/Origin.h>

#ifndef CGAL_M_PI_2
#define CGAL_M_PI_2 1.57079632679489661923
#endif
#ifndef CGAL_M_PI_4
#define CGAL_M_PI_4 0.785398163397448309616
#endif

namespace CGAL {
  namespace Shape_detection {
  namespace internal {
    template<class PointAccessor>
    class Octree;
  }

    /*!
     \ingroup PkgShapeDetectionRANSACShapes
     \brief Base class for shape types that defines an interface to construct a
            shape from a set of points and to compute the point distance and normal
            deviation from the surface normal. It is used during detection to
            identify the inliers from the input data and to extract the largest
            connected component in the inlier points.
     */
  template <class Traits>
  class Shape_base {
    /// \cond SKIP_IN_MANUAL
    template <class T>
    friend class Efficient_RANSAC;
    template <class T>
    friend class Region_growing_depr;
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

    typedef typename Traits::FT FT; ///< Number type.
    typedef typename Traits::Point_3 Point_3; ///< Point type.
    typedef typename Traits::Vector_3 Vector_3; ///< Vector type.

    // \todo The property maps should be passed here instead of `compute`
    Shape_base() :
    m_is_valid(false),
      m_lower_bound((std::numeric_limits<FT>::min)()),
      m_upper_bound((std::numeric_limits<FT>::min)()),
      m_score(0),
      m_sum_expected_value(0),
      m_nb_subset_used(0) {
    }

    virtual ~Shape_base() {}

    /*!
      Returns the indices of the points in the input range assigned to this shape.
     */
    const std::vector<std::size_t> & indices_of_assigned_points() const {
      return m_indices;
    }

    /*!
      Returns a string containing the shape type
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
    virtual std::size_t connected_component(
      std::vector<std::size_t> &indices, FT cluster_epsilon) {
      if (indices.size() == 0)
        return 0;

      if (!this->supports_connected_component())
        return connected_component_kdTree(indices, cluster_epsilon);

      // Fetching parameters
      FT min[] = {0,0}, max[] = {0,0};

      std::vector<std::pair<FT, FT> > parameter_space;
      parameter_space.resize(indices.size());

      parameters(m_indices, parameter_space, cluster_epsilon, min, max);

      // Determine required size of bitmap
      std::size_t u_extent = std::size_t(ceil((max[0] - min[0]) / cluster_epsilon));
      std::size_t v_extent = std::size_t(ceil((max[1] - min[1]) / cluster_epsilon));

      // Handle singular case
      u_extent = (u_extent == 0) ? 1 : u_extent;
      v_extent = (v_extent == 0) ? 1 : v_extent;

      std::vector<unsigned int> bitmap;
      bitmap.resize(u_extent * v_extent, 0);

      // Fill bitmap
      for (std::size_t i = 0;i<parameter_space.size();i++) {
        int u = (int)((parameter_space[i].first - min[0]) / cluster_epsilon);
        int v = (int)((parameter_space[i].second - min[1]) / cluster_epsilon);

        u = (u < 0) ? 0 : (((std::size_t)u >= u_extent) ? (int)u_extent - 1 : u);
        v = (v < 0) ? 0 : (((std::size_t)v >= v_extent) ? (int)v_extent - 1 : v);

        bitmap[size_t(v) * u_extent + size_t(u)] = true;
      }

      // Iterate through the bitmap
      std::vector<unsigned int> map;
      map.reserve(64);
      map.resize(2);

      for (std::size_t y = 0;y<v_extent;y++) {
        for (std::size_t x = 0;x<u_extent;x++) {
          if (!bitmap[y * u_extent + x])
            continue;

          unsigned int w = (x > 0) ?
            bitmap[y  * u_extent + x - 1] : 0;

          unsigned int n = (y > 0) ?
            bitmap[(y - 1) * u_extent + x] : 0;

          unsigned int nw = (x > 0 && y > 0) ?
            bitmap[(y - 1) * u_extent + x - 1] : 0;

          unsigned int ne = ((x + 1 < u_extent) && y > 0) ?
            bitmap[(y - 1) * u_extent + x + 1] : 0;

          // Find smallest set label;
          unsigned int curLabel = static_cast<unsigned int>(map.size());

          curLabel = (w != 0) ?
            (std::min<unsigned int>)(curLabel, w) : curLabel;

          curLabel = (n != 0) ?
            (std::min<unsigned int>)(curLabel, n) : curLabel;

          curLabel = (nw != 0) ?
            (std::min<unsigned int>)(curLabel, nw) : curLabel;

          curLabel = (ne != 0) ?
            (std::min<unsigned int>)(curLabel, ne) : curLabel;

          // Update merge map.
          if (curLabel != map.size()) {
            if (w > curLabel) update_label(map, w, curLabel);
            if (nw > curLabel) update_label(map, nw, curLabel);
            if (n > curLabel) update_label(map, n, curLabel);
            if (ne > curLabel) update_label(map, ne, curLabel);
          }
          else map.push_back(static_cast<unsigned int>(map.size()));

          bitmap[y * u_extent + x] = curLabel;
        }
      }

      // post_wrap to handle boundaries in different shape types.
      if (map.size() > 3)
        post_wrap(bitmap, u_extent, v_extent, map);

      // Propagate label changes
      for (unsigned int j = 3; j < static_cast<unsigned int>(map.size()); j++)
        update_label(map, j, map[j]);

      // Update labels
      for (std::size_t y = 0;y<v_extent;y++)
        for (std::size_t x = 0;x<u_extent;x++) {
          unsigned int label = bitmap[y * u_extent + x];

          if (!label)
            continue;

          if (map[label] != label)
            bitmap[y * u_extent + x] = map[label];
        }

      // Count points per label.
      std::vector<unsigned int> count(map.size(), 0);

      for (std::size_t i = 0;i<parameter_space.size();i++) {
        int u = (int)((parameter_space[i].first - min[0]) / cluster_epsilon);
        int v = (int)((parameter_space[i].second - min[1]) / cluster_epsilon);

        u = (u < 0) ? 0 : (((std::size_t)u >= u_extent) ? (int)u_extent - 1 : u);
        v = (v < 0) ? 0 : (((std::size_t)v >= v_extent) ? (int)v_extent - 1 : v);

        count[bitmap[size_t(v) * u_extent + size_t(u)]]++;
      }

      // Find largest component. Start at index 2 as 0/1 are reserved for
      // basic free/occupied bitmap labels.
      unsigned int largest = 2;
      for (unsigned int i = 3; i < static_cast<unsigned int>(count.size()); i++)
        largest = (count[largest] < count[i]) ? i : largest;

      // Extract sought-after indices.
      std::vector<std::size_t> comp_indices;
      comp_indices.reserve(count[largest]);

      for (std::size_t i = 0;i<parameter_space.size();i++) {
        int u = (int)((parameter_space[i].first - min[0]) / cluster_epsilon);
        int v = (int)((parameter_space[i].second - min[1]) / cluster_epsilon);

        u = (u < 0) ? 0 : (((std::size_t)u >= u_extent) ? (int)u_extent - 1 : u);
        v = (v < 0) ? 0 : (((std::size_t)v >= v_extent) ? (int)v_extent - 1 : v);

        if (bitmap[size_t(v) * u_extent + size_t(u)] == largest)
          comp_indices.push_back(indices[i]);
      }

      indices = comp_indices;

      return m_score = indices.size();
    }

    /*!
      Determines the largest cluster with a point-to-point
      distance not larger than `cluster_epsilon`. This general version performs
      a region growing within the inliers using a Kd-tree.
     */
    std::size_t connected_component_kdTree(std::vector<std::size_t> &indices,
                                           FT cluster_epsilon) {
      typedef boost::tuple<Point_3, std::size_t> Point_and_size_t;
      typedef CGAL::Search_traits_adapter<Point_and_size_t,
        CGAL::Nth_of_tuple_property_map<0, Point_and_size_t>,
        typename Traits::Search_traits> Search_traits_adapter;

      typedef CGAL::Kd_tree<Search_traits_adapter> Kd_Tree;
      typedef CGAL::Fuzzy_sphere<Search_traits_adapter> Fuzzy_sphere;

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
                                  std::vector<FT> &distances) const = 0;

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

    virtual void post_wrap(const std::vector<unsigned int> &bitmap,
                   const std::size_t &u_extent,
                   const std::size_t &v_extent,
                   std::vector<unsigned int> &labels) const {
      // Avoid compiler warnings about unused parameters.
      (void) bitmap;
      (void) u_extent;
      (void) v_extent;
      (void) labels;
    }

    // return last computed score, or -1 if no score yet
    FT inline score() const {
      return FT(m_score);
    }

    int inline subsets() const {
      return m_nb_subset_used;
    }

    // sorting is performed by expected value
    operator FT() const {
      return expected_value();
    }

    void inline update_label(std::vector<unsigned int> &labels, unsigned int i,
                             unsigned int &new_value) const {
      if (labels[i] != i)
        update_label(labels, labels[i], new_value);

      if (new_value < labels[i])
        labels[i] = new_value;
      else
        new_value = labels[i];
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

    // Two shapes are considered the same if from 18 points randomly selected
    // from the shapes at least 12 match both shapes
    bool is_same(const Shape_base *other) const {
      if (!other)
        return false;
      if (other->m_indices.size() == 0)
        return true;

      const std::size_t num = 9;
      int score = 0;
      std::vector<std::size_t> indices(num);
      for (std::size_t i = 0;i<num;i++)
        indices[i] = m_indices[get_default_random()(
          static_cast<unsigned int>(m_indices.size()))];

      std::vector<FT> dists(num), angles(num);
      other->squared_distance(indices, dists);
      other->cos_to_normal(indices, angles);

      for (std::size_t i = 0;i<num;i++)
        if (dists[i] <= m_epsilon && angles[i] > m_normal_threshold)
          score++;

      if (score < 3)
        return false;

      for (std::size_t i = 0;i<num;i++)
        indices[i] = other->m_indices[get_default_random()(
          static_cast<unsigned int>(other->m_indices.size()))];

      this->squared_distance(indices, dists);
      this->cos_to_normal(indices, angles);

      for (std::size_t i = 0;i<num;i++)
        if (dists[i] <= m_epsilon && angles[i] > m_normal_threshold)
          score++;

      return (score >= 12);
    }

    virtual void parameters(const std::vector<std::size_t>& indices,
                            std::vector<std::pair<FT, FT> >& parameter_space,
                            FT &cluster_epsilon,
                            FT min[2],
                            FT max[2]) const {
      // Avoid compiler warnings about unused parameters.
      (void)indices;
      (void)parameter_space;
      (void)cluster_epsilon;
      (void)min;
      (void)max;
    }

    void compute(const std::vector<std::size_t>& indices,
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

      create_shape(indices);
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

    void compute_bound(const std::size_t num_evaluated_points,
                       const std::size_t num_available_points) {
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
      const FT xn = FT(double(x) * double(n));
      const FT q = FT(xn * double(UN - x) * (UN - n) / (UN - 1));
      const FT sq = CGAL::sqrt(q);
      low  = (xn - sq) / UN;
      high = (xn + sq)/UN;

      if (!is_finite<FT>(low) || !is_finite<FT>(high)) {
        low = high = 0;
      }
    }

    virtual bool supports_connected_component() const {
      return false;
    };

    // ------------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------------
    FT get_x(const Vector_3& v) const { return m_traits.compute_x_3_object()(v); }
    FT get_y(const Vector_3& v) const { return m_traits.compute_y_3_object()(v); }
    FT get_z(const Vector_3& v) const { return m_traits.compute_z_3_object()(v); }
    FT get_x(const Point_3& p) const { return m_traits.compute_x_3_object()(p); }
    FT get_y(const Point_3& p) const { return m_traits.compute_y_3_object()(p); }
    FT get_z(const Point_3& p) const { return m_traits.compute_z_3_object()(p); }

    Point_3 constr_pt() const
    { return m_traits.construct_point_3_object()(ORIGIN); }
    Point_3 constr_pt(FT x, FT y, FT z) const
    { return m_traits.construct_point_3_object()(x, y, z); }
    Vector_3 constr_vec() const
    { return m_traits.construct_vector_3_object()(NULL_VECTOR); }
    Vector_3 constr_vec(const Point_3& p, const Point_3& q) const
    { return m_traits.construct_vector_3_object()(p, q); }

    FT sqlen(const Vector_3& v) const
    { return m_traits.compute_squared_length_3_object()(v); }
    Vector_3 scale(const Vector_3& v, FT scale) const
    { return m_traits.construct_scaled_vector_3_object()(v, scale); }
    Vector_3 sum_vectors(const Vector_3& u, const Vector_3& v) const
    { return m_traits.construct_sum_of_vectors_3_object()(u, v); }
    Point_3 transl(const Point_3& p, const Vector_3 &v) const
    { return m_traits.construct_translated_point_3_object()(p, v); }
    FT scalar_pdct(const Vector_3& u, const Vector_3& v) const
    { return m_traits.compute_scalar_product_3_object()(u, v); }
    Vector_3 cross_pdct(const Vector_3& u, const Vector_3& v) const
    { return m_traits.construct_cross_product_vector_3_object()(u, v); }

  protected:
    /// \endcond
    //
    /// \cond SKIP_IN_MANUAL
    /*!
      Contains indices of the inliers of the candidate, access
      to the point and normal data is provided via property maps.
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

    Input_iterator m_first;

    Traits m_traits;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;
    /// \endcond
  };
}
}
#endif
