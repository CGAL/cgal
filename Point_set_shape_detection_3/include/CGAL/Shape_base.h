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

namespace CGAL {
  namespace internal {
    template<class PointAccessor>
    class Octree;
  }
    
    /*!
     \brief Base class of shape types. Provides access to assigned points.
     */
  template <class Sd_traits>
  class Shape_base {
      /// \cond SKIP_IN_MANUAL
    template <class T>
    friend class Shape_detection_3;
    template<class PointAccessor>
    friend class internal::Octree;
      /// \endcond

  public:
    /// \name Types 
    /// @{

    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point; ///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector; ///< vector type.
    typedef typename Sd_traits::Point_pmap Point_pmap;  ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap; ///< property map to access the unoriented normal of an input point.

    typedef Shape_base<Sd_traits> Shape; ///< own type.

    Shape_base() :
    m_isValid(true),
      m_lower_bound((std::numeric_limits<FT>::min)()),
      m_upper_bound((std::numeric_limits<FT>::min)()),
      m_score(0),
      m_sum_expected_value(0),
      m_nb_subset_used(0),
      m_has_connected_component(false) {
    }

    virtual ~Shape_base() {}
      
    /*!
      Indices into the input data of all points assigned to this shape.
    */

    const std::vector<int> &assigned_points() {
      return m_indices;
    }
      
    /*!
      Helper function writing shape type and numerical parameters into a string.
     */
    
    virtual std::string info() const = 0;

  protected:
      /// \cond SKIP_IN_MANUAL
      struct Compare_by_max_bound {
          bool operator() (Shape *a, Shape *b) {
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

    //return last computed score, or -1 if no score yet
    FT inline score() const {
      return m_score;
    } 

    int inline subsets() const {
      return m_nb_subset_used;
    }

    //so we can sort by expected value
    operator FT() const {
      return expected_value();
    }

    void update_points(const std::vector<int> &shapeIndex) {
      if (!m_indices.size())
        return;
      int start = 0, end = m_indices.size() - 1;
      while (start < end) {
        while (shapeIndex[m_indices[start]] == -1 && start < end) start++;
        while (shapeIndex[m_indices[end]] != -1 && start < end) end--;
        if (shapeIndex[m_indices[start]] != -1 && shapeIndex[m_indices[end]] == -1 && start < end) {
          unsigned int tmp = m_indices[start];
          m_indices[start] = m_indices[end];
          m_indices[end] = tmp;
        }
      }
      m_indices.resize(end);
      m_score = m_indices.size();
    }

    bool is_valid() const {
      return m_isValid;
    }

    virtual FT squared_distance(const Point &p) const = 0;
    virtual void squared_distance(std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) = 0;
    virtual FT cos_to_normal(const Point &p, const Vector &n) const = 0;
    virtual void cos_to_normal(std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const = 0;

    virtual void parameter_extend(const Point &center, FT width, FT min[2], FT max[2]) const = 0;
    virtual void parameters(std::vector<std::pair<FT, FT> > &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const = 0;

    unsigned int connected_component(FT m_bitmapEpsilon, const Point &center, FT width) {
      if (m_indices.size() == 0)
        return 0;

      //       if (m_hasconnected_component)
      //         return m_score;

      m_has_connected_component = true;
      if (!supports_connected_component())
        return m_indices.size();

      //ccCount++;
      //clock_t s, e;
      //s = clock();

      FT min[2], max[2];
      //parameterExtend(center, width, min, max);
      std::vector<std::pair<FT, FT> > parameterSpace;
      parameterSpace.resize(m_indices.size());

      parameters(parameterSpace, m_indices, min, max);
      int iMin[2], iMax[2];
      iMin[0] = min[0] / m_bitmapEpsilon;
      iMin[1] = min[1] / m_bitmapEpsilon;
      iMax[0] = max[0] / m_bitmapEpsilon;
      iMax[1] = max[1] / m_bitmapEpsilon;

      int uExtend = abs(iMax[0] - iMin[0]) + 2;
      int vExtend = abs(iMax[1] - iMin[1]) + 2;

      std::vector<std::vector<int> > bitmap;
      std::vector<bool> visited;
      bitmap.resize(uExtend * vExtend);
      visited.resize(uExtend * vExtend, false);

      bool wrapU = wraps_u();
      bool wrapV = wraps_v();

      for (unsigned int i = 0;i<parameterSpace.size();i++) {
        int u = (parameterSpace[i].first - min[0]) / m_bitmapEpsilon;
        int v = (parameterSpace[i].second - min[1]) / m_bitmapEpsilon;
        if (u < 0 || u >= uExtend) {
          if (wrapU) {
            while (u < 0) u += uExtend;
            while (u >= uExtend) u-= uExtend;
          }
          else {
            std::cout << "cc: u out of bounds: " << u << std::endl;
            u = (u < 0) ? 0 : (u >= uExtend) ? uExtend - 1 : u;
          }
        }
        if (v < 0 || v >= vExtend) {
          if (wrapV) {
            while (v < 0) v += vExtend;
            while (v >= vExtend) v-= vExtend;
          }
          else {
            std::cout << "cc: v out of bounds: " << u << std::endl;
            v = (v < 0) ? 0 : (v >= vExtend) ? vExtend - 1 : v;
          }
        }
        bitmap[v * uExtend + u].push_back(m_indices[i]);
      }

      std::vector<std::vector<int> > cluster;
      for (unsigned int i = 0;i<(uExtend * vExtend);i++) {
        cluster.push_back(std::vector<int>());
        if (bitmap[i].empty())
          continue;
        if (visited[i])
          continue;

        std::stack<int> fields;
        fields.push(i);
        while (!fields.empty()) {
          int f = fields.top();
          fields.pop();
          if (visited[f])
            continue;
          visited[f] = true;
          if (bitmap[f].empty())
            continue;

          // copy indices
          std::copy(bitmap[f].begin(), bitmap[f].end(), std::back_inserter(cluster.back()));

          // grow 8-neighborhood
          int vIndex = f / uExtend;
          int uIndex = f % uExtend;
          bool upperBorder = vIndex == 0;
          bool lowerBorder = vIndex == (vExtend - 1);
          bool leftBorder = uIndex == 0;
          bool rightBorder = uIndex == (uExtend - 1);

          int n;
          if (!upperBorder) {
            n = f - uExtend;
            if (!visited[n])
              fields.push(n);
          }
          else if (wrapV) {
            n = f + (vExtend - 1) * uExtend;
            if (!visited[n]) fields.push(n);
          }

          if (!leftBorder) {
            n = f - 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f + uExtend - 1;
            if (!visited[n]) fields.push(n);
          }

          if (!lowerBorder) {
            n = f + uExtend;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapV) {
            n = f - (vExtend - 1) * uExtend;
            if (!visited[n]) fields.push(n);
          }

          if (!rightBorder) {
            n = f + 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f - uExtend + 1;
            if (!visited[n]) fields.push(n);
          }
        }
      }

      int maxCluster = 0;
      for (unsigned int i = 1;i<cluster.size();i++) {
        if (cluster[i].size() > cluster[maxCluster].size()) {
          maxCluster = i;
        }
      }

      m_indices = cluster[maxCluster];

      //e = clock();
      //ccTime += e - s;

      return m_score = m_indices.size();
    }

    void compute(const std::set<int> &indices, Input_iterator first, Point_pmap point_pmap, Normal_pmap normal_pmap, FT epsilon, FT normal_threshold) {
      if (indices.size() < required_samples())
        return;

      m_first = first;
      m_point_pmap = point_pmap;
      m_normal_pmap = normal_pmap;
      m_epsilon = epsilon;
      m_normal_threshold = normal_threshold;

      std::vector<int> output(indices.begin(), indices.end());

      create_shape(output);
    }

    virtual void create_shape(const std::vector<int> &indices) = 0;

    inline bool operator<(const Shape &c) const {
      return expected_value() < c.expected_value();
    }

    unsigned int cost_function(const std::vector<int> &shapeIndex, FT epsilon, FT normal_threshold, const std::vector<unsigned int> &indices) {
      std::vector<FT> dists, angles;
      dists.resize(indices.size());
      squared_distance(dists, shapeIndex, indices);
      angles.resize(indices.size());
      cos_to_normal(angles, shapeIndex, indices);

      unsigned int scoreBefore = m_indices.size();

      FT eps = epsilon * epsilon;
      for (unsigned int i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          if (dists[i] <= eps && angles[i] > normal_threshold)
            m_indices.push_back(indices[i]);
        }
      }

      return m_indices.size() - scoreBefore;
    }

    template<typename T> bool is_finite(T arg) {
      return arg == arg && 
        arg != std::numeric_limits<T>::infinity() &&
        arg != -std::numeric_limits<T>::infinity();
    }

    void compute_bound(const int sizeS1, const int sizeP) {
      hypergeometrical_dist(-2 - sizeS1, -2 - sizeP, -1 - signed(m_indices.size()), m_lower_bound, m_upper_bound);
      m_lower_bound = -1 - m_lower_bound;
      m_upper_bound = -1 - m_upper_bound;
    }

    void hypergeometrical_dist(const int UN, const int x, const FT n, FT &low, FT &high) {
      if (UN == 1 || UN == 0)
        printf("something wrong here, denominator is zero (UN %d)!! \n", UN);
      if (x > UN)
        printf("SizeP1 smaller than sizeP, something wrong (and sqrt may be negative !!!!");
      FT sq = sqrtf(x * n * (UN- x) * (UN - n) / (UN - 1));
      low  = (x * n - sq) / UN;
      high = (x * n + sq)/UN;

      if (!is_finite<FT>(low) || !is_finite<FT>(high)) {
        low = high = 0;
      }
    }

    virtual int required_samples() const = 0;

    virtual bool supports_connected_component() const = 0;
    virtual bool wraps_u() const = 0;
    virtual bool wraps_v() const = 0;

  protected:
    FT m_epsilon;
    FT m_normal_threshold;	 //deviation of normal, used during first check of the 3 normal

    bool m_isValid;
    FT m_lower_bound;
    FT m_upper_bound;

    unsigned int m_score;

    FT m_sum_expected_value;
    int m_nb_subset_used;		//count the number of subset used so far for the score, and thus indicate the next one to use
    bool m_has_connected_component;

    std::vector<int> m_indices;	//indices of the points fitting to the candidate
    Input_iterator m_first;
    Point_pmap m_point_pmap;
      Normal_pmap m_normal_pmap;
      /// \endcond
  };

  namespace internal {
    class Shape_factory_base {
    public:
    virtual ~Shape_factory_base() {}
      virtual void *create() = 0;
    };
  }
    
    /*!
     \brief Template class for creating a factory for a shape type.
     */
  template<class Shape>
  class Shape_factory : public internal::Shape_factory_base {
  public:
      /*!
       Returns a new instance of the shape type.
       */
    virtual void *create() {
      return new Shape;
    }
  };
}
#endif
