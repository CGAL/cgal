#ifndef CGAL_SHAPE_DETECTION_3_SHAPE_BASE_H
#define CGAL_SHAPE_DETECTION_3_SHAPE_BASE_H

#include <vector>
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
    /// \cond SKIP_IN_MANUAL
    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    /// \endcond

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

    const std::vector<size_t> &assigned_points() const {
      return m_indices;
    }
      
    /*!
      Helper function writing shape type and numerical parameters into a string.
     */

    virtual std::string info() const {
      return std::string();
    }

    /*!
      Provides the squared Euclidean distance of the point to the shape.
     */

    virtual FT squared_distance(const Point &p) const = 0;

  protected:
      
    /*!
      Constructs the shape based on a minimal set of samples from the input data.
     */
    virtual void create_shape(const std::vector<size_t> &indices) = 0;
    
    /*!
      Determines the largest cluster of with a point to point distance not larger than cluster_epsilon.
     */
    size_t connected_component(FT cluster_epsilon) {
      if (m_indices.size() == 0)
        return 0;

      //       if (m_hasconnected_component)
      //         return m_score;

      m_has_connected_component = true;
      if (!supports_connected_component())
        return m_indices.size();

      FT min[2], max[2];
      //parameterExtend(center, width, min, max);
      std::vector<std::pair<FT, FT> > parameterSpace;
      parameterSpace.resize(m_indices.size());

      parameters(parameterSpace, m_indices, min, max);
      int iMin[2], iMax[2];
      iMin[0] = (int) (min[0] / cluster_epsilon);
      iMin[1] = (int) (min[1] / cluster_epsilon);
      iMax[0] = (int) (max[0] / cluster_epsilon);
      iMax[1] = (int) (max[1] / cluster_epsilon);

      size_t uExtent = abs(iMax[0] - iMin[0]) + 2;
      size_t vExtent = abs(iMax[1] - iMin[1]) + 2;

      std::vector<std::vector<int> > bitmap;
      std::vector<bool> visited;
      bitmap.resize(uExtent * vExtent);
      visited.resize(uExtent * vExtent, false);

      bool wrapU = wraps_u();
      bool wrapV = wraps_v();

      for (size_t i = 0;i<parameterSpace.size();i++) {
        int u = (parameterSpace[i].first - min[0]) / cluster_epsilon;
        int v = (parameterSpace[i].second - min[1]) / cluster_epsilon;
        if (u < 0 || u >= uExtent) {
          if (wrapU) {
            while (u < 0) u += uExtent;
            while (u >= uExtent) u-= uExtent;
          }
          else {
            std::cout << "cc: u out of bounds: " << u << std::endl;
            u = (u < 0) ? 0 : (u >= uExtent) ? uExtent - 1 : u;
          }
        }
        if (v < 0 || v >= vExtent) {
          if (wrapV) {
            while (v < 0) v += vExtent;
            while (v >= vExtent) v-= vExtent;
          }
          else {
            std::cout << "cc: v out of bounds: " << u << std::endl;
            v = (v < 0) ? 0 : (v >= vExtent) ? vExtent - 1 : v;
          }
        }
        bitmap[v * uExtent + u].push_back(m_indices[i]);
      }

      std::vector<std::vector<size_t> > cluster;
      for (size_t i = 0;i<(uExtent * vExtent);i++) {
        cluster.push_back(std::vector<size_t>());
        if (bitmap[i].empty())
          continue;
        if (visited[i])
          continue;

        std::stack<size_t> fields;
        fields.push(i);
        while (!fields.empty()) {
          size_t f = fields.top();
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
          bool lowerBorder = vIndex == (vExtent - 1);
          bool leftBorder = uIndex == 0;
          bool rightBorder = uIndex == (uExtent - 1);

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
      for (size_t i = 1;i<cluster.size();i++) {
        if (cluster[i].size() > cluster[maxCluster].size()) {
          maxCluster = i;
        }
      }

      m_indices = cluster[maxCluster];

      return m_score = m_indices.size();
    }
 
    /*!
      Provides the squared Euclidean distance of a set of points.
     */
    virtual void squared_distance(std::vector<FT> &dists,
      const std::vector<size_t> &indices) = 0;   

    /*!
      Provides the deviation of the point normal from the surface normal at the projected point in form of the dot product.
     */
    virtual void cos_to_normal(std::vector<FT> &angles,
      const std::vector<size_t> &indices) const = 0;

    /*!
      Defines the minimal number of samples required for construction.
     */
    virtual size_t required_samples() const = 0;

    /*!
      Retrieves the point for an index.
     */
    const Point_3 &get_point(size_t i) const {
      return get(this->m_point_pmap, *(this->m_first + i));
    }
    
    /*!
      Retrieves the normal vector for an index.
     */
    const Vector_3 &get_normal(size_t i) const {
      return get(this->m_normal_pmap, *(this->m_first + i));
    }
    
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
      size_t start = 0, end = m_indices.size() - 1;
      while (start < end) {
        while (shapeIndex[m_indices[start]] == -1
          && start < end) start++;

        while (shapeIndex[m_indices[end]] != -1
          && start < end) end--;

        if (shapeIndex[m_indices[start]] != -1
          && shapeIndex[m_indices[end]] == -1
          && start < end) {
          size_t tmp = m_indices[start];
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

    virtual void squared_distance(std::vector<FT> &dists,
                                  const std::vector<int> &shapeIndex,
                                  const std::vector<size_t> &indices) = 0;

    //virtual FT cos_to_normal(const Point &p, const Vector &n) const = 0;

    virtual void cos_to_normal(std::vector<FT> &angles,
                               const std::vector<int> &shapeIndex,
                               const std::vector<size_t> &indices) const = 0;

    virtual void parameters(std::vector<std::pair<FT, FT> > &parameterSpace,
                            const std::vector<size_t> &indices,
                            FT min[2],
                            FT max[2]) const {
    }

    void compute(const std::set<size_t> &indices,
                 Input_iterator first,
                 Point_pmap point_pmap,
                 Normal_pmap normal_pmap,
                 FT epsilon,
                 FT normal_threshold) {
      if (indices.size() < required_samples())
        return;

      m_first = first;
      m_point_pmap = point_pmap;
      m_normal_pmap = normal_pmap;
      m_epsilon = epsilon;
      m_normal_threshold = normal_threshold;

      std::vector<size_t> output(indices.begin(), indices.end());

      create_shape(output);
    }

    inline bool operator<(const Shape &c) const {
      return expected_value() < c.expected_value();
    }

    size_t cost_function(const std::vector<int> &shapeIndex,
                         FT epsilon,
                         FT normal_threshold,
                         const std::vector<size_t> &indices) {
      std::vector<FT> dists, angles;
      dists.resize(indices.size());
      squared_distance(dists, shapeIndex, indices);
      angles.resize(indices.size());
      cos_to_normal(angles, shapeIndex, indices);

      size_t scoreBefore = m_indices.size();

      FT eps = epsilon * epsilon;
      for (size_t i = 0;i<indices.size();i++) {
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
      if (UN == 1 || UN == 0)
        printf("something wrong here, \
               denominator is zero (UN %d)!! \n", UN);
      if (x > UN)
        printf("SizeP1 smaller than sizeP, something wrong \
               (and sqrt may be negative !!!!");
      FT sq = sqrt(x * n * (UN- x) * (UN - n) / (UN - 1));
      low  = (x * n - sq) / UN;
      high = (x * n + sq)/UN;

      if (!is_finite<FT>(low) || !is_finite<FT>(high)) {
        low = high = 0;
      }
    }


    virtual bool supports_connected_component() {
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
    /*!
      Contains indices of the points supporting the candidate. Access to the point and normal data is provided via property maps.
     */
    std::vector<size_t> m_indices;
    /// \cond SKIP_IN_MANUAL

    FT m_epsilon;

    //deviation of normal, used during first check of the 3 normal
    FT m_normal_threshold;

    bool m_isValid;
    FT m_lower_bound;
    FT m_upper_bound;

    size_t m_score;

    FT m_sum_expected_value;

    //count the number of subset used so far for the score,
    //and thus indicate the next one to use
    size_t m_nb_subset_used;
    bool m_has_connected_component;

    //indices of the points fitting to the candidate
    std::vector<size_t> m_indices;
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
