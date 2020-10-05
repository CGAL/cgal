// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Florent Lafarge, Simon Giraudot

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_DEPRECATED_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_DEPRECATED_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <functional>

// Shape detection includes.
#include <CGAL/Shape_detection/Efficient_RANSAC/Plane.h>

// CGAL includes.
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/squared_distance_3.h>

// CGAL boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>

// Boost includes.
#include <boost/make_shared.hpp>
#include <boost/iterator/filter_iterator.hpp>

// Deprecated -->
#define CGAL_DEPRECATED_HEADER "<CGAL/Shape_detection/deprecated/Region_growing.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Shape_detection/Region_growing/Region_growing.h>"

#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "CGAL::Shape_detection_3::Region_growing<> has been replaced by the class "\
  "CGAL::Shape_detection::Region_growing<>."

#include <CGAL/internal/deprecation_warning.h>

namespace CGAL {
namespace Shape_detection {

  #ifdef DOXYGEN_NS
    namespace deprecated {
  #endif

  /*!
  \ingroup PkgShapeDetectionDEPR
  \brief A shape detection algorithm using a region growing method.

  Given a point set in 3D space with unoriented normals, sampled on surfaces,
  this class enables to detect subsets of connected points lying on the surface of primitive shapes.
  Each input point is assigned to either none or at most one detected primitive
  shape. The implementation follows \cgalCite{cgal:lm-clscm-12}.

  \tparam Traits a model of `EfficientRANSACTraits`
  */
  template<class Traits>
  class Region_growing_depr {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Input_range::iterator Input_iterator;
    typedef typename Traits::FT FT; ///< number type.
    typedef typename Traits::Point_3 Point; ///< point type.
    typedef typename Traits::Vector_3 Vector; ///< vector type.
    typedef typename Traits::Plane_3 Plane; ///< plane type.
    /// \endcond

    typedef typename Traits::Input_range Input_range;
    ///< Model of the concept `Range` with random access iterators, providing input points and normals
    /// through the following two property maps.

    typedef typename Traits::Point_map Point_map;
    ///< property map to access the location of an input point.

    typedef typename Traits::Normal_map Normal_map;
    ///< property map to access the unoriented normal of an input point.

    typedef CGAL::Shape_detection::Shape_base<Traits> Shape; ///< shape type.

    typedef CGAL::Shape_detection::Plane<Traits> Plane_shape; ///< shape type.

  #ifdef DOXYGEN_RUNNING

    typedef unspecified_type Shape_range;
    ///< An `Iterator_range` with a bidirectional constant iterator type with value type `boost::shared_ptr<Shape>`.

    typedef unspecified_type Plane_range;
    ///< An `Iterator_range` with a bidirectional constant iterator type with value type `boost::shared_ptr<Plane_shape>`.

  #else

    struct Shape_range : public Iterator_range<
      typename std::vector<boost::shared_ptr<Shape> >::const_iterator> {
      typedef Iterator_range<
        typename std::vector<boost::shared_ptr<Shape> >::const_iterator> Base;

      Shape_range(boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > >
      extracted_shapes) :
      Base(make_range(extracted_shapes->begin(),
      extracted_shapes->end())), m_extracted_shapes(extracted_shapes) { }

    private:
      boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > >
      m_extracted_shapes; // keeps a reference to the shape vector
    };

    struct Plane_range : public Iterator_range<
      typename std::vector<boost::shared_ptr<Plane_shape> >::const_iterator> {
      typedef Iterator_range<
        typename std::vector<boost::shared_ptr<Plane_shape> >::const_iterator> Base;

      Plane_range(boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > >
        extracted_shapes) : Base(make_range(extracted_shapes->begin(),
        extracted_shapes->end())), m_extracted_shapes(extracted_shapes) { }

    private:
      boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > >
        m_extracted_shapes; // keeps a reference to the shape vector
    };

  #endif

    /// \cond SKIP_IN_MANUAL
    struct Filter_unassigned_points {
      Filter_unassigned_points() : m_shape_index(dummy) {}
      Filter_unassigned_points(const std::vector<int> &shapeIndex)
        : m_shape_index(shapeIndex) { }

      bool operator()(std::size_t x) {
        if (x < m_shape_index.size())
          return m_shape_index[x] == -1;
        else return true; // to prevent infinite incrementing
      }
      const std::vector<int>&  m_shape_index;
      std::vector<int> dummy;
    };

    typedef boost::filter_iterator<Filter_unassigned_points,
      boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t> >
      Point_index_iterator;
    ///< iterator for indices of points.
    /// \endcond

#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Point_index_range;
    ///< An `Iterator_range` with a bidirectional iterator with value type `std::size_t`
    ///  as indices into the input data that has not been assigned to a shape.
    ///  As this range class has no `size()` method, the method
    ///  `Efficient_RANSAC::number_of_unassigned_points()` is provided.
#else
    typedef Iterator_range<Point_index_iterator> Point_index_range;
#endif

    /// @}

    /// \name Parameters
    /// @{

    /*!
      %Parameters for the shape detection algorithm. They are explained in detail
      in Section \ref Shape_detection_RANSACParameters of the User Manual.
    */
    struct Parameters {
      Parameters()
        : min_points((std::numeric_limits<std::size_t>::max)())
        , epsilon(-1)
        , normal_threshold((FT) 0.9)
        , cluster_epsilon(-1)
      { }

      std::size_t min_points; ///< Minimum number of points of a shape. %Default value: 1% of total number of input points.
      FT epsilon;             ///< Maximum tolerance Euclidean distance from a point and a shape. %Default value: 1% of bounding box diagonal.
      FT normal_threshold;          ///< Maximum tolerance normal deviation from a point's normal to the normal on shape at projected point. %Default value: 0.9 (around 25 degrees).
      FT cluster_epsilon;            ///< Maximum distance between points to be considered connected. %Default value: 1% of bounding box diagonal.
    };

    /// @}

  private:

    class My_point_map {

      Input_iterator input_iterator_first;
      Point_map point_map;

    public:
      typedef typename Point_map::value_type value_type;
      typedef typename Point_map::reference reference;
      typedef std::size_t key_type;
      typedef boost::lvalue_property_map_tag category;

      My_point_map () { }
      My_point_map (Input_iterator first, Point_map map)
        : input_iterator_first (first), point_map (map) { }

      reference operator[](key_type k) const {
        return get(point_map, *(input_iterator_first + k));
      }

      friend reference get (const My_point_map& pmap, key_type idx) {
        return pmap[idx];
      }
    };

    typedef typename Traits::Search_traits Search_traits_base;
    typedef Search_traits_adapter <std::size_t, My_point_map, Search_traits_base> Search_traits;
    typedef CGAL::Fuzzy_sphere<Search_traits> Sphere;
    typedef CGAL::Kd_tree<Search_traits> Tree;
    typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_sphere;
    typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
    typedef typename Neighbor_search::Distance Distance;

    class Sort_by_planarity
    {
      Input_iterator m_first;
      Point_map m_point_map;
      My_point_map m_index_map;
      Tree& m_tree;
      FT m_cluster_epsilon;
      mutable boost::shared_ptr<std::vector<std::size_t> > nb_points;
      mutable boost::shared_ptr<std::vector<FT> > score;

    public:
      Sort_by_planarity (Input_iterator first,
                         Point_map point_map,
                         My_point_map index_map,
                         std::size_t size, Tree& tree, FT cluster_epsilon)
        : m_first (first)
        , m_point_map (point_map)
        , m_index_map (index_map)
        , m_tree (tree)
        , m_cluster_epsilon (cluster_epsilon)
        , nb_points (new std::vector<std::size_t>(size, 0))
        , score (new std::vector<FT>(size, -1.)) { }

      void compute_score (const std::size_t& idx) const
      {
        static std::vector<std::size_t> neighbors;

        neighbors.clear();
        Sphere fs (get(m_point_map, *(m_first + idx)), m_cluster_epsilon * 2, 0, m_tree.traits());
        m_tree.search (std::back_inserter (neighbors), fs);

        (*nb_points)[idx] = neighbors.size();

        Plane dummy;
        (*score)[idx]
          = CGAL::linear_least_squares_fitting_3
          (boost::make_transform_iterator (neighbors.begin(),
                                           CGAL::Property_map_to_unary_function<My_point_map>(m_index_map)),
           boost::make_transform_iterator (neighbors.end(),
                                           CGAL::Property_map_to_unary_function<My_point_map>(m_index_map)),
           dummy,
           CGAL::Dimension_tag<0>());
      }

      bool operator() (const std::size_t& a, const std::size_t& b) const
      {
        if ((*score)[a] == -1.)
          compute_score(a);
        if ((*score)[b] == -1.)
          compute_score(b);

        // if ((*nb_points)[a] != (*nb_points)[b])
        //   return (*nb_points)[a] > (*nb_points)[b];
        return (*score)[a] > (*score)[b];
      }
    };

    // Creates a function pointer for instancing shape instances.
    template <class ShapeT>
    static Shape *factory() {
      return new ShapeT;
    }

    Parameters m_options;

    // Traits class.
    Traits m_traits;

    // Maps index into points to assigned extracted primitive.
    std::vector<int> m_shape_index;
    std::size_t m_num_available_points;
    std::size_t m_num_total_points;

    // Give the index of the subset of point i.
    std::vector<int> m_index_subsets;

    boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > > m_extracted_shapes;

    std::vector<Shape *(*)()> m_shape_factories;

    // Iterators of input data.
    bool m_valid_iterators;
    Input_iterator m_input_iterator_first, m_input_iterator_beyond;
    Point_map m_point_map;
    Normal_map m_normal_map;
    My_point_map m_index_map;
    Tree* m_tree;

  public:

    /// \name Initialization
    /// @{

    /*!
      Constructs an empty shape detection object.
    */
    Region_growing_depr (Traits t = Traits())
      : m_traits(t)
      , m_num_available_points(0)
      , m_num_total_points(0)
      , m_valid_iterators(false)
      , m_tree (nullptr)
    {}

    /*!
      Releases all memory allocated by this instances including shapes.
    */
    ~Region_growing_depr() {
      clear();
    }

    /*!
      Retrieves the traits class.
     */
    const Traits&
    traits() const
    {
      return m_traits;
    }

    /*!
      Retrieves the point property map.
    */
    const Point_map& point_map() const { return m_point_map; }

    /*!
      Retrieves the normal property map.
    */
    const Normal_map& normal() const { return m_normal_map; }

    Input_iterator input_iterator_first() const
    {
      return m_input_iterator_first;
    }

    Input_iterator input_iterator_beyond() const
    {
      return m_input_iterator_beyond;
    }

    /*!
      Sets the input data. The range must stay valid
      until the detection has been performed and the access to the
      results is no longer required. This function first calls `clear()`.
    */
    void set_input(
      Input_range& input_range,
      ///< range of input data.
      Point_map point_map = Point_map(),
      ///< property map to access the position of an input point.
      Normal_map normal_map = Normal_map()
      ///< property map to access the normal of an input point.
      )
    {
        m_point_map = point_map;
        m_normal_map = normal_map;

        m_input_iterator_first = input_range.begin();
        m_input_iterator_beyond = input_range.end();

        m_index_map = My_point_map (m_input_iterator_first, m_point_map);
        clear();

        m_extracted_shapes =
          boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();

        m_num_available_points = m_num_total_points = std::distance(
          m_input_iterator_first, m_input_iterator_beyond);

        m_valid_iterators = true;
    }
    /*!
      Registers in the detection engine the shape type `ShapeType` that must inherit from `Shape_base`.
      For example, for registering a plane as detectable shape you should call
      `region_growing.add_shape_factory< Shape_detection::Plane<Traits> >();`.

      Note that if your call is within a template, you should add the
      `template` keyword just before `add_shape_factory`:
      `region_growing.template add_shape_factory<
      Shape_detection::Plane<Traits> >();`.

      \note So far, region growing algorithm only supports plane detection.
    */
    template <class Shape_type>
    void add_shape_factory() {
      CGAL_static_assertion_msg ((boost::is_convertible<Shape_type, Plane_shape>::value),
                                 "Region growing only supports Plane shapes.");

      m_shape_factories.push_back(factory<Shape_type>);
    }

    /*!
      Constructs internal data structures required for the shape detection.
      These structures only depend on the input data, i.e. the points and
      normal vectors. This method is called by `detect()`, if it was not called
      before by the user.
    */
    bool preprocess() {
      if (m_num_total_points == 0)
        return false;

      m_tree = new Tree (boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(0),
                         boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(m_num_total_points),
                         typename Tree::Splitter(),
                         m_index_map);
      m_tree->build();
      return true;
    }
    /// @}

    /// \name Memory Management
    /// @{
    /*!
      Removes all shape types registered for detection.
     */
    void clear_shape_factories() {
      m_shape_factories.clear();
    }

    /*!
      Removes all detected shapes.
      All internal structures are cleaned, including formerly detected shapes.
      Thus iterators and ranges retrieved through `shapes()`, `planes()` and `indices_of_unassigned_points()`
      are invalidated.
    */
    void clear()
    {
      // If there is no data yet, there are no data structures.
      if (!m_valid_iterators)
        return;

      if (m_tree != nullptr)
      {
        delete m_tree;
        m_tree = nullptr;
      }

      std::vector<int>().swap(m_shape_index);

      m_extracted_shapes =
        boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();

      m_num_available_points = m_num_total_points;
    }

    /// @}

    /// \name Detection
    /// @{

    /*!
      Performs the shape detection.

      \param options %Parameters for shape detection.

      \param callback can be omitted if the algorithm should be run
      without any callback. It is called regularly when the algorithm
      is running: the current advancement (between 0. and 1.) is
      passed as parameter. If it returns `true`, then the algorithm
      continues its execution normally; if it returns `false`, the
      algorithm is stopped. Note that this interruption may leave the
      class in an invalid state.

      \return `true` if shape types have been registered and
              input data has been set. Otherwise, `false` is returned.
    */
    bool detect(const Parameters &options = Parameters(),
                const std::function<bool(double)>& callback
                = std::function<bool(double)>())
    {
      // No shape types for detection or no points provided, exit.
      if (m_shape_factories.size() == 0 ||
          (m_input_iterator_beyond - m_input_iterator_first) == 0)
        return false;

      if (!preprocess())
        return false;

      // Reset data structures possibly used by former search.
      m_extracted_shapes =
        boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();
      m_num_available_points = m_num_total_points;

      m_options = options;

      Bbox_3 bbox = CGAL::bbox_3
        (boost::make_transform_iterator (m_input_iterator_first,
                                         CGAL::Property_map_to_unary_function<Point_map>(m_point_map)),
         boost::make_transform_iterator (m_input_iterator_beyond,
                                         CGAL::Property_map_to_unary_function<Point_map>(m_point_map)));

      FT bbox_diagonal = (FT) CGAL::sqrt(
          (bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
        + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
        + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));

      // Epsilon or cluster_epsilon have been set by the user?
      // If not, derive from bounding box diagonal.
      m_options.epsilon = (m_options.epsilon < 0)
        ? bbox_diagonal * (FT) 0.01 : m_options.epsilon;

      m_options.cluster_epsilon = (m_options.cluster_epsilon < 0)
        ? bbox_diagonal * (FT) 0.01 : m_options.cluster_epsilon;

      // Minimum number of points has been set?
      m_options.min_points =
        (m_options.min_points >= m_num_available_points) ?
          (std::size_t)((FT)0.01 * m_num_available_points) :
          m_options.min_points;
      m_options.min_points = (m_options.min_points < 10) ? 10 : m_options.min_points;

      // Initializing the shape index
      m_shape_index.assign(m_num_available_points, -1);

      Distance tr_dist (m_index_map);

      // Initialization structures.
      int class_index = -1;

      if (callback && !callback(0.))
        return false;

      std::vector<std::size_t> neighbors;

      std::vector<std::size_t> sorted_indices (m_num_total_points);
      std::copy (boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(0),
                 boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(m_num_total_points),
                 sorted_indices.begin());
#define SORT_INDICES
#ifdef SORT_INDICES
      std::sort (sorted_indices.begin(),
                 sorted_indices.end(),
                 Sort_by_planarity(m_input_iterator_first, m_point_map, m_index_map,
                                   m_num_total_points, *m_tree, m_options.cluster_epsilon));
#endif
      std::vector<double> fits(m_num_total_points, -1.);
      m_num_available_points = m_num_total_points;
      std::vector<std::size_t> index_container;
      std::vector<std::size_t> index_container_former_ring;
      std::set<std::size_t> index_container_current_ring;

      std::size_t num_to_test = 1 + m_num_total_points - m_options.min_points;
      for (std::size_t I = 0; I < num_to_test; ++ I)
      {
        if (callback && !callback(1. - (m_num_available_points / double(m_num_total_points))))
          return false;

        std::size_t i = sorted_indices[I];

        Input_iterator it = m_input_iterator_first + i;

        if (m_shape_index[i] != -1)
          continue;

        m_shape_index[i] = ++ class_index;

        int conti = 0; // to accelerate least_square fitting

        Vector plane_normal = get(m_normal_map, *it);
        plane_normal = plane_normal / std::sqrt(plane_normal * plane_normal);

        Plane optimal_plane(get(m_point_map, *it),
                            plane_normal);

        // initialization containers

        index_container.clear();
        index_container.push_back(i);

        index_container_former_ring.clear();
        index_container_former_ring.push_back(i);

        index_container_current_ring.clear();

        // propagation
        bool propagation = true;
        do
        {
          propagation = false;

          for (typename std::vector<std::size_t>::iterator icfrit = index_container_former_ring.begin();
               icfrit != index_container_former_ring.end(); ++ icfrit)
          {
            std::size_t point_index = *icfrit;

            Input_iterator pit = m_input_iterator_first + point_index;

            neighbors.clear();
            Sphere fs (get(m_point_map, *pit), m_options.cluster_epsilon, 0, m_tree->traits());
            m_tree->search (std::back_inserter (neighbors), fs);

            for (std::size_t nb = 0; nb < neighbors.size(); ++ nb)
            {
              std::size_t neighbor_index = neighbors[nb];
              Input_iterator nbit = m_input_iterator_first + neighbor_index;

              if (m_shape_index[neighbor_index] != -1)
                continue;

              const Point& neighbor = get(m_point_map, *nbit);
              double distance = CGAL::squared_distance(neighbor, optimal_plane);

              Vector normal = get(m_normal_map, *nbit);
              normal = normal / std::sqrt (normal * normal);

              if (distance > m_options.epsilon * m_options.epsilon
                  || std::fabs(normal * plane_normal) < m_options.normal_threshold)
                continue;

              m_shape_index[neighbor_index] = class_index;
              propagation = true;
              index_container_current_ring.insert(neighbor_index);
            }
          }

          // update containers
          index_container_former_ring.clear();
          index_container_former_ring.reserve (index_container_current_ring.size());
          index_container.reserve (index_container.size() + index_container_current_ring.size());
          for (typename std::set<std::size_t>::iterator lit = index_container_current_ring.begin();
               lit != index_container_current_ring.end(); ++lit)
          {
            index_container_former_ring.push_back(*lit);
            index_container.push_back(*lit);
          }
          index_container_current_ring.clear();

          conti++;
          if (index_container.size() < 5)
            continue;

          if ((conti < 10) || (conti<50 && conti % 10 == 0) || (conti>50 && conti % 500 == 0))
          {
            std::list<Point> listp;
            for (typename std::vector<std::size_t>::iterator icit = index_container.begin();
                 icit != index_container.end(); ++ icit)
              listp.push_back(get (m_point_map, *(m_input_iterator_first + *icit)));

            Plane reajusted_plane;
            CGAL::linear_least_squares_fitting_3(listp.begin(),
                                                 listp.end(),
                                                 reajusted_plane,
                                                 CGAL::Dimension_tag<0>());
            optimal_plane = reajusted_plane;
            plane_normal = optimal_plane.orthogonal_vector();
            plane_normal = plane_normal / std::sqrt(plane_normal * plane_normal);
          }
        }
        while (propagation);

        if (index_container.size() >= m_options.min_points)
        {
          Shape *p = (Shape *) (*(m_shape_factories.begin()))();

          // BUG something fishy around here.
          m_num_available_points -= index_container.size();

          p->compute (index_container,  m_input_iterator_first, m_traits,
                      m_point_map, m_normal_map, m_options.epsilon, m_options.normal_threshold);
          p->m_indices.clear();
          std::copy (index_container.begin(), index_container.end(),
                     std::back_inserter (p->m_indices));

          Plane_shape* ps = dynamic_cast<Plane_shape*>(p);
          CGAL_assume (ps != nullptr);
          ps->update (optimal_plane);
          m_extracted_shapes->push_back (boost::shared_ptr<Shape>(p));
        }
        else
        {
          class_index--;
          m_shape_index[i] = -1;
          for (typename std::vector<std::size_t>::iterator icit = index_container.begin();
               icit != index_container.end(); ++ icit)
            m_shape_index[*icit] = -1;
        }
      }
      return true;
    }

    /// @}

    /// \name Access
    /// @{
    /*!
      Returns an `Iterator_range` with a bidirectional iterator with value type
      `boost::shared_ptr<Shape>` over the detected shapes in the order of detection.

      \note So far, region growing algorithm only supports plane
      detection, so this method is equivalent to `planes()` except
      that it returns planes with the abstract type `Shape`.
    */
    Shape_range shapes() const {
      return Shape_range(m_extracted_shapes);
    }

    /*!
      Returns an `Iterator_range` with a bidirectional iterator with
      value type `boost::shared_ptr<Plane_shape>` over only the
      detected planes in the order of detection.
    */
    Plane_range planes() const {
      boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > > planes
        = boost::make_shared<std::vector<boost::shared_ptr<Plane_shape> > >();

      for (std::size_t i = 0; i < m_extracted_shapes->size(); ++ i)
      {
        boost::shared_ptr<Plane_shape> pshape
          = boost::dynamic_pointer_cast<Plane_shape>((*m_extracted_shapes)[i]);

        // Ignore all shapes other than plane.
        if (pshape != boost::shared_ptr<Plane_shape>())
          planes->push_back (pshape);
      }
      return Plane_range(planes);
    }

    /*!
      Number of points not assigned to a shape.
    */
    std::size_t number_of_unassigned_points() {
      return m_num_available_points;
    }

    /*!
      Returns an `Iterator_range` with a bidirectional iterator with the value type
      `std::size_t` as indices into the input data that has not been assigned to a shape.
    */
    Point_index_range indices_of_unassigned_points() {
      Filter_unassigned_points fup(m_shape_index);

      Point_index_iterator p1 =
        boost::make_filter_iterator<Filter_unassigned_points>(
        fup,
        boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(0),
        boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(m_shape_index.size()));

      return make_range(p1, Point_index_iterator(p1.end()));
    }
    /// @}
  };

  #ifdef DOXYGEN_NS
    } // namespace deprecated
  #endif

} // namespace Shape_detection
} // namespace CGAL

// --<

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_DEPRECATED_H
