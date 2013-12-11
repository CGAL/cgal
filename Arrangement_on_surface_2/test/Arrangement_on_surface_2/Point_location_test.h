#ifndef CGAL_POINT_LOCATION_TEST_H
#define CGAL_POINT_LOCATION_TEST_H

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_specified_points_generator.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arr_triangulation_point_location.h>

#include "IO_test.h"

/*! Point location test */
template <typename Geom_traits_T, typename Topol_traits_T>
class Point_location_test : public IO_test<Geom_traits_T> {
public:
  typedef Geom_traits_T                                 Geom_traits;
  typedef Topol_traits_T                                Topol_traits;
  typedef IO_test<Geom_traits>                          Base;

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename Base::Points_vector                  Points_vector;
  typedef typename Base::Xcurves_vector                 Xcurves_vector;
  typedef typename Base::Curves_vector                  Curves_vector;

  typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                        Arrangement;

  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;

  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;

  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;

  typedef typename Points_vector::iterator              Point_iterator;

  typedef std::vector<CGAL::Object>                     Objects_vector;
  typedef Objects_vector::iterator                      Object_iterator;

  typedef typename boost::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>    Cell_handle;
  typedef std::vector<Cell_handle>                      Variants_vector;
  typedef typename Variants_vector::iterator            Variant_iterator;

protected:
  typedef typename CGAL::Arr_naive_point_location<Arrangement>
                                                    Naive_point_location;
  typedef typename CGAL::Arr_simple_point_location<Arrangement>
                                                    Simple_point_location;
  typedef typename CGAL::Arr_walk_along_line_point_location<Arrangement>
                                                    Walk_point_location;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement>
                                                    Lm_point_location;
  typedef typename CGAL::Arr_random_landmarks_generator<Arrangement>
                                                    Random_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Random_lm_generator>
                                                    Lm_random_point_location;
  typedef typename CGAL::Arr_grid_landmarks_generator<Arrangement>
                                                    Grid_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Grid_lm_generator>
                                                    Lm_grid_point_location;
  typedef typename CGAL::Arr_halton_landmarks_generator<Arrangement>
                                                    Halton_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Halton_lm_generator>
                                                    Lm_halton_point_location;
  typedef typename CGAL::Arr_middle_edges_landmarks_generator<Arrangement>
                                                    Middle_edges_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Middle_edges_generator>
    Lm_middle_edges_point_location;
  typedef typename CGAL::Arr_landmarks_specified_points_generator<Arrangement>
    Specified_points_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Specified_points_generator>
    Lm_specified_points_point_location;
  typedef typename CGAL::Arr_trapezoid_ric_point_location<Arrangement>
    Trapezoid_ric_point_location;

  typedef CGAL::Arr_triangulation_point_location<Arrangement>
    Triangulation_point_location;

  // ===> Add new point location type here <===

protected:
  /*! The geometry traits */
  const Geom_traits& m_geom_traits;

  /*! The arrangement */
  Arrangement* m_arr;

  /*! The input data file of the query points*/
  std::string m_filename_queries;

  /*! The query points */
  Points_vector m_query_points;

  Naive_point_location* m_naive_pl;                              // 0
  Simple_point_location* m_simple_pl;                            // 1
  Walk_point_location* m_walk_pl;                                // 2
  Lm_point_location* m_lm_pl;                                    // 3
  Lm_random_point_location* m_random_lm_pl;                      // 4
  Lm_grid_point_location* m_grid_lm_pl;                          // 5
  Lm_halton_point_location* m_halton_lm_pl;                      // 6
  Lm_middle_edges_point_location* m_middle_edges_lm_pl;          // 7
  Lm_specified_points_point_location* m_specified_points_lm_pl;  // 8
  Triangulation_point_location* m_triangulation_pl;              // 9
  Trapezoid_ric_point_location* m_trapezoid_ric_pl;              // 10
  Trapezoid_ric_point_location* m_trapezoid_ric_no_grnt_pl;      // 11

  Random_lm_generator* m_random_g;
  Grid_lm_generator* m_grid_g;
  Halton_lm_generator* m_halton_g;
  Middle_edges_generator* m_middle_edges_g;
  Specified_points_generator* m_specified_points_g;

  // // ===> Change the number of point-location startegies
  // //      when a new point location is added. <===
  #define MAX_NUM_POINT_LOCATION_STRATEGIES 12

  int verify(Objects_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES],
             size_t size, size_t pls_num);

  int verify(Variants_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES],
             size_t size, size_t pls_num);

public:
  /*! Constructor */
  Point_location_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Point_location_test()
  {
    deallocate_arrangement();
    clear();
    deallocate_pl_strategies();
  }

  void set_filenames(const char* points_filename, const char* xcurves_filename,
                     const char* curves_filename, const char* queries_filename);

  /*! Initialize the data structures */
  virtual bool init();

  /*! Perform the test */
  virtual bool perform();

  /*! Clear the data structures */
  virtual void clear();

  bool allocate_arrangement();

  void deallocate_arrangement();

  bool construct_arrangement();

  void clear_arrangement();

  bool allocate_pl_strategies();

  bool construct_pl_strategies();

  bool attach_pl_strategies();

  void deallocate_pl_strategies();

  template <typename Point_location,
            typename InputIterator, typename OutputIterator>
  void query(Point_location& pl, const char* type,
             InputIterator begin, InputIterator end, OutputIterator oi)
  {
    typedef InputIterator                       Input_iterator;

    CGAL::Timer timer;
    timer.reset(); timer.start();
    Input_iterator piter;
    for (piter = begin; piter != end; ++piter) {
      Point_2 q = (*piter);
      *oi++ = pl.locate(q);
    }
    timer.stop();
    std::cout << type << " location took " << timer.time() << std::endl;
  }
};

/*!
 * Constructor.
 */
template <typename Geom_traits_T, typename Topol_traits_T>
Point_location_test<Geom_traits_T, Topol_traits_T>::
Point_location_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_arr(NULL),
  m_naive_pl(NULL),
  m_simple_pl(NULL),
  m_walk_pl(NULL),
  m_lm_pl(NULL),
  m_random_lm_pl(NULL),
  m_grid_lm_pl(NULL),
  m_halton_lm_pl(NULL),
  m_middle_edges_lm_pl(NULL),
  m_specified_points_lm_pl(NULL),
  m_triangulation_pl(NULL),
  m_trapezoid_ric_pl(NULL),
  m_trapezoid_ric_no_grnt_pl(NULL),
  m_random_g(NULL),
  m_grid_g(NULL),
  m_halton_g(NULL),
  m_middle_edges_g(NULL),
  m_specified_points_g(NULL)
{}

/*! Set the file names */
template <typename Geom_traits_T, typename Topol_traits_T>
void Point_location_test<Geom_traits_T, Topol_traits_T>::
set_filenames(const char* points_filename,
              const char* xcurves_filename,
              const char* curves_filename,
              const char* queries_filename)
{
  Base::set_filenames(points_filename, xcurves_filename, curves_filename);
  m_filename_queries.assign(queries_filename);
}

/*! Initialize the data structures */
template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::init()
{
  if (!Base::init()) return false;

  // Read the query points
  if (!this->read_points(m_filename_queries.c_str(), m_query_points))
    return false;

  return true;
}

/*! Clear the data structures */
template <typename Geom_traits_T, typename Topol_traits_T>
void Point_location_test<Geom_traits_T, Topol_traits_T>::clear()
{
  Base::clear();
  m_query_points.clear();
  m_filename_queries.clear();
}

/*! Clear the data structures */
template <typename Geom_traits_T, typename Topol_traits_T>
void Point_location_test<Geom_traits_T, Topol_traits_T>::
deallocate_pl_strategies()
{
  if (m_naive_pl) {
    delete m_naive_pl;
    m_naive_pl = NULL;
  }
  if (m_simple_pl) {
    delete m_simple_pl;
    m_simple_pl = NULL;
  }
  if (m_walk_pl) {
    delete m_walk_pl;
    m_walk_pl = NULL;
  }

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)
  if (m_lm_pl) {
    delete m_lm_pl;
    m_lm_pl = NULL;
  }
  if (m_random_lm_pl) {
    delete m_random_lm_pl;
    m_random_lm_pl = NULL;
  }
  if (m_grid_lm_pl) {
    delete m_grid_lm_pl;
    m_grid_lm_pl = NULL;
  }
  if (m_halton_lm_pl) {
    delete m_halton_lm_pl;
    m_halton_lm_pl = NULL;
  }
  if (m_middle_edges_lm_pl) {
    delete m_middle_edges_lm_pl;
    m_middle_edges_lm_pl = NULL;
  }
  if (m_specified_points_lm_pl) {
    delete m_specified_points_lm_pl;
    m_specified_points_lm_pl = NULL;
  }

  // Free Generators
  if (m_random_g) {
    delete m_random_g;
    m_random_g = NULL;
  }
  if (m_grid_g) {
    delete m_grid_g;
    m_grid_g = NULL;
  }
  if (m_halton_g) {
    delete m_halton_g;
    m_halton_g = NULL;
  }
  if (m_middle_edges_g) {
    delete m_middle_edges_g;
    m_middle_edges_g = NULL;
  }
  if (m_specified_points_g) {
    delete m_specified_points_g;
    m_specified_points_g = NULL;
  }
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  if (m_triangulation_pl) {
    delete m_triangulation_pl;
    m_triangulation_pl = NULL;
  }
#endif

  if (m_trapezoid_ric_pl) {
    delete m_trapezoid_ric_pl;
    m_trapezoid_ric_pl = NULL;
  }
  if (m_trapezoid_ric_no_grnt_pl) {
    delete m_trapezoid_ric_no_grnt_pl;
    m_trapezoid_ric_no_grnt_pl = NULL;
  }
}

template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::allocate_arrangement()
{
  if (!(m_arr = new Arrangement(&m_geom_traits))) return false;
  return true;
}

template <typename Geom_traits_T, typename Topol_traits_T>
void Point_location_test<Geom_traits_T, Topol_traits_T>::
deallocate_arrangement()
{
  if (m_arr) {
    delete m_arr;
    m_arr = NULL;
  }
}

template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::
construct_arrangement()
{
  // Insert all into the arrangement
  CGAL::insert(*m_arr, this->m_xcurves.begin(), this->m_xcurves.end());
  // insert(*m_arr, m_points.begin(), m_points.end());
  CGAL::insert(*m_arr, this->m_curves.begin(), this->m_curves.end());

  // Print the size of the arrangement.
  std::cout << "V = " << m_arr->number_of_vertices()
            << ",  E = " << m_arr->number_of_edges()
            << ",  F = " << m_arr->number_of_faces() << std::endl;

  return true;
}

template <typename Geom_traits_T, typename Topol_traits_T>
void Point_location_test<Geom_traits_T, Topol_traits_T>::clear_arrangement()
{
  if (m_arr) m_arr->clear();
}

template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::
allocate_pl_strategies()
{
  // Allocate all point location strategies.
  if (!(m_naive_pl = new Naive_point_location())) return false;
  if (!(m_simple_pl = new Simple_point_location())) return false;
  if (!(m_walk_pl = new Walk_point_location())) return false;

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  if (!(m_lm_pl = new Lm_point_location())) return false;
  if (!(m_random_lm_pl = new Lm_random_point_location())) return false;
  if (!(m_grid_lm_pl = new Lm_grid_point_location())) return false;
  if (!(m_halton_lm_pl = new Lm_halton_point_location())) return false;
  if (!(m_middle_edges_lm_pl = new Lm_middle_edges_point_location()))
    return false;
  if (!(m_specified_points_lm_pl = new Lm_specified_points_point_location()))
    return false;
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  if (!(m_triangulation_pl = new Triangulation_point_location())) return false;
#endif

  if (!(m_trapezoid_ric_pl = new Trapezoid_ric_point_location())) return false;
  if (!(m_trapezoid_ric_no_grnt_pl = new Trapezoid_ric_point_location()))
    return false;

  // ===> Add new point location instance here. <===
  return true;
}

template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::
construct_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  m_naive_pl = new Naive_point_location(*m_arr);                        // 0
  m_simple_pl = new Simple_point_location(*m_arr);                      // 1
  m_walk_pl = new Walk_point_location(*m_arr);                          // 2

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  timer.reset(); timer.start();
  m_lm_pl = new Lm_point_location(*m_arr);                              // 3
  timer.stop();
  std::cout << "Lm (vert) construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_random_g = new Random_lm_generator(*m_arr);
  m_random_lm_pl = new Lm_random_point_location(*m_arr, m_random_g);    // 4
  timer.stop();
  std::cout << "Random lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_grid_g = new Grid_lm_generator(*m_arr);
  m_grid_lm_pl = new Lm_grid_point_location(*m_arr, m_grid_g);          // 5
  timer.stop();
  std::cout << "Grid lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_halton_g = new Halton_lm_generator(*m_arr);
  m_halton_lm_pl = new Lm_halton_point_location(*m_arr, m_halton_g);    // 6
  timer.stop();
  std::cout << "Halton lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_middle_edges_g = new Middle_edges_generator(*m_arr);
  m_middle_edges_lm_pl =
    new Lm_middle_edges_point_location(*m_arr, m_middle_edges_g);       // 7
  timer.stop();
  std::cout << "Middle edges lm construction took " << timer.time()
            << std::endl;

  timer.reset(); timer.start();
  m_specified_points_g = new Specified_points_generator(*m_arr);
  m_specified_points_lm_pl =
    new Lm_specified_points_point_location(*m_arr, m_specified_points_g); // 8
  timer.stop();
  std::cout << "Specified_points lm construction took "
            << timer.time() << std::endl;
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  timer.reset(); timer.start();
  m_triangulation_pl = new Triangulation_point_location(*m_arr);        // 9
  timer.stop();
  std::cout << "Triangulation lm construction took "
            << timer.time() << std::endl;
#endif

  timer.reset(); timer.start();
  m_trapezoid_ric_pl = new Trapezoid_ric_point_location(*m_arr);        // 10
  timer.stop();
  std::cout << "Trapezoid RIC construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_trapezoid_ric_no_grnt_pl =
    new Trapezoid_ric_point_location(*m_arr, false);                    // 11
  timer.stop();
  std::cout << "Trapezoid RIC without-guarantees construction took "
            << timer.time() << std::endl;

  std::cout << std::endl;

  // ===> Add new point location instance here. <===
  return true;
}

template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::attach_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  m_naive_pl->attach(*m_arr);
  m_simple_pl->attach(*m_arr);
  m_walk_pl->attach(*m_arr);

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  timer.reset(); timer.start();
  m_lm_pl->attach(*m_arr);
  timer.stop();
  std::cout << "Lm (vert) construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_random_g = new Random_lm_generator(*m_arr);
  m_random_lm_pl->attach(*m_arr, m_random_g);
  timer.stop();
  std::cout << "Random lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_grid_g = new Grid_lm_generator(*m_arr);
  m_grid_lm_pl->attach(*m_arr, m_grid_g);
  timer.stop();
  std::cout << "Grid lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_halton_g = new Halton_lm_generator(*m_arr);
  m_halton_lm_pl->attach(*m_arr, m_halton_g);
  timer.stop();
  std::cout << "Halton lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_middle_edges_g = new Middle_edges_generator(*m_arr);
  m_middle_edges_lm_pl->attach(*m_arr, m_middle_edges_g);
  timer.stop();
  std::cout << "Middle edges lm construction took " << timer.time()
            << std::endl;

  timer.reset(); timer.start();
  m_specified_points_g = new Specified_points_generator(*m_arr);
  m_specified_points_lm_pl->attach(*m_arr, m_specified_points_g);
  timer.stop();
  std::cout << "Specified_points lm construction took "
            << timer.time() << std::endl;
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  timer.reset(); timer.start();
  m_triangulation_pl->attach(*m_arr);
  timer.stop();
  std::cout << "Triangulation lm construction took "
            << timer.time() << std::endl;
#endif

  timer.reset(); timer.start();
  m_trapezoid_ric_pl->attach(*m_arr);
  timer.stop();
  std::cout << "Trapezoid RIC construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  m_trapezoid_ric_no_grnt_pl->with_guarantees(false);
  m_trapezoid_ric_no_grnt_pl->attach(*m_arr);

  timer.stop();
  std::cout << "Trapezoid RIC without-guarantees construction took "
            << timer.time() << std::endl;

  std::cout << std::endl;

  // ===> Add new point location instance here. <===
  return true;
}

// Perform the test
template <typename Geom_traits_T, typename Topol_traits_T>
bool Point_location_test<Geom_traits_T, Topol_traits_T>::perform()
{
#if ((CGAL_ARR_POINT_LOCATION_VERSION < 2) || \
     defined(CGAL_ARR_POINT_LOCATION_CONVERSION))
  Objects_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES];
#else
  Variants_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES];
#endif

  // Locate the points in the list using all point location strategies.

  // std::cout << "Time in seconds" << std::endl;
  std::cout << std::endl;

  size_t pl_index = 0;

  query(*m_naive_pl, "Naive", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Naive

  query(*m_simple_pl, "Simple", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Simple

  query(*m_walk_pl, "Walk", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Walk

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  query(*m_lm_pl, "Landmarks (vertices)",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks (vertices)

  query(*m_random_lm_pl, "Landmarks random",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks random

  query(*m_grid_lm_pl, "Landmarks grid",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks grid

  query(*m_halton_lm_pl, "Landmarks Halton",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks Halton

  query(*m_middle_edges_lm_pl, "Landmarks middle edges",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks middle edges

  query(*m_specified_points_lm_pl, "Landmarks specified points",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Landmarks specified points
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  query(*m_triangulation_pl, "Triangulation",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Triangulation
#endif

  query(*m_trapezoid_ric_pl, "Trapezoidal RIC",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Trapezoidal RIC

  // Trapezoidal RIC without guarantees
  query(*m_trapezoid_ric_no_grnt_pl, "Trapezoidal RIC without guarantees",
        m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));

  std::cout << std::endl;

  // ===> Add a call to operate the new point location. <===

  // Number of point location strategies used.
  size_t pls_num = pl_index;
  std::cout << "Number of strategies is " << pls_num << std::endl;

  // End Location

  int result = 0;

  // get size of objects
  size_t size = objs[0].size();
  std::cout << "size is " << size << std::endl;

  for (pl_index = 0; pl_index < pls_num; ++pl_index) {
    if (size != objs[pl_index].size()) {
      std::cout << "Error: size of pl number " << pl_index << " is "
                << objs[pl_index].size() << std::endl;
      result += -1;
      break;
    }
  }

  result += verify(objs, size, pls_num);
  return (result == 0);
}

// Verify the results
template <typename Geom_traits_T, typename Topol_traits_T>
int Point_location_test<Geom_traits_T, Topol_traits_T>::
verify(Objects_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES],
       size_t size, size_t pls_num)
{
  Vertex_const_handle vh_ref, vh_cur;
  Halfedge_const_handle hh_ref, hh_cur;
  Face_const_handle fh_ref, fh_cur;

  int result = 0;

  // Assign and check results
  size_t qi; //qi is the query point index
  for (qi = 0; qi < size; ++qi) {
    // Assign object to a face
    if (CGAL::assign(fh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          if (fh_cur != fh_ref) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expecte: a face." << std::endl;
            std::cout << "Actual: a different face" << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expecte: a face." << std::endl;
        result += -1;
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      //if (fh_ref->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
      continue;
    }

    // Assign object to a halfedge
    if (CGAL::assign (hh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          if ((hh_cur != hh_ref) && (hh_cur->twin() != hh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a halfedge, " << hh_ref->curve()
                      << std::endl;
            std::cout << "Actual: a different halfedge: " << hh_cur->curve()
                      << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a halfedge, " << hh_ref->curve() << std::endl;
        result += -1;
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      continue;
    }

    // Assign object to a vertex
    if (CGAL::assign(vh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          if (vh_cur != vh_ref) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a vertex, "<< vh_ref->point() << std::endl;
            std::cout << "Actual: a different vertex, "<< vh_cur->point()
                      << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a vertex, "<< vh_ref->point() << std::endl;
        result += -1;
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknown object." << std::endl;
      }
      continue;
    }

    std::cout << "Error: Unknown!" << std::endl;
    result += -1;
  }
  return result;
}

// Verify the results
template <typename Geom_traits_T, typename Topol_traits_T>
int Point_location_test<Geom_traits_T, Topol_traits_T>::
verify(Variants_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES],
       size_t size, size_t pls_num)
{
  const Vertex_const_handle* vh_ref;
  const Halfedge_const_handle* hh_ref;
  const Face_const_handle* fh_ref;

  const Vertex_const_handle* vh_cur;
  const Halfedge_const_handle* hh_cur;
  const Face_const_handle* fh_cur;

  int result = 0;

  // Assign and check results
  size_t qi; //qi is the query point index
  for (qi = 0; qi < size; ++qi) {
    // Assign object to a face
    fh_ref = boost::get<Face_const_handle>(&(objs[0][qi]));
    if (fh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          if ((*fh_cur) != (*fh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a face." << std::endl;
            std::cout << "Actual: a different face." << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a face." << std::endl;
        result += -1;
	hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
	if (hh_cur) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
	vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      //if ((*fh_ref)->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
      continue;
    }

    // Assign object to a halfedge
    hh_ref = boost::get<Halfedge_const_handle>(&(objs[0][qi]));
    if (hh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
        if (hh_cur) {
          if (((*hh_cur) != (*hh_ref)) && ((*hh_cur)->twin() != (*hh_ref))) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a halfedge, " << (*hh_ref)->curve()
                      << std::endl;
            std::cout << "Actual: a different halfedge, " << (*hh_cur)->curve()
                      << std::endl;
            result += -1;
          }
          continue;
        }
        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a halfedge, " << (*hh_ref)->curve()
                  << std::endl;
        result += -1;
	fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
	vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      continue;
    }

    // Assign object to a vertex
    vh_ref = boost::get<Vertex_const_handle>(&(objs[0][qi]));
    if (vh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          if ((*vh_cur) != (*vh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a vertex: "<< (*vh_ref)->point()
                      << std::endl;
            std::cout << "Actual: a different vertex: "<< (*vh_cur)->point()
                      << std::endl;
            result += -1;
          }
          continue;
        }
        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a vertex: "<< (*vh_ref)->point() << std::endl;
        result += -1;
	fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
	hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
        if (hh_cur) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknown object." << std::endl;
      }
      continue;
    }

    std::cout << "Error: Unknown!" << std::endl;
    result += -1;
  }
  return result;
}

#endif
