#ifndef CGAL_POINT_LOCATION_TEST_H
#define CGAL_POINT_LOCATION_TEST_H

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Arrangement_2.h>
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
//#include <CGAL/Arr_triangulation_point_location.h>

#include "IO_test.h"

/*! Point location test */
template <typename T_Traits>
class Point_location_test : public IO_test<T_Traits> {
private:
  typedef T_Traits                                      Traits;
  typedef IO_test<Traits>                               Base;

public:
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename Base::Points_vector                  Points_vector;
  typedef typename Base::Xcurves_vector                 Xcurves_vector;
  typedef typename Base::Curves_vector                  Curves_vector;
  
  typedef CGAL::Arrangement_2<Traits>                   Arrangement_2;

  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;

  typedef typename Points_vector::iterator              Point_iterator;
  typedef std::vector<CGAL::Object>                     Objects_vector;
  typedef Objects_vector::iterator                      Object_iterator;
  
protected:
  typedef typename CGAL::Arr_naive_point_location<Arrangement_2>     
                                                    Naive_point_location;
  typedef typename CGAL::Arr_simple_point_location<Arrangement_2>     
                                                    Simple_point_location;
  typedef typename CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
                                                    Walk_point_location;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2> 
                                                    Lm_point_location;
  typedef typename CGAL::Arr_random_landmarks_generator<Arrangement_2>
                                                    Random_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2,
                                                      Random_lm_generator> 
                                                    Lm_random_point_location;
  typedef typename CGAL::Arr_grid_landmarks_generator<Arrangement_2>
                                                    Grid_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2,
                                                      Grid_lm_generator> 
                                                    Lm_grid_point_location;
  typedef typename CGAL::Arr_halton_landmarks_generator<Arrangement_2>
                                                    Halton_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2,
                                                      Halton_lm_generator> 
                                                    Lm_halton_point_location;
  typedef typename CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
                                                    Middle_edges_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2,
                                                      Middle_edges_generator> 
    Lm_middle_edges_point_location;
  typedef typename CGAL::Arr_landmarks_specified_points_generator<Arrangement_2>
    Specified_points_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2,
                                                      Specified_points_generator> 
    Lm_specified_points_point_location;
  typedef typename CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> 
    Trapezoid_ric_point_location;

  //   typedef CGAL::Arr_triangulation_point_location<Arrangement_2> 
  //     Triangulation_point_location;

  // ===> Add new point location type here <===

  Naive_point_location* m_naive_pl;                              // 0
  Simple_point_location* m_simple_pl;                            // 1
  Walk_point_location* m_walk_pl;                                // 2
  Lm_point_location* m_lm_pl;                                    // 3
  Lm_random_point_location* m_random_lm_pl;                      // 4
  Lm_grid_point_location* m_grid_lm_pl;                          // 5
  Lm_halton_point_location* m_halton_lm_pl;                      // 6
  Lm_middle_edges_point_location* m_middle_edges_lm_pl;          // 7
  Lm_specified_points_point_location* m_specified_points_lm_pl;  // 8
  // Triangulation_point_location m_triangulation_pl;            // 9
  Trapezoid_ric_point_location* m_trapezoid_ric_pl;              // 10
  Trapezoid_ric_point_location* m_trapezoid_ric_no_grnt_pl;      // 11

  Random_lm_generator* m_random_g;
  Grid_lm_generator* m_grid_g;
  Halton_lm_generator* m_halton_g;
  Middle_edges_generator* m_middle_edges_g;
  Specified_points_generator* m_specified_points_g;
  
  // // ===> Change the number of point-location startegies
  // //      when a new point location is added. <===
  #define MAX_NUM_POINT_LOCATION_STRATEGIES 11
  
public:
  /*! Constructor */
  Point_location_test();

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
      CGAL::Object obj = pl.locate(q);
      *oi++ = obj;
    }
    timer.stop();
    std::cout << type << " location took " << timer.time() << std::endl;
  }

  /*! The arrangement */
  Arrangement_2* m_arr;  

private:
  /*! The input data file of the query points*/
  std::string m_filename_queries;

  /*! The query points */
  Points_vector m_query_points;
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <typename T_Traits>
Point_location_test<T_Traits>::Point_location_test() :
  m_naive_pl(NULL), 
  m_simple_pl(NULL), 
  m_walk_pl(NULL), 
  m_lm_pl(NULL), 
  m_random_lm_pl(NULL), 
  m_grid_lm_pl(NULL), 
  m_halton_lm_pl(NULL), 
  m_middle_edges_lm_pl(NULL), 
  m_specified_points_lm_pl(NULL), 
  // m_triangulation_pl(NULL), 
  m_trapezoid_ric_pl(NULL), 
  m_trapezoid_ric_no_grnt_pl(NULL), 
  m_random_g(NULL), 
  m_grid_g(NULL), 
  m_halton_g(NULL), 
  m_middle_edges_g(NULL), 
  m_specified_points_g(NULL),
  m_arr(NULL)
{}

/*! Set the file names */
template <typename T_Traits>
void Point_location_test<T_Traits>::set_filenames(const char* points_filename,
                                                  const char* xcurves_filename,
                                                  const char* curves_filename,
                                                  const char* queries_filename)
{
  Base::set_filenames(points_filename, xcurves_filename, curves_filename);
  m_filename_queries.assign(queries_filename);
}

/*! Initialize the data structures */
template <typename T_Traits>
bool Point_location_test<T_Traits>::init()
{
  if (!Base::init()) return false;

  // Read the query points
  if (!this->read_points(m_filename_queries.c_str(), m_query_points))
    return false;

  return true;
}

/*! Clear the data structures */
template<class T_Traits>
void Point_location_test<T_Traits>::clear()
{
  Base::clear();
  m_query_points.clear();
  m_filename_queries.clear();
}

/*! Clear the data structures */
template<class T_Traits>
void Point_location_test<T_Traits>::deallocate_pl_strategies()
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

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)
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

  // if (m_triangulation_pl) {
  //   delete m_triangulation_pl;
  //   m_triangulation_pl = NULL;
  // }
  if (m_trapezoid_ric_pl) {
    delete m_trapezoid_ric_pl;
    m_trapezoid_ric_pl = NULL;
  }
  if (m_trapezoid_ric_no_grnt_pl) {
    delete m_trapezoid_ric_no_grnt_pl;
    m_trapezoid_ric_no_grnt_pl = NULL;
  }
}

template <typename T_Traits>
bool Point_location_test<T_Traits>::allocate_arrangement()
{
  if (!(m_arr = new Arrangement_2())) return false;
  return true;
}

template <typename T_Traits>
void Point_location_test<T_Traits>::deallocate_arrangement()
{
  if (m_arr) {
    delete m_arr;
    m_arr = NULL;
  }
}

template <typename T_Traits>
bool Point_location_test<T_Traits>::construct_arrangement()
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

template <typename T_Traits>
void Point_location_test<T_Traits>::clear_arrangement()
{
  if (m_arr) m_arr->clear();
}

template <typename T_Traits>
bool Point_location_test<T_Traits>::allocate_pl_strategies()
{
  // Allocate all point location strategies.
  if (!(m_naive_pl = new Naive_point_location())) return false;
  if (!(m_simple_pl = new Simple_point_location())) return false;
  if (!(m_walk_pl = new Walk_point_location())) return false;

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

  if (!(m_lm_pl = new Lm_point_location())) return false;
  if (!(m_random_lm_pl = new Lm_random_point_location())) return false;
  if (!(m_grid_lm_pl = new Lm_grid_point_location())) return false;
  if (!(m_halton_lm_pl = new Lm_halton_point_location())) return false;
  if (!(m_middle_edges_lm_pl = new Lm_middle_edges_point_location()))
    return false;
  if (!(m_specified_points_lm_pl = new Lm_specified_points_point_location()))
    return false;

  // if (!(m_triangulation_pl = new Triangulation_point_location()))
  //   return false;

#endif
  
  if (!(m_trapezoid_ric_pl = new Trapezoid_ric_point_location())) return false;
  if (!(m_trapezoid_ric_no_grnt_pl = new Trapezoid_ric_point_location()))
    return false;

  // ===> Add new point location instance here. <===
  return true;
}

template <typename T_Traits>
bool Point_location_test<T_Traits>::construct_pl_strategies()
{
  typedef T_Traits Traits;

  // Initialize all point location strategies.
  CGAL::Timer timer;

  m_naive_pl = new Naive_point_location(*m_arr);                        // 0
  m_simple_pl = new Simple_point_location(*m_arr);                      // 1
  m_walk_pl = new Walk_point_location(*m_arr);                          // 2

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

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

  // timer.reset(); timer.start();
  // m_triangulation_pl = new Triangulation_point_location(*m_arr);     // 9
  // timer.stop(); 
  // std::cout << "Triangulation lm construction took "
  //           << timer.time() << std::endl;

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

template <typename T_Traits>
bool Point_location_test<T_Traits>::attach_pl_strategies()
{
  typedef T_Traits Traits;

  // Initialize all point location strategies.
  CGAL::Timer timer;

  m_naive_pl->attach(*m_arr);
  m_simple_pl->attach(*m_arr);
  m_walk_pl->attach(*m_arr);

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

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

  // timer.reset(); timer.start();
  // m_location triangulation_lm_pl->attach(*m_arr);
  // timer.stop(); 
  // std::cout << "Triangulation lm construction took "
  //           << timer.time() << std::endl;

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
template <typename T_Traits>
bool Point_location_test<T_Traits>::perform()
{
  Objects_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES];
  typename Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
  typename Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
  typename Arrangement_2::Face_const_handle      fh_ref, fh_curr;
  
  // Locate the points in the list using all point location strategies.

  // std::cout << "Time in seconds" << std::endl;
  std::cout << std::endl;

  unsigned int pl_index = 0;

  query(*m_naive_pl, "Naive", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Naive

  query(*m_simple_pl, "Simple", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Simple
  
  query(*m_walk_pl, "Walk", m_query_points.begin(), m_query_points.end(),
        std::back_inserter(objs[pl_index++]));  // Walk

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

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
  
  // Triangulation
  // query(*m_triangulation_pl, "Triangulation",
  //       m_query_points.begin(), m_query_points.end(),
  //       std::back_inserter(objs[pl_index++]));  // Triangulation

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
  unsigned int pls_num = pl_index;
  std::cout << "Number of strategies is " << pls_num << std::endl;  

  // End Location

  int result = 0;

  //init all obejct iterators
  Object_iterator ob_iter[MAX_NUM_POINT_LOCATION_STRATEGIES];
  for (pl_index = 0; pl_index < pls_num; ++pl_index)
    ob_iter[pl_index] = objs[pl_index].begin();

  // get size of objects
  unsigned int size = objs[0].size();
  std::cout << "size is " << size << std::endl;

  for (pl_index = 0; pl_index < pls_num; ++pl_index) {
    if (size != objs[pl_index].size()) {
      std::cout << "Error: size of pl number " << pl_index << " is "
                << objs[pl_index].size() << std::endl;
      result = -1;
    }
  }

  // Assign and check results
  unsigned int qi; //qi is the query point index
  for (qi = 0; qi < size; ++qi) {
    // Assign object to a face
    if (CGAL::assign(fh_ref, ob_iter[0][qi])) {
      for (unsigned int pl_index = 1; pl_index < pls_num; ++pl_index) {
        if (! CGAL::assign(fh_curr, ob_iter[pl_index][qi])) {
          std::cout << "Error in point location number " << pl_index;
          if (CGAL::assign(fh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", an halfedge returned instead of a face"
                      << std::endl;
          }
          else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", a vertex returned instead of a face"
                      << std::endl;
          }
          else {
            std::cout << ", an unknowen object returned instead of a face"
                      << std::endl;
          }
          result = -1;
        }
        else if (fh_curr != fh_ref) {
          std::cout << "Error: point location number " 
                    << pl_index << " return a different face" << std::endl;
          result = -1;
        }
      }  
      //if (fh_ref->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
    }

    // Assign object to a halfedge
    else if (CGAL::assign (hh_ref, ob_iter[0][qi])) {
      std::cout << "Halfedge: " << hh_ref->curve() << std::endl;
      for (unsigned int pl_index = 1; pl_index < pls_num; ++pl_index) {
        if (! CGAL::assign(hh_curr, ob_iter[pl_index][qi])) {
          std::cout << "Error in point location number " << pl_index;
          if (CGAL::assign(fh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", a face returned instead of an halfedge"
                      << std::endl;
          }
          else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", a vertex returned instead of an halfedge"
                      << std::endl;
          }
          else {
            std::cout << ", an unknowen object returned instead of an halfedge"
                      << std::endl;
          }
          result = -1;
        }
        else if ((hh_curr != hh_ref) && (hh_curr->twin() != hh_ref)) {
          std::cout << "Error: point location number " 
                    << pl_index << " return a different halfedge" << std::endl;
          std::cout << "Halfedge (curr): " << hh_curr->curve() << std::endl;
          result = -1;
        }
      }
    }

    // Assign object to a vertex
    else if (CGAL::assign(vh_ref, ob_iter[0][qi])) {
      for (unsigned int pl_index = 1; pl_index < pls_num; ++pl_index) {
        if (! CGAL::assign(vh_curr, ob_iter[pl_index][qi])) {
          std::cout << "Error in point location number " << pl_index;
          if (CGAL::assign(fh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", a face returned instead of a vertex"<< std::endl;
          }
          else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi])) {
            std::cout << ", an halfedge returned instead of a vertex"
                      << std::endl;
          }
          else {
            std::cout << ", an unknown object returned instead of a vertex"
                      << std::endl;
          }
          result = -1;
        }
        else if (vh_curr != vh_ref) {
          std::cout << "Error: point location number " 
                    << pl_index << " return a different vertex"<< std::endl;
          result = -1;
        }
      }
      std::cout << "Vertex: "<< vh_ref->point() << std::endl;
    }

    else {
      std::cout << "Illegal point-location result." << std::endl;    
      result = -1;
    }
  }
  return (result == 0);
}

#endif
