#include <vector>
#include <cstdlib>
#include <sstream>
#include  <string>

#include <CGAL/basic.h>

#include "test_configuration.h"
#include "test_kernel.h"

#include <CGAL/Cartesian.h>
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

#if TEST_TRAITS == SEGMENT_TRAITS
#include <CGAL/Arr_segment_traits_2.h>
#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS
#include <CGAL/Arr_circle_segment_traits_2.h>
#elif TEST_TRAITS == LINEAR_TRAITS
#include <CGAL/Arr_linear_traits_2.h>
#endif

#if TEST_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>          Traits;
#define TRAITS_TYPE "point_location_Segments"

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS
typedef CGAL::Arr_circle_segment_traits_2<Kernel>   Traits;
typedef CGAL::Cartesian<Number_type>                Rat_kernel;
typedef Rat_kernel::FT                              Rat_nt;
typedef Rat_kernel::Circle_2                        Circle_2;
typedef Rat_kernel::Line_2                          Line_2;
typedef Rat_kernel::Segment_2                       Segment_2;
typedef Rat_kernel::Point_2                         Rat_point_2;
#define TRAITS_TYPE "point_location_circle Segments"

#elif TEST_TRAITS == LINEAR_TRAITS
typedef CGAL::Arr_linear_traits_2<Kernel>           Traits;
#endif

typedef Traits::Point_2                             Point_2;
typedef Traits::Curve_2                             Curve_2;
typedef Traits::X_monotone_curve_2                  X_monotone_curve_2;

#include "IO_test.h"

typedef CGAL::Arrangement_2<Traits>                 Arrangement_2;
typedef Arrangement_2::Halfedge_handle              Halfedge_handle;
typedef Arrangement_2::Edge_const_iterator          Edge_const_iterator;
typedef Arrangement_2::Vertex_const_iterator        Vertex_const_iterator;

typedef CGAL::Arr_naive_point_location<Arrangement_2>     
                                                    Naive_point_location;
typedef CGAL::Arr_simple_point_location<Arrangement_2>     
                                                    Simple_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
                                                    Walk_point_location;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2> 
                                                    Lm_point_location;
typedef CGAL::Arr_random_landmarks_generator<Arrangement_2>
                                                    Random_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Random_lm_generator> 
                                                    Lm_random_point_location;
typedef CGAL::Arr_grid_landmarks_generator<Arrangement_2>
                                                    Grid_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Grid_lm_generator> 
                                                    Lm_grid_point_location;
typedef CGAL::Arr_halton_landmarks_generator<Arrangement_2>
                                                    Halton_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Halton_lm_generator> 
                                                    Lm_halton_point_location;
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
                                                    Middle_edges_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2,
                                           Middle_edges_generator> 
  Lm_middle_edges_point_location;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> 
  Trapezoid_ric_point_location;

typedef CGAL::Arr_landmarks_specified_points_generator<Arrangement_2>
  Specified_points_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2,
                                           Specified_points_generator> 
  Lm_specified_points_point_location;

// typedef CGAL::Arr_triangulation_point_location<Arrangement_2> 
//   Lm_triangulation_point_location;

// ===> Add new point location type here <===

typedef std::vector<Point_2>                              Points_vector;
typedef Points_vector::iterator                           Point_iterator;
typedef std::vector<X_monotone_curve_2>                   Xcurves_vector;
typedef std::vector<Curve_2>                              Curves_vector;
typedef std::vector<CGAL::Object>                         Objects_vector;
typedef Objects_vector::iterator                          Object_iterator;

// ===> Change the number of point-location startegies
//      when a new point location is added. <===
#define MAX_NUM_POINT_LOCATION_STRATEGIES 11

int remove(Arrangement_2& arr, const X_monotone_curve_2& xcv)
{
  int rc = -1;          // be pasimistic, assume nothing is removed.
  
  const Traits* traits = arr.geometry_traits();
  Traits::Equal_2 equal = traits->equal_2_object();

  Arrangement_2::Edge_iterator eit;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    const X_monotone_curve_2& xcv_arr = eit->curve();
    if (equal(xcv, xcv_arr)) {
      arr.remove_edge(eit);
      rc = 0;           // found a curve to remove.
      break;
    }
  }
  return rc;
}

int read_perform_opts(Arrangement_2& arr,
                      const Xcurves_vector& xcurves, 
                      const Points_vector& points, 
                      const Curves_vector& curves, 
                      const char* ops_filename)
{
  int rc = 0;

  CGAL::Timer timer;
  timer.reset(); 
  timer.start();

  std::ifstream op_stream(ops_filename);
  if (!op_stream.is_open()) {
    std::cerr << "Cannot open file " << ops_filename << "!" << std::endl;
    return -1;
  }
  std::string sline;
  while (std::getline(op_stream, sline)) {
    std::istringstream line(sline);
    char cmd;
    line >> cmd;

    if (cmd == 'a') {
      // Insert all into the arrangement
      insert(arr, xcurves.begin(), xcurves.end());
      // insert(arr, points.begin(), points.end());
      insert(arr, curves.begin(), curves.end());
      continue;
    }
    
    unsigned int id;
    line >> id;
    if (id >= xcurves.size()) {
      std::cerr << "Index of x-monotone curve " << id << " is out of range ("
                << xcurves.size() << ") in " << ops_filename << "!"
                << std::endl;
      rc = -1;
      continue;
    }
    if (cmd == 'i') CGAL::insert(arr, xcurves[id]);

    if (cmd == 'd') {
      if (remove(arr, xcurves[id]) < 0)
        rc = -1;
    }
  }
  op_stream.close();
  timer.stop(); ///END
  std::cout << "Arrangement aggregate construction took " 
            << timer.time() << std::endl;  

  return rc;
}


/*! */
int construct_and_query(Arrangement_2& arr,
                        const Xcurves_vector& xcurves, 
                        const Points_vector& points, 
                        const Curves_vector& curves, 
                        const char* ops_filename,  
                        Points_vector& query_points)
{
  //init - all point locations
  CGAL::Timer timer;

  Naive_point_location naive_pl(arr);                           // 0
  Simple_point_location simple_pl(arr);                         // 1
  Walk_point_location walk_pl(arr);                             // 2

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

  timer.reset(); timer.start();
  Lm_point_location lm_pl(arr);                                 // 3
  timer.stop(); 
  std::cout << "Lm (vert) construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  Random_lm_generator random_g(arr);
  Lm_random_point_location random_lm_pl(arr, &random_g);        // 4
  timer.stop(); 
  std::cout << "Random lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  Grid_lm_generator grid_g(arr);      
  Lm_grid_point_location grid_lm_pl(arr, &grid_g);              // 5
  timer.stop(); 
  std::cout << "Grid lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  Halton_lm_generator halton_g(arr);
  Lm_halton_point_location halton_lm_pl(arr, &halton_g);        // 6
  timer.stop(); 
  std::cout << "Halton lm construction took " << timer.time() << std::endl;

  timer.reset(); timer.start();
  Middle_edges_generator middle_edges_g(arr);
  Lm_middle_edges_point_location middle_edges_lm_pl(arr, &middle_edges_g); // 7
  timer.stop(); 
  std::cout << "Middle edges lm construction took " << timer.time()
            << std::endl;

  Specified_points_generator::Points_set spc_points;
  //spc_points.push_back(Point_2(1, 1));
  //spc_points.push_back(Point_2(2, 2));
  timer.reset(); timer.start();
  //Specified_points_generator specified_points_g(arr,spc_points);
  Specified_points_generator specified_points_g(arr);
  Lm_specified_points_point_location
    specified_points_lm_pl(arr, &specified_points_g);                   // 8
  timer.stop(); 
  std::cout << "Specified_points lm construction took "
            << timer.time() << std::endl;

  // timer.reset(); timer.start();
  // Lm_triangulation_point_location triangulation_lm_pl(arr);          // 9
  // timer.stop(); 
  // std::cout << "Triangulation lm construction took "
  //           << timer.time() << std::endl;

#endif
  
  timer.reset(); timer.start();
  Trapezoid_ric_point_location trapezoid_ric_pl_grnt(arr);              // 9
  timer.stop(); 
  std::cout << "Trapezoid RIC with-guarantees construction took " 
            << timer.time() << std::endl;
  
  timer.reset(); timer.start();
  Trapezoid_ric_point_location trapezoid_ric_pl_no_grnt(arr,false);     // 10
  timer.stop(); 
  std::cout << "Trapezoid RIC without-guarantees construction took " 
            << timer.time() << std::endl;
  
  std::cout << std::endl;

  // ===> Add new point location instance here. <===

  if (read_perform_opts(arr, xcurves, points, curves, ops_filename) < 0)
    return -1;
 
  // Print the size of the arrangement.
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  Objects_vector objs[MAX_NUM_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
  Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
  Arrangement_2::Face_const_handle      fh_ref, fh_curr;
  
  // LOCATE the points in the list using all PL strategies

  // std::cout << "Time in seconds" << std::endl;
  std::cout << std::endl;

  Point_iterator piter;
  unsigned int pl_index = 0;

  // Naive
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = naive_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Naive location took " << timer.time() << std::endl;

  // Simple
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = simple_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Simple location took " << timer.time() << std::endl;

  // Walk
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = walk_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Walk location took " << timer.time() << std::endl;

#if (TEST_TRAITS == SEGMENT_TRAITS) || (TEST_TRAITS == LINEAR_TRAITS)

  // Landmarks (vertices)
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Landmarks (vertices) location took " << timer.time()
            << std::endl;

  // Landmarks random
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = random_lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Random LM location took " << timer.time() << std::endl;

  // Landmarks grid
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = grid_lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Grid LM location took " << timer.time() << std::endl;

  // Landmarks Halton
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = halton_lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Halton LM location took " << timer.time() << std::endl;

  // Landmarks Middle edges
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = middle_edges_lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Middle edges LM location took " << timer.time() <<std::endl;

  // Landmarks specified points
  timer.reset(); timer.start();
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = specified_points_lm_pl.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop();
  ++pl_index;
  std::cout << "Specified points LM location took " << timer.time()
            << std::endl;

  // Triangulation
  // timer.reset(); timer.start();
  // for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
  //   Point_2 q = (*piter);
  //   CGAL::Object obj = triangulation_lm_pl.locate(q);
  //   objs[pl_index].push_back(obj);
  // }
  // timer.stop();
  // ++pl_index;
  // std::cout << "Triangulation LM location took " << timer.time() << std::endl;

#endif

  // Trapezoidal RIC
  timer.reset(); timer.start(); //START
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = trapezoid_ric_pl_grnt.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop(); ///END
  ++pl_index;
  std::cout << "Trapezoidal RIC with-guarantees location took " << timer.time() << std::endl;

  // Trapezoidal RIC without guarantees
  timer.reset(); timer.start(); //START
  for (piter = query_points.begin(); piter != query_points.end(); ++piter) {
    Point_2 q = (*piter);
    CGAL::Object obj = trapezoid_ric_pl_no_grnt.locate(q);
    objs[pl_index].push_back(obj);
  }
  timer.stop(); ///END
  ++pl_index;
  std::cout << "Trapezoidal RIC without-guarantees location took " << timer.time() << std::endl;
    
  std::cout << std::endl;

  // ===> Add a call to operate the the new point location. <===

  //END LOCATION 
  unsigned int pls_num = pl_index;
  std::cout << "Number of strategies is " << pl_index << std::endl;
  int result = 0;

  //init all obejct iterators
  Object_iterator ob_iter[MAX_NUM_POINT_LOCATION_STRATEGIES];
  for (pl_index = 0; pl_index < pls_num; ++pl_index)
    ob_iter[pl_index] = objs[pl_index].begin();

  //get size of objects
  unsigned int size = objs[0].size();
  std::cout << "size is " << size << std::endl;

  for (pl_index = 0; pl_index < pls_num; ++pl_index) {
    if (size != objs[pl_index].size()) {
      std::cout << "Error: size of pl number " << pl_index << " is "
                << objs[pl_index].size() << std::endl;
      result = -1;
    }
  }

  //assign and check results
  unsigned int qi; //qi is the query point index
  for (qi = 0; qi < size; ++qi) {
    //assign object to a face
    if (CGAL::assign(fh_ref, ob_iter[0][qi])) {
      for (int pl_index = 1; pl_index < pls_num; ++pl_index) {
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

    //assign object to a halfedge
    else if (CGAL::assign (hh_ref, ob_iter[0][qi])) {
      std::cout << "Halfedge: " << hh_ref->curve() << std::endl;
      for (int pl_index = 1; pl_index < pls_num; ++pl_index) {
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

    //assign object to a vertex
    else if (CGAL::assign(vh_ref, ob_iter[0][qi])) {
      for (int pl_index = 1; pl_index < pls_num; ++pl_index) {
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
  return (result);
}

/*! */
int read_points(const char* points_filename, Points_vector& points)
{
  // read points from file into associative container
  std::ifstream point_stream(points_filename);
  if (!point_stream.is_open()) {
    std::cerr << "Cannot open file " << points_filename << "!" << std::endl;
    return (-1);
  }

  int points_count = 0;
  point_stream >> points_count;
  points.resize(points_count);
  IO_test<Traits> io_test;
  for (int i = 0; i < points_count; ++i)
    io_test.read_point(point_stream, points[i]);

  return 0;
}

/*! */
int read_input(const char* curves_filename,
               Xcurves_vector& xcurves, Points_vector& points,
               Curves_vector& curves)
{
  // read curves from file into associative container
  std::ifstream in_stream(curves_filename);
  if (!in_stream.is_open()) {
    std::cerr << "Cannot open file " << curves_filename << "!" << std::endl;
    return (-1);
  }

  IO_test<Traits> io_test;

  unsigned int xcurves_count = 0;
  in_stream >> xcurves_count;
  xcurves.resize(xcurves_count);
  for (int i = 0; i < xcurves_count; ++i)
    if (!io_test.read_xcurve(in_stream, xcurves[i]))
      return (-1);

  unsigned int points_count = 0;
  in_stream >> points_count;
  points.resize(points_count);
  for (int i = 0; i < points_count; ++i)
    if (!io_test.read_point(in_stream, points[i]))
      return (-1);
  
  unsigned int curves_count = 0;
  in_stream >> curves_count;
  curves.resize(curves_count);
  for (int i = 0; i < curves_count; ++i)
    if (!io_test.read_curve(in_stream, curves[i]))
      return (-1);

  return 0;
}

int check_point_location(Arrangement_2& arr,
                         const Xcurves_vector& xcurves, 
                         const Points_vector& points, 
                         const Curves_vector& curves, 
                         const char* ops_filename, 
                         const char* queries_filename)
{
  // Read point and insert them into a list of points
  Points_vector query_points;
  if (read_points(queries_filename, query_points) < 0) {
    std::cerr << "ERROR in read_points."<< std::endl << std::endl;
    return -1;
  }

  // Check point location of points
  if (construct_and_query(arr, xcurves, points, curves, 
                          ops_filename, query_points) < 0) 
  {
    std::cerr << "ERROR in check_point_location."<< std::endl << std::endl;
    return -1;
  }
  std::cout << std::endl;
  return 0;
}

// int remove(Arrangement_2& arr, const X_monotone_curve_2& xcv)
// {
//   int rc = -1;          // be pasimistic, assume nothing is removed.
  
//   const Traits* traits = arr.geometry_traits();
//   Traits::Equal_2 equal = traits->equal_2_object();

//   Arrangement_2::Edge_iterator eit;
//   for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
//     const X_monotone_curve_2& xcv_arr = eit->curve();
//     if (equal(xcv, xcv_arr)) {
//       arr.remove_edge(eit);
//       rc = 0;           // found a curve to remove.
//     }
//   }
//   return rc;
// }

// int read_perform_opts(Arrangement_2& arr,
//                       const Xcurves_vector& xcurves, 
//                       const Points_vector& points, 
//                       const Curves_vector& curves, 
//                       const char* ops_filename)
// {
//   int rc = 0;

//   CGAL::Timer timer;
//   timer.reset(); 
//   timer.start();

//   std::ifstream op_stream(ops_filename);
//   if (!op_stream.is_open()) {
//     std::cerr << "Cannot open file " << ops_filename << "!" << std::endl;
//     return -1;
//   }
//   std::string sline;
//   while (std::getline(op_stream, sline)) {
//     std::istringstream line(sline);
//     char cmd;
//     line >> cmd;

//     if (cmd == 'a') {
//       // Insert all into the arrangement
//       insert(arr, xcurves.begin(), xcurves.end());
//       // insert(arr, points.begin(), points.end());
//       insert(arr, curves.begin(), curves.end());
//       continue;
//     }
    
//     unsigned int id;
//     line >> id;
//     if (id >= xcurves.size()) {
//       std::cerr << "Index of x-monotone curve " << id << " is out of range ("
//                 << xcurves.size() << ") in " << ops_filename << "!"
//                 << std::endl;
//       rc = -1;
//       continue;
//     }
//     if (cmd == 'i') CGAL::insert(arr, xcurves[id]);

//     if (cmd == 'd') {
//       if (remove(arr, xcurves[id]) < 0)
//         rc = -1;
//     }
//   }
//   op_stream.close();
//   timer.stop(); ///END
//   std::cout << "Arrangement aggregate construction took " 
//             << timer.time() << std::endl;  

//   return rc;
// }

int test(const char* curves_filename, const char* ops_filename,
         const char* queries_filename)
{
  // Read curves 
  Xcurves_vector xcurves;
  Curves_vector curves;
  Points_vector points;
  if (read_input(curves_filename, xcurves, points, curves) < 0)
    return -1;

  // Read and perform operations  
  Arrangement_2 arr;
  //if (read_perform_opts(arr, xcurves, points, curves, ops_filename) < 0)
  //  return -1;
  
  // Issue point location queries.
  if (check_point_location(arr,xcurves, points, curves, 
                           ops_filename, queries_filename) < 0)
  {
    return -1;
  }
  return 0;
}

int main(int argc, char* argv[])
{
   // Obtain arguments
  std::cout << "argc: " << argc << std::endl;
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " input_file op_file query_file"
              << std::endl;
    std::cout << "input_file  - the input curves file" << std::endl;
    std::cout << "op_file     - the input operations points" << std::endl;
    std::cout << "query_file  - the input query points" << std::endl;
    return -1;
  }

  int success = 0;
  for (int i = 1; i < argc; i += 3) {
    const char* curves_filename = argv[i];
    const char* ops_filename = argv[i+1];
    const char* queries_filename = argv[i+2];
    if (test(curves_filename, ops_filename, queries_filename) < 0) {
      std::cerr << "ERROR : " << argv[0] << " " << argv[i] << " "
                << argv[i+1] << " " << argv[i+2] << std::endl;
      success = -1;
    }
  }
  return success;
}
