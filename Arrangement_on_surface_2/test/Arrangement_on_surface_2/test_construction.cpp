#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP

  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Number_type;

#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <cstring> // for std::strncmp

#include "utils.h"

template <class Traits_>
class Point_equal {
public:
  typedef Traits_                      Traits_2;
  typedef typename Traits_2::Point_2   Point_2;

  Point_equal(const Traits_2& traits) : m_traits(traits) {}
  bool operator()(const Point_2& p1, const Point_2& p2)
  { return (m_traits.equal_2_object()(p1, p2)); }

private:
  const Traits_2& m_traits;
};

template <class Traits_>
class Curve_equal {
public:
  typedef Traits_                                  Traits_2;
  typedef typename Traits_2::X_monotone_curve_2    X_monotone_curve_2;

  Curve_equal(const Traits_2& traits) : m_traits(traits) {}
  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  {
    return (m_traits.equal_2_object()(c1, c2) && (c1.data() == c2.data()));
  }

private:
  const Traits_2& m_traits;
};

template <class Arrangement>
bool
are_same_results(Arrangement& arr, 
                 unsigned int n_v,
                 unsigned int n_e,
                 unsigned int n_f,
                 const std::vector<typename Arrangement::Point_2>&
                   pts_from_file,
                 const std::vector<typename Arrangement::X_monotone_curve_2>&
                   subcurves_from_file)
{
  typedef typename Arrangement::Traits_2              Traits_2;
  typedef typename Arrangement::Point_2               Point_2;
  typedef typename Arrangement::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Arrangement::Vertex_iterator       Vertex_iterator;
  typedef typename Arrangement::Edge_iterator         Edge_iterator;

  const Traits_2* traits = arr.traits();

  if (arr.number_of_vertices() != n_v) return false;
  if (arr.number_of_edges() != n_e) return false;
  if (arr.number_of_faces() != n_f) return false;

  std::vector<Point_2>  pts(arr.number_of_vertices());
  unsigned int i = 0;
  Vertex_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit, ++i)
    pts[i] = vit->point();
  std::sort(pts.begin(), pts.end(), Point_compare<Traits_2>());
  Point_equal<Traits_2>  point_eq(*traits);
  if (!std::equal(pts.begin(), pts.end(), pts_from_file.begin(), point_eq))
    return false;

  std::vector<X_monotone_curve_2> curves_res(arr.number_of_edges());
  i = 0;
  Edge_iterator eit;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i)
    curves_res[i] = eit->curve();    
  std::sort(curves_res.begin(), curves_res.end(), Curve_compare<Traits_2>());

  Curve_equal<Traits_2> curve_eq(*traits);
  if (! std::equal(curves_res.begin(), curves_res.end(),
                   subcurves_from_file.begin(), curve_eq))
    return false;

  return true;
}

typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Base_traits_2;
typedef Base_traits_2::Point_2                        Base_point_2;
typedef Base_traits_2::Curve_2                        Base_curve_2;
typedef Base_traits_2::X_monotone_curve_2             Base_x_monotone_curve_2;
typedef CGAL::Arr_curve_data_traits_2<Base_traits_2,
                                      unsigned int,
                                      std::plus<unsigned int> >  
                                                      Traits_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_iterator                Vertex_iterator;
typedef Arrangement_2::Edge_iterator                  Edge_iterator;
typedef std::vector<Curve_2>                          CurveContainer;

bool test_one_file(std::ifstream& in_file, bool verbose)
{ 
  // Read the input curves and isolated vertices.
  CurveContainer        curves;
  std::vector<Point_2>  iso_verts;

  unsigned int num_of_curves;
  in_file >> num_of_curves;
  unsigned int i;
  for (i = 0; i < num_of_curves; ++i) { 
    Base_point_2 source, target;
    in_file >> source >> target;
    
    Kernel ker;
    if (ker.compare_xy_2_object()(source, target) != CGAL::EQUAL) {
      Base_curve_2 base_cv(source, target);
      curves.push_back(Curve_2(base_cv, 1));
    }
    else iso_verts.push_back(source);
  }

  // Read the point and curves that correspond to the arrangement vertices
  // and edges, respectively.
  unsigned int n_vertices;
  unsigned int n_edges;
  unsigned int n_faces;

  in_file >> n_vertices >> n_edges >> n_faces;

  std::vector<Point_2>  pts_from_file(n_vertices);
  for (i = 0; i < n_vertices ; ++i)
    in_file >> pts_from_file[i];

  std::vector<X_monotone_curve_2>  subcurves_from_file(n_edges);
  for (i = 0; i < n_edges ; ++i) { 
    Base_curve_2 base_cv;
    unsigned int k;
    in_file >> base_cv >> k;
    subcurves_from_file[i] = Curve_2(base_cv, k);
  }

  Arrangement_2 arr;

  ////////////////////////////////////////////////////////////
  // test incremental construction
  CurveContainer::const_iterator it;
  for (it = curves.begin(); it != curves.end(); ++it)
    CGAL::insert(arr, *it);
  std::vector<Point_2>::const_iterator poit;
  for (poit = iso_verts.begin(); poit != iso_verts.end(); ++poit)
    CGAL::insert_point(arr, *poit);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cerr << "ERROR : The incremental insertion test failed." << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(1) Passed incremental insertion." << std::endl;
  arr.clear();

  ////////////////////////////////////////////////////////////
  // test aggregated construction

  CGAL::insert(arr, curves.begin(), curves.end());
  // when creating insert_points, this call should be fixed to insert_points.
  for (poit = iso_verts.begin(); poit != iso_verts.end(); ++poit)
    CGAL::insert_point(arr, *poit);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cerr << "ERROR : The aggregated construction test failed." << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(2) Passed aggregated construction." << std::endl;
  arr.clear();

  /////////////////////////////////////////////////////////////
  // insert half of the curves aggregatley and than insert the rest
  // aggregatley (test the addition visitor)
  CGAL::insert(arr, curves.begin(), curves.begin() + (num_of_curves/2));
  CGAL::insert(arr, curves.begin() + (num_of_curves/2), curves.end());
  // when creating insert_points, this call should be fixed to insert_points.
  for (poit = iso_verts.begin(); poit != iso_verts.end(); ++poit)
    CGAL::insert_point(arr, *poit);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cerr << "ERROR : The aggregated insertion test failed." << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(3) Passed aggregated insertion." << std::endl;
  arr.clear();

  //////////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves incrementally with  insert_x_monotone_curve
  for (i = 0; i < n_edges; ++i)
    CGAL::insert(arr, subcurves_from_file[i]);
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cerr << "ERROR : The incremental x-monotone test failed." << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(4) Passed incremental x-monotone insertion." << std::endl;
  arr.clear();

  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves aggregatley with  insert_x_monotone_curves
  CGAL::insert(arr, subcurves_from_file.begin(), subcurves_from_file.end());
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cerr << "ERROR : The aggregated x-monotone construction test failed." 
              << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(5) Passed aggregated x-monotone construction." << std::endl;
  arr.clear();


  /////////////////////////////////////////////////////////////////////
  // insert half of the disjoint subcurves aggregatley and than insert the
  // rest aggregatley with insert_x_monotone_curves(test the addition visitor)
  CGAL::insert(arr, subcurves_from_file.begin(),
               subcurves_from_file.begin() + (n_edges/2));
  CGAL::insert(arr, subcurves_from_file.begin() + (n_edges/2),
               subcurves_from_file.end());
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout << "ERROR : The aggregated x-monotone inertion test failed." 
              << std::endl;
     return false;
  }

  if (verbose)
    std::cout << "(6) Passed aggregated x-monotone insertion." << std::endl;
  arr.clear();

  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves incrementally with 
  // insert_non_intersecting_curve
  for (i = 0; i < n_edges; ++i)
    CGAL::insert_non_intersecting_curve(arr, subcurves_from_file[i]);
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout << "ERROR : The incremental non-intersecting test failed." 
              << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(7) Passed incremental non-intersecting insertion." 
              << std::endl;
  arr.clear();

  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves aggregatley with
  // insert_non_intersecting_curves
  CGAL::insert_non_intersecting_curves(arr, subcurves_from_file.begin(),
                                       subcurves_from_file.end());
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout 
      <<"ERROR : The aggregated non-intersecting construction test failed." 
      << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(8) Passed aggregated non-intersecting construction." 
              << std::endl;
  arr.clear();

  // insert half of the disjoint subcurves aggregatley and than insert the
  // rest aggregatley with insert_non_intersecting_curves (test the addition
  // visitor)
  CGAL::insert_non_intersecting_curves(arr, subcurves_from_file.begin(),
                                       subcurves_from_file.begin() +
                                         (n_edges/2));
  CGAL::insert_non_intersecting_curves(arr,
                                       subcurves_from_file.begin() +
                                         (n_edges/2),
                                       subcurves_from_file.end());
  for (i = 0; i < iso_verts.size(); ++i)
    CGAL::insert_point(arr, iso_verts[i]);
  if (! are_same_results(arr, n_vertices, n_edges, n_faces,
                         pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout 
      <<"ERROR : The aggregated non-intersecting insertion test failed." 
      << std::endl;
    return false;
  }

  if (verbose)
    std::cout << "(9) Passed aggregated non-intersecting insertion." 
              << std::endl;
  arr.clear();

  return true;
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    std::exit(-1);
  }

  bool  verbose = false;
  
  if (argc > 2 && std::strncmp(argv[2], "-v", 2) == 0)
    verbose = true;

  int success = 0;
  for (int i = 1; i < argc; ++i) {
    std::string str(argv[i]);
    if (str.empty()) continue;

    std::string::iterator itr = str.end();
    --itr;
    while (itr != str.begin()) {
      std::string::iterator tmp = itr;
      --tmp;
      if (*itr == 't')  break;
      
      str.erase(itr);
      itr = tmp;
    }
    if (str.size() <= 1) continue;
    std::ifstream inp(str.c_str());
    if (!inp.is_open()) {
      std::cerr << "Failed to open " << str << std::endl;
      return (-1);
    }
    if (! test_one_file(inp, verbose)) {
      inp.close();
      std::cerr << str << ": ERROR" << std::endl;
      success = -1;
    }
    else std::cout <<str << ": succeeded" << std::endl;
    inp.close();
  }
  
  return success;
}
