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
#include<CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cstdlib>

#include "utils.h"

#include <CGAL/assertions.h>

CGAL::Failure_function prev_error_handler;
CGAL::Failure_function prev_warning_handler;

void failure_handler(const char * type, const char * expr, const char * file,
                        int line, const char * msg)
{
    std::cout << "type " << type << std::endl;
    std::cout << "expr " << expr << std::endl;
    std::cout << "file " << file << std::endl;
    std::cout << "line " << line << std::endl;
}

template <class Traits>
class Point_equal
{
  typedef typename Traits::Point_2    Point_2;
public:

  bool operator()(const Point_2& p1, const Point_2& p2)
  {
    Traits tr;
    return (tr.equal_2_object()(p1, p2));
  }
};

template <class Traits>
class Curve_equal
{
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
public:

  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  {
    Traits tr;
    return (tr.equal_2_object()(c1, c2) && (c1.data() == c2.data()));
  }
};

template <class Arrangement>
bool are_same_results
    (Arrangement& arr, 
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

  if(arr.number_of_vertices() != n_v)
    return false;

  if(arr.number_of_edges() != n_e)
    return false;

  if(arr.number_of_faces() != n_f)
    return false;

  
  std::vector<Point_2>  pts(arr.number_of_vertices());
  unsigned int i=0;
  for(Vertex_iterator vit = arr.vertices_begin();
      vit != arr.vertices_end();
      ++vit, ++i)
  {
    pts[i] = vit->point();
  }
  std::sort(pts.begin(), pts.end(), Point_compare<Traits_2>());
  Point_equal<Traits_2>  point_eq;
  if(!std::equal(pts.begin(), pts.end(), pts_from_file.begin(), point_eq))
    return false;

  std::vector<X_monotone_curve_2> curves_res(arr.number_of_edges());
  i=0;
  for(Edge_iterator eit = arr.edges_begin();
      eit != arr.edges_end();
      ++eit, ++i)
  {
    curves_res[i] = eit->curve();    
  }
  std::sort(curves_res.begin(), curves_res.end(), Curve_compare<Traits_2>());

  Curve_equal<Traits_2> curve_eq;
  if (! std::equal (curves_res.begin(), curves_res.end(),
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



bool  test_one_file(std::ifstream& in_file)
{ 
  CurveContainer curves;
  std::vector<Point_2>  iso_verts;

  unsigned int num_of_curves;
  in_file >> num_of_curves;
  unsigned int i;
  for(i=0 ; i < num_of_curves ; i++)
  { 
    Base_point_2 source, target;
    in_file >> source >> target;
    
    Kernel ker;
    if (ker.compare_xy_2_object() (source, target) != CGAL::EQUAL)
    {
      Base_curve_2 base_cv(source, target);
      curves.push_back(Curve_2(base_cv, 1));
    }
    else
    {
      iso_verts.push_back(source);
    }
  }

  unsigned int n_vertices;
  unsigned int n_edges;
  unsigned int n_faces;

  in_file >> n_vertices >> n_edges >> n_faces;
  std::vector<Point_2>  pts_from_file(n_vertices);
  for(i=0 ; i < n_vertices ; i++)
  { 
    in_file >> pts_from_file[i];
  }

  std::vector<X_monotone_curve_2>  subcurves_from_file(n_edges);
  for(i=0 ; i < n_edges ; i++)
  { 
    Base_curve_2 base_cv;
    unsigned int k;
    in_file >> base_cv >> k;
    subcurves_from_file[i] = Curve_2(base_cv, k);
  }



  Arrangement_2 arr;

  // test incremental construction
  for(CurveContainer::const_iterator it = curves.begin(); 
      it != curves.end(); ++it)
  {
    CGAL::insert_curve(arr, *it);
  }

  std::vector<Point_2>::const_iterator poit;
  for (poit = iso_verts.begin();
       poit != iso_verts.end(); ++poit)
  {
    CGAL::insert_point(arr, *poit);
  }

  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_curve failed\n";
    return false;
  }

  arr.clear();

  ////////////////////////////////////////////////////////////
  // test aggregate construction

  CGAL::insert_curves(arr, curves.begin(), curves.end());
  // when creating insert_points, this call should be fixed to insert_points.
  for (poit = iso_verts.begin();
       poit != iso_verts.end(); ++poit)
  {
    CGAL::insert_point(arr, *poit);
  }

  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_curves failed\n";
    return false;
  }
  arr.clear();

  /////////////////////////////////////////////////////////////
  // insert half of the curves aggregatley and than insert the rest
  // aggregatley (test the addition visitor)
  CGAL::insert_curves(arr, curves.begin(), curves.begin() + (curves.size()/2));
  CGAL::insert_curves(arr, curves.begin() + (curves.size()/2), curves.end());
  // when creating insert_points, this call should be fixed to insert_points.
  for (poit = iso_verts.begin();
       poit != iso_verts.end(); ++poit)
  {
    CGAL::insert_point(arr, *poit);
  }

  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_curves (addition) failed\n";
    return false;
  }
  arr.clear();

  //////////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves incrementally with  insert_x_monotone_curve

  for(i=0; i<n_edges; ++i)
  {
    CGAL::insert_x_monotone_curve(arr, subcurves_from_file[i]);
  }
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }
  
  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_x_monotone_curve failed\n";
    return false;
  }
  arr.clear();
  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves aggregatley with  insert_x_monotone_curves
  CGAL::insert_x_monotone_curves (arr,
                                  subcurves_from_file.begin(),
                                  subcurves_from_file.end());
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }
  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_x_monotone_curves failed\n";
    return false;
  }
  arr.clear();


  /////////////////////////////////////////////////////////////////////
  // insert half of the disjoint subcurves aggregatley and than insert the
  // rest aggregatley with insert_x_monotone_curves(test the addition visitor)
  CGAL::insert_x_monotone_curves (arr,
                                  subcurves_from_file.begin(),
                                  subcurves_from_file.begin() + (n_edges/2));
  CGAL::insert_x_monotone_curves (arr,
                                  subcurves_from_file.begin() + (n_edges/2),
                                  subcurves_from_file.end());
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }

  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_x_monotone_curves (addition) failed\n";
    return false;
  }
  arr.clear();

  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves incrementally with 
  // insert_non_intersecting_curve
  for(i=0; i<n_edges; ++i)
  {
    CGAL::insert_non_intersecting_curve(arr, subcurves_from_file[i]);
  }
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }
  
  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_non_intersecting_curve failed\n";
    return false;
  }
  arr.clear();

  /////////////////////////////////////////////////////////////////////
  // insert the disjoint subcurves aggregatley with
  // insert_non_intersecting_curves
  CGAL::insert_non_intersecting_curves (arr,
                                        subcurves_from_file.begin(),
                                        subcurves_from_file.end());
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }
  
  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_non_intersecting_curves failed\n";
    return false;
  }
  arr.clear();

  

  // insert half of the disjoint subcurves aggregatley and than insert the
  // rest aggregatley with insert_non_intersecting_curves (test the addition
  // visitor)
  CGAL::insert_non_intersecting_curves 
    (arr,
     subcurves_from_file.begin(),
     subcurves_from_file.begin() + (n_edges/2));
  
  CGAL::insert_non_intersecting_curves
    (arr,
     subcurves_from_file.begin() + (n_edges/2),
     subcurves_from_file.end());
  
  for(i=0; i<iso_verts.size(); ++i)
  {
    CGAL::insert_point(arr, iso_verts[i]);
  }

  if (! are_same_results (arr,
                          n_vertices, n_edges, n_faces,
                          pts_from_file, subcurves_from_file) ||
      ! CGAL::is_valid(arr))
  {
    std::cout<<"ERROR : insert_non_intersecting_curves (addition) failed\n";
    return false;
  }
  arr.clear();

  return true;
}

int main(int argc, char **argv)
{
  CGAL::set_error_behaviour(CGAL::CONTINUE);
  CGAL::set_warning_behaviour(CGAL::CONTINUE);
  prev_error_handler = CGAL::set_error_handler(failure_handler);
  prev_warning_handler = CGAL::set_warning_handler(failure_handler);
  if(argc < 2)
  {
    std::cerr<<"Missing input file\n";
    std::exit (-1);
  }

  int success = 0;
  for(int i=1; i<argc; ++i)
  {
    std::string str(argv[i]);
    if(str.empty())
      continue;

    std::string::iterator itr = str.end();
    --itr;
    while(itr != str.begin())
    {
      std::string::iterator tmp = itr;
      --tmp;
      if(*itr == 't')
        break;
      
      str.erase(itr);
      itr = tmp;
      
    }
    if(str.size() <= 1)
      continue;
    std::ifstream inp(str.c_str());
    if(!inp.is_open())
    {
      std::cerr<<"Failed to open " <<str<<"\n";
      return (-1);
    }
    try
    {
      if (! test_one_file(inp))
      {
        inp.close();
        std::cout<<str<<": ERROR\n";
        success = -1;
      }
      else
      {
        std::cout<<str<<": succeeded\n";
      }
    }
    catch (std::bad_cast e)
    {
      inp.close();
      std::cout<<str<<": ERROR problem with casting\n";
      std::cout<< e.what() << std::endl;
      success = -1;			    
    }
    catch (std::bad_alloc e)
    {
      inp.close();
      std::cout<<str<<": ERROR problem with memory allocation\n";
      std::cout<< e.what() << std::endl;
      success = -1;
    }
    catch (std::bad_typeid e)
    {
      inp.close();
      std::cout<<str<<": ERROR problem with typeid\n";
      std::cout<< e.what() << std::endl;
      success = -1;
    }
    catch (std::runtime_error e)
    {
      inp.close();
      std::cout<<str<<": ERROR runtime error\n";
      std::cout<< e.what() << std::endl;
      success = -1;
    }
    catch (std::exception e)
    {
      inp.close();
      std::cout<<str<<": ERROR file exists but general exception is thrown\n";
      std::cout<< e.what() << std::endl;
      success = -1;
    }
    inp.close();
  }
  CGAL::set_error_handler(prev_error_handler);
  CGAL::set_warning_handler(prev_warning_handler);
  return success;
}
