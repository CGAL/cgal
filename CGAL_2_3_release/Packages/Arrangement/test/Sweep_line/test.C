#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define Cartesian Cr
#define Arr_polyline_traits APT
#define Quotient Qt
#define Planar_ma2p_traits_wrap PMTW
#define _List_iterator Lit
#define bidirectional_iterator_tag Bitt
#define Planar_map_2 Pm2
#define Arr_2_default_dcel A2dd
#define Point_2 Pt2
#define allocator Altr
#define Td_X_trapezoid TdXt
#define Td_traits Tdt
#endif

#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#ifndef CGAL_PM_DEFAULT_DCEL_H
#include <CGAL/Pm_default_dcel.h>
#endif
#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif
#ifndef CGAL_ARR_PMWX_H
#include <CGAL/Arr_pmwx.h>
#endif
#ifndef SWEEP_TO_CONSTRUCT_PLANAR_MAP_H
#include <CGAL/sweep_to_construct_planar_map.h>
#endif
#ifndef SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H
#include <CGAL/sweep_to_produce_planar_map_subcurves.h>
#endif
//#ifndef CGAL_SWEEP_LINE_H
//#include <CGAL/Sweep_line.h>
//#endif

#define CGAL_SEGMENT_TRAITS        1
#define CGAL_SEGMENT_LEDA_TRAITS   2
#define CGAL_POLYLINE_TRAITS      11
#define CGAL_POLYLINE_LEDA_TRAITS 12

// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_LEDA_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS)


int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

// Choose traits

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
#include <CGAL/Arr_segment_exact_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
#error Currently not supported (July 2000)
#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_polyline_traits.h>
#else
#error No traits defined for test
#endif

// Picking a default  point location strategy
// See comment above.
#ifndef CGAL_ARR_TEST_POINT_LOCATION
//#define CGAL_ARR_TEST_POINT_LOCATION 1
#define CGAL_ARR_TEST_POINT_LOCATION 2
//#define CGAL_ARR_TEST_POINT_LOCATION 3
#endif

#if CGAL_ARR_TEST_POINT_LOCATION == 1
  // By default we use Trapezoidal Decomposition
#elif CGAL_ARR_TEST_POINT_LOCATION == 2
#include <CGAL/Pm_naive_point_location.h>  
#elif CGAL_ARR_TEST_POINT_LOCATION == 3
#include <CGAL/Pm_walk_along_line_point_location.h>
#else
#error No point location strategy defined for test
#endif
 
// Using my own temporary version which include arr.is_valid()
//#include "Arrangement_2_Debug.h"
//#include "Arrangement_2_Shai.h"

#ifndef CGAL_PM_DEFAULT_DCEL_H
#include <CGAL/Pm_default_dcel.h>
#endif
#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

//#ifndef CGAL_ARR_PMWX_H
//#include <CGAL/Arr_pmwx.h>
//#endif

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>

//#define  CGAL_PMWX_TEST_SWEEP 

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
  typedef CGAL::Quotient<int>                  NT;
//typedef leda_rational                        NT;
  typedef CGAL::Cartesian<NT>                  R;
  typedef CGAL::Arr_segment_exact_traits<R>    Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                        NT;
  typedef CGAL::Arr_leda_segment_exact_traits  Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  typedef CGAL::Quotient<int>                  NT;
  typedef CGAL::Cartesian<NT>                  R;
  typedef CGAL::Arr_polyline_traits<R>         Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
  typedef leda_rational                        NT;
  typedef CGAL::Arr_leda_polyline_traits<>     Traits;

#endif

typedef Traits::Point                          Point;
typedef Traits::X_curve                        X_curve;
typedef Traits::Curve                          Curve;
typedef CGAL::Pm_default_dcel<Traits>          Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>       PM;
typedef CGAL::Arr_pmwx<PM>                     Pmwx;
typedef Pmwx::Pmwx_change_notification         Notifier;
 
// we use the namespace std for compatability with MSVC
typedef std::list<Point>                     Point_list;

// Defining IO operators for polyline curves.
#if (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS)


CGAL_BEGIN_NAMESPACE

std::ostream&  operator<<(std::ostream& os,  
			  const Curve& cv)
{
  typedef Curve::const_iterator       Points_iterator;
  
  os<<cv.size()<<std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os<<" "<<*points_iter;

  return os;
}


std::istream&  operator>>(std::istream& in,  
			  Curve& cv)
{
  typedef Curve::value_type           Point;

  std::size_t  size;

  in >> size;

  for (unsigned int i = 0; i < size; i++){
    Point  p;
    
    in >> p;
    
    cv.push_back(p);  
  }
  
  return in;
}

CGAL_END_NAMESPACE
#endif

using namespace std;


template <class PM>
class Pm_polyline_traits_test
{
  PM  pm;  

public:
#if CGAL_ARR_TEST_POINT_LOCATION == 3  
  Pm_polyline_traits_test() : pm(new CGAL::Pm_walk_along_line_point_location<PM>) {};

#elif CGAL_ARR_TEST_POINT_LOCATION == 2
  Pm_polyline_traits_test() : pm(new CGAL::Pm_naive_point_location<PM>) {};
#else
  // None
#endif     

  /****************************
   * Class Implementation
   ****************************/
private:
  
  int                      num_curves;

  Point_list               all_points_list;
  Point_list               test_point_list;
  std::list<Pmwx::Locate_type> exp_type_list;
  
  unsigned                 expected_num_vertices, expected_num_edges, expected_num_faces;
    //expected_num_overlaps,
    //actual_num_overlaps;
 
  // count overlap references in arrangement,
  // that is every reference from a halfedge of an overlapped edge
  // to its overlapping curves.
    /*unsigned count_overlaps(Pmwx & pmwx)
      {
      Pmwx::Halfedge_iterator hit;
      Pmwx::Overlap_circulator oe;
      unsigned count, counted_overlaps = 0;
      
      std::cout << "halfedge: overlapping edges" << std::endl;
      for (hit = pmwx.halfedges_begin(); hit != pmwx.halfedges_end(); ++hit, ++hit) 
      {
      std::cout << (*hit).vertex()->point();
      std::cout << (*hit).opposite()->vertex()->point() << ": " << std::endl;
      oe=hit->overlap_edges();
      // we count how many edges refer to this halfedge
      // there is always at least one.
      // if there is more than one, there is an overlap
      count = 0;
      do {
      std::cout << "     ";
      std::cout << (*oe).halfedge()->vertex()->point();
      std::cout << (*oe).halfedge()->opposite()->vertex()->point() << std::endl;;
      count ++;
      } while (++oe != hit->overlap_edges());
      // we substract 1 from edges refering to this halfedge, see above
      counted_overlaps = counted_overlaps + (count - 1);
      std::cout << std::endl;
      }
      
      return counted_overlaps;
      } */
 
  void print_vertices(PM& pm)
    {
      typename PM::Vertex_const_iterator vit;

      std::cout << "Vertices in Pmwx:" << std::endl;
      for(vit = pm.vertices_begin(); vit != pm.vertices_end(); vit++)
	{
	  std::cout << (*vit).point() << " , ";
	}
      std::cout << std::endl;
   } 
 
  void print_kind_of_location(Pmwx::Locate_type &lt)
    {
      switch (lt) {
      case  PM::VERTEX:
	std::cout << "Vertex ";
	break;
      case  PM::EDGE:
	std::cout<< "Edge ";
	break;
      case  PM::FACE:
	std::cout<< "Face ";
	break;
      case  PM::UNBOUNDED_VERTEX:
	std::cout<< "UnBounded Vertex ";
	break;
      case  PM::UNBOUNDED_EDGE:
	std::cout<< "UnBounded Edge ";
	break;
      case  PM::UNBOUNDED_FACE:
	std::cout<< "UnBounded Face ";
	break;
      }
      std::cout << std::endl;
    }
  
  bool point_is_in_expected_place(PM& pm, Point &pnt, typename PM::Locate_type exp_lt)
    {
      typename PM::Locate_type    location_of_vertex;
      
      pm.locate(pnt ,location_of_vertex);
      print_kind_of_location(location_of_vertex);
      return (location_of_vertex == exp_lt);
    }
  
  void check_that_vertices_are_in_arrangement(PM & pm, Point_list & all_points_list)
    {
      Point_list::iterator pit;
      
      for (pit = all_points_list.begin(); pit != all_points_list.end(); pit++)
	{
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  std::cout << (*pit).xcoord() << " " << (*pit).ycoord() << "*** ";
#else
	  std::cout << (*pit).x() << " " << (*pit).y() << "*** ";
#endif
	  CGAL_assertion(point_is_in_expected_place(pm, *pit, PM::VERTEX) ||
                         point_is_in_expected_place(pm, *pit, PM::EDGE));
	}
    }
  
  void points_in_expected_place(PM &                    pm,
                                Point_list &            point_list,
                                std::list<typename PM::Locate_type> & lt_list)
    {
      Point_list::iterator              pit;
      typename std::list<typename PM::Locate_type>::iterator lt_it;
      
      for (pit = point_list.begin(), lt_it = lt_list.begin();
	   pit != point_list.end();
	   pit++, lt_it++)
	{
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  //#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  //CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || 
          //    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  std::cout << (*pit).xcoord() << " " << (*pit).ycoord() << "*** ";
#else
	  std::cout << (*pit).x() << " " << (*pit).y() << "*** ";
#endif
	  CGAL_assertion(point_is_in_expected_place(pm, *pit, *lt_it));
	}
    }

  void show_comparison()
    {
      std::cout << "expected # of vertices: ";
      std::cout << expected_num_vertices << std::endl;
      std::cout << "  actual # of vertices: ";
      std::cout << pm.number_of_vertices() << std::endl;

      std::cout << "expected # of edges: ";
      std::cout << expected_num_edges << std::endl;
      std::cout << "  actual # of edges: ";
      std::cout << pm.number_of_halfedges() / 2<< std::endl;      

      std::cout << "expected # of faces: ";
      std::cout << expected_num_faces << std::endl;
      std::cout << "  actual # of faces: ";
      std::cout << pm.number_of_faces() << std::endl;

      //std::cout << "expected # of overlaping edges: ";
      //std::cout << expected_num_overlaps << std::endl;
      //std::cout << "  actual # of overlaping edges: ";
      //std::cout << actual_num_overlaps << std::endl;
    }

  NT get_next_num(std::ifstream& file)
    {
      CGAL::Quotient<int> num;
      NT            result(INT_MAX);
      std::string   s;
      char          c = 0;

      //file.set_ascii_mode();
      while ( file && (result == NT(INT_MAX) ))
	{
	  // try to convert next token to integer
	  file >> c;

	  if (c=='#') // comment
	    {
	      std::getline(file, s);
	    }
	  else
	    {
	      file.putback(c);
	      file >> num;
              result = NT(num.numerator(), num.denominator());
	    }
	}

      // convertion failed, data file format error
      CGAL_assertion(result != NT(INT_MAX));

      return result;
    }

  int get_next_int(std::ifstream& file)
    {

#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
      // The to_long precondition is that number is indeed long
      // is supplied here since input numbers are small.
      return get_next_num(file).numerator().to_long();
#else
      return get_next_num(file).numerator();
#endif
    }

#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS

  Curve read_segment_curve(std::ifstream& file, bool reverse_order = false)
  {
    Curve segment;
    NT    x,y; 

    // read two segment points
    x = get_next_num(file); y = get_next_num(file);
    Point p1(x,y);
    x = get_next_num(file); y = get_next_num(file);
    Point p2(x,y);
    
    all_points_list.push_back(p1);
    all_points_list.push_back(p2);
    
    if (reverse_order)
      segment = Curve(p1,p2);
    else
      segment = Curve(p2,p1);

    return segment;
  }

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS

Curve read_polyline_curve(std::ifstream& file, bool reverse_order = false)
  {
    Curve                 polyline;
    NT                    x,y; 
    int                   num_x_curves ;
    Point_list            point_list;
    Point_list::iterator  plit;

    num_x_curves = get_next_int(file);
    
    while (num_x_curves--) {
      x = get_next_num(file); y = get_next_num(file);
      Point s(x,y);
      if (reverse_order)
	point_list.push_front(s);
      else
	point_list.push_back(s);
    }
    for (plit = point_list.begin(); //, cit = polyline.begin();
	 plit != point_list.end();
	 plit++) 
      {
	//std::cout << *plit << std::endl;
	polyline.push_back(*plit);
      }
    
    all_points_list.splice(all_points_list.end(), point_list); 
    return polyline;
}

#else
  #error No curve read function defined
#endif

  //#ifdef CGAL_PM_TEST_SWEEP
  void read_file_build_pm(std::ifstream& file, bool sweep_to_subcurves)
  { 
    NT          x,y;
    list<Curve> curves;
    Curve curr_curve;
    //SWEEP_LINE         sweep_line;
    // 1. read polylines and build arrangement
    
    // read number of polylines
    num_curves = get_next_int(file);
    
    // read curves (test specific)
    while (num_curves--) {
#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
      
      curr_curve = read_segment_curve(file);
      
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
      
      curr_curve = read_polyline_curve(file);
#endif
      
      //cout<<curr_curve<<endl;
      curves.push_back(curr_curve);
    }

    if (sweep_to_subcurves)
      CGAL::sweep_to_construct_planar_map(curves.begin(), curves.end(), pm);
    else{
      std::list<X_curve>   subcurves;
      Traits traits;
      CGAL::sweep_to_produce_planar_map_subcurves(curves.begin(), curves.end(), traits, subcurves);
      CGAL::sweep_to_construct_planar_map(subcurves.begin(), subcurves.end(), pm);
      
      //for (std::list<X_curve>::iterator scv_iter = subcurves.begin(); scv_iter != subcurves.end(); scv_iter++)
      //  pm.insert(*scv_iter);
    }

    // 2. read test vertices
    int num_test_points, exp_type;
    
    // read no. of test vertices
    num_test_points = get_next_int(file);
    
    while (num_test_points--) {
      x = get_next_num(file); y = get_next_num(file);
      std::cout << x << "," << y << std::endl;
      Point s(x,y);
      test_point_list.push_back(s);
      
      exp_type = get_next_int(file);
      exp_type_list.push_back( (typename PM::Locate_type) exp_type);
    }
    
    // 3. read expected arrangement properties
    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    expected_num_vertices = get_next_int(file);
    cout<<"The expected_num_vertices is "<<expected_num_vertices<<endl;
    
    //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
    expected_num_edges = get_next_int(file);
    cout<<"The expected_num_edge is "<<expected_num_edges<<endl;

    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    expected_num_faces = get_next_int(file);
    cout<<"The expected_num_faces is "<<expected_num_faces<<endl;
    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    //expected_num_overlaps = get_next_int(file);
    
  }

  /*#else
  void read_file_build_pm(std::ifstream& file, bool reverse_order)
  {
      NT    x,y; 
      Curve curr_curve;

      // 1. read polylines and build arrangement

      // read number of polylines
      num_curves = get_next_int(file);

      // read curves (test specific)
      while (num_curves--) {
#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS

        curr_curve = read_segment_curve(file, reverse_order);

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS

        curr_curve = read_polyline_curve(file, reverse_order);
#endif

	pm.insert(curr_curve);
      }

      // 2. read test vertices
      int num_test_points, exp_type;

      // read no. of test vertices
      num_test_points = get_next_int(file);

      while (num_test_points--) {
	x = get_next_num(file); y = get_next_num(file);
	std::cout << x << "," << y << std::endl;
	Point s(x,y);
	test_point_list.push_back(s);
	
	exp_type = get_next_int(file);
	exp_type_list.push_back( (typename PM::Locate_type) exp_type);
      }

      // 3. read expected arrangement properties
      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_vertices = get_next_int(file);

      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_edges = get_next_int(file);

      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_faces = get_next_int(file);

      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      //expected_num_overlaps = get_next_int(file);
      
    }
    #endif*/

  /****************************
   * Class Interface
   ****************************/
public:
  
  void start(char * filename, bool sweep_to_subcurves)
  {
    // Read data from file. Build Arrangement.
    std::ifstream file(filename);
    read_file_build_pm(file, sweep_to_subcurves);
    
    // Check validity of arrangement after insertion
    CGAL_assertion(pm.is_valid());
    
    // Check that vertices read are indeed in the arrangement
    check_that_vertices_are_in_arrangement(pm, all_points_list);
    
    // count overlaps
    //actual_num_overlaps = count_overlaps(pm); 
    
    show_comparison();   
    
    CGAL_assertion (pm.number_of_vertices()  == expected_num_vertices);
    // verify that test points are as located in the arrangemet as expected
    points_in_expected_place(pm, test_point_list, exp_type_list);
    CGAL_assertion (pm.number_of_halfedges() == expected_num_edges * 2);
    CGAL_assertion (pm.number_of_faces()     == expected_num_faces);
    //CGAL_assertion (actual_num_overlaps       == expected_num_overlaps);
  }
};

int main(int argc, char* argv[])
{
  
  Pm_polyline_traits_test<PM>   test;
  bool                          sweep_to_subcurves = false;

  if (argc < 2 || argc > 3) {
    std::cout << "usage: test data_file [subcurves]" << std::endl;
    exit(1);
  }

  //sweep_to_subcurves = (argc == 3 && 0 == strcmp(argv[2], "subcurves"));
  if (argc == 3) {
    std::string second_par(argv[2]);
    if (second_par.compare("subcurves") == 0) {
      sweep_to_subcurves = true;
    }
  }

  test.start(argv[1], sweep_to_subcurves);
  return 0;
}

#endif //CGAL_ARR_TEST_LEDA_CONFLICT
