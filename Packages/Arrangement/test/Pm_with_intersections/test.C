#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define Cartesian Cr
#define Arr_polyline_traits APT
#define Quotient Qt
#define Planar_map_traits_wrap PMTW
#define _List_iterator Lit
#define bidirectional_iterator_tag Bitt
#define Planar_map_2 Pm2
#define Arr_2_default_dcel A2dd
#define Point_2 Pt2
#define allocator Altr
#define Td_X_trapezoid TdXt
#define Td_traits Tdt
#endif

#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
//#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Pm_default_dcel.h>


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
  #define CGAL_ARR_TEST_POINT_LOCATION 1
  //#define CGAL_ARR_TEST_POINT_LOCATION 2
  //#define CGAL_ARR_TEST_POINT_LOCATION 3
#endif

#if CGAL_ARR_TEST_POINT_LOCATION == 1
  // By default we use Trapezoidal Decomposition
  #include <CGAL/Pm_default_point_location.h>  
#elif CGAL_ARR_TEST_POINT_LOCATION == 2
  #include <CGAL/Pm_naive_point_location.h>  
#elif CGAL_ARR_TEST_POINT_LOCATION == 3
  #include <CGAL/Pm_walk_along_line_point_location.h>
#elif CGAL_ARR_TEST_POINT_LOCATION == 4
  #include <CGAL/Pm_simple_point_location.h>
#else
  #error No point location strategy defined for test
#endif
 
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
  typedef CGAL::Quotient<int>                  NT;
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

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef Traits::Curve                        Curve;
typedef CGAL::Arr_base_node<X_curve>         Base_node;
//typedef CGAL::Arr_2_default_dcel<Traits>     Dcel;
typedef CGAL::Pm_default_dcel<Traits>        Dcel;

typedef CGAL::Planar_map_2<Dcel,Traits>                    Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>  Subdivision;
//typedef Subdivision::Planar_map                           Planar_map;
 

// we use the namespace std for compatability with MSVC
typedef std::list<Point>                     Point_list;

class Arr_polyline_traits_test
{
  Subdivision subd;  

public:
#if CGAL_ARR_TEST_POINT_LOCATION == 3  
  Arr_polyline_traits_test() : 
    subd(new CGAL::Pm_walk_along_line_point_location<Planar_map>) {};

#elif CGAL_ARR_TEST_POINT_LOCATION == 4  
  Arr_polyline_traits_test() : 
    subd(new CGAL::Pm_simple_point_location<Planar_map>) {};

#elif CGAL_ARR_TEST_POINT_LOCATION == 2
  Arr_polyline_traits_test() : 
    subd(new CGAL::Pm_naive_point_location<Planar_map>) {};
#else
  // CGAL_ARR_TEST_POINT_LOCATION == 1
  Arr_polyline_traits_test() : 
    subd(new CGAL::Pm_default_point_location<Planar_map>) {};
  // None
#endif

  /****************************
   * Class Implementation
   ****************************/
private:
  
  int                      num_polylines;

  Point_list               all_points_list;
  Point_list               test_point_list;
  std::list<Subdivision::Locate_type> exp_type_list;
  
  unsigned                 expected_num_vertices,
                           expected_num_edges,
                           expected_num_faces,
                           expected_num_overlaps,
                           actual_num_overlaps;
 
  void print_vertices(Subdivision & subd)
    {
      Subdivision::Vertex_const_iterator vit;

      std::cout << "Vertices in Arrangement:" << std::endl;
      for(vit = subd.vertices_begin(); vit != subd.vertices_end(); vit++)
	{
	  std::cout << (*vit).point() << " , ";
	}
      std::cout << std::endl;
   } 
 
  void print_kind_of_location(Subdivision::Locate_type &lt)
    {
      switch (lt) {
      case Subdivision::VERTEX:
	std::cout << "Vertex ";
	break;
      case Subdivision::EDGE:
	std::cout<< "Edge ";
	break;
      case Subdivision::FACE:
	std::cout<< "Face ";
	break;
      case Subdivision::UNBOUNDED_VERTEX:
	std::cout<< "UnBounded Vertex ";
	break;
      case Subdivision::UNBOUNDED_EDGE:
	std::cout<< "UnBounded Edge ";
	break;
      case Subdivision::UNBOUNDED_FACE:
	std::cout<< "UnBounded Face ";
	break;
      }
      std::cout << std::endl;
    }
  
  bool point_is_in_expected_place(Subdivision & subd, Point &pnt, Subdivision::Locate_type exp_lt)
    {
      Subdivision::Locate_type location_of_vertex;
      
      subd.locate(pnt ,location_of_vertex);
      print_kind_of_location(location_of_vertex);
      return (location_of_vertex == exp_lt);
    }
  
  void check_that_vertices_are_in_arrangement(Subdivision & subd, Point_list & all_points_list)
    {
      Point_list::iterator pit;
      
      for (pit = all_points_list.begin(); pit != all_points_list.end(); pit++)
	{
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  std::cout << (*pit).xcoord() << " " << (*pit).ycoord() << "*** ";
#else
	  std::cout << (*pit).x() << " " << (*pit).y() << "*** ";
#endif
	  CGAL_assertion(point_is_in_expected_place(subd, *pit, Subdivision::VERTEX) ||
                         point_is_in_expected_place(subd, *pit, Subdivision::EDGE));
	}
    }
  
  void points_in_expected_place(Subdivision & subd,
				Point_list &  point_list,
				std::list<Subdivision::Locate_type> & lt_list)
    {
      Point_list::iterator                          pit;
      std::list<Subdivision::Locate_type>::iterator lt_it;
      
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
	  CGAL_assertion(point_is_in_expected_place(subd, *pit, *lt_it));
	}
      std::cout << std::endl;
    }

  void show_comparison()
    {
      std::cout << "expected # of vertices: ";
      std::cout << expected_num_vertices << std::endl;
      std::cout << "  actual # of vertices: ";
      std::cout << subd.number_of_vertices() << std::endl;

      std::cout << "expected # of edges: ";
      std::cout << expected_num_edges << std::endl;
      std::cout << "  actual # of edges: ";
      std::cout << subd.number_of_halfedges() / 2<< std::endl;      

      std::cout << "expected # of faces: ";
      std::cout << expected_num_faces << std::endl;
      std::cout << "  actual # of faces: ";
      std::cout << subd.number_of_faces() << std::endl;

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

  Curve read_segment_curve(std::ifstream& file, bool reverse_order)
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

Curve read_polyline_curve(std::ifstream& file, bool reverse_order)
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

  void read_file_build_arrangement(std::ifstream& file, bool reverse_order)
  {
      NT    x,y; 
      Curve curr_curve;

      // 1. read polylines and build arrangement

      // read number of polylines
      num_polylines = get_next_int(file);

      // read curves (test specific)
      while (num_polylines--) {
#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS

        curr_curve = read_segment_curve(file, reverse_order);

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS

        curr_curve = read_polyline_curve(file, reverse_order);
#endif

	subd.insert(curr_curve);
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
	exp_type_list.push_back( (Subdivision::Locate_type) exp_type);
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

      // Shai: There are no overlaps in pmwx.      
      //expected_num_overlaps = get_next_int(file);
      
    }

  /****************************
   * Class Interface
   ****************************/
public:
  
  void start(char * filename, bool reverse_order)
    {
      // Read data from file. Build Arrangement.
      std::ifstream file(filename);
      read_file_build_arrangement(file, reverse_order);
      
      // DEBUG
      //print_vertices(subd);

      // debug
//       Subdivision::Face_handle    f  = subd->unbounded_face();
//       Subdivision::Holes_iterator it = f->holes_begin(),end=f->holes_end();
//       Subdivision::Ccb            c  = *it;
      //const X_curve& cv = curr->curve();

      // Check validity of arrangement after insertion
      CGAL_assertion(subd.is_valid());
            
      // Check that vertices read are indeed in the arrangement
      check_that_vertices_are_in_arrangement(subd, all_points_list);

      // count overlaps
      // Shai: There are no overlaps in pmwx.
      //actual_num_overlaps = count_overlaps(subd); 

      show_comparison();   
      
      CGAL_assertion (subd.number_of_vertices()  == expected_num_vertices);
      // verify that test points are as located in the arrangemet as expected
      points_in_expected_place(subd, test_point_list, exp_type_list);
      CGAL_assertion (subd.number_of_halfedges() == expected_num_edges * 2);
      CGAL_assertion (subd.number_of_faces()     == expected_num_faces);
      // Shai: There are no overlaps in pmwx.
      //CGAL_assertion (actual_num_overlaps       == expected_num_overlaps);

    }
};

int main(int argc, char* argv[])
{
  
  Arr_polyline_traits_test test;
  bool                     reverse_order = false;

  if (argc < 2 || argc > 3) {
    std::cout << "usage: test data_file [reverse]" << std::endl;
    exit(1);
  }

  //reverse_order = (argc == 3 && 0 == strcmp(argv[2], "reverse"));
  if (argc == 3) {
    std::string second_par(argv[2]);
    if (second_par.compare("reverse") == 0) {
      reverse_order = true;
    }
  }

  test.start(argv[1], reverse_order);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT
