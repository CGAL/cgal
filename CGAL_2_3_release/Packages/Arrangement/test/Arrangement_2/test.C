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

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <fstream>

#define CGAL_SEGMENT_TRAITS        1
#define CGAL_SEGMENT_LEDA_TRAITS   2
#define CGAL_POLYLINE_TRAITS      11
#define CGAL_POLYLINE_LEDA_TRAITS 12
#define CGAL_SEGMENT_CIRCLE_TRAITS 21

// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
  //#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
  //#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
  //#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_TRAITS
  //#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_LEDA_TRAITS
  #define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_CIRCLE_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS)

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
//#error Currently not supported (July 2000)
  #include <CGAL/leda_rational.h>
  #include <CGAL/Arr_leda_polyline_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
  #include <CGAL/leda_real.h>
  #include <CGAL/Arr_segment_circle_traits.h>
#else
  #error No traits defined for test
#endif

// Picking a default  point location strategy
// See comment above.
#ifndef CGAL_ARR_TEST_POINT_LOCATION
  //#define CGAL_ARR_TEST_POINT_LOCATION 1 // Trapezoidal Decomposition
  #define CGAL_ARR_TEST_POINT_LOCATION 2 // Naive
  //#define CGAL_ARR_TEST_POINT_LOCATION 3 // Walk
  //#define CGAL_ARR_TEST_POINT_LOCATION 4 // Simple 
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
 
#include <CGAL/Arrangement_2.h>

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

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
  typedef leda_real                            NT;
  typedef CGAL::Arr_segment_circle_traits<NT>  Traits;
  typedef Traits::Segment                      Segment;
  typedef Traits::Circle                       Circle;

#endif

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef Traits::Curve                        Curve;
typedef CGAL::Arr_base_node<X_curve>         Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>     Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;
typedef Arr_2::Planar_map                           Planar_map;
 
// we use the namespace std for compatability with MSVC
typedef std::list<Point>                     Point_list;

class Arr_polyline_traits_test
{
  Arr_2 arr;  

public:
#if CGAL_ARR_TEST_POINT_LOCATION == 4
  Arr_polyline_traits_test() : arr(new CGAL::Pm_simple_point_location<Planar_map>) {};

#elif CGAL_ARR_TEST_POINT_LOCATION == 3  
  Arr_polyline_traits_test() : arr(new CGAL::Pm_walk_along_line_point_location<Planar_map>) {};

#elif CGAL_ARR_TEST_POINT_LOCATION == 2
  Arr_polyline_traits_test() : arr(new CGAL::Pm_naive_point_location<Planar_map>) {};
#else
  // CGAL_ARR_TEST_POINT_LOCATION == 1
  Arr_polyline_traits_test() : arr(new CGAL::Pm_default_point_location<Planar_map>) {};
  // None
#endif

  /****************************
   * Class Implementation
   ****************************/
private:
  
  int                      num_polylines;

  Point_list               all_points_list;
  Point_list               test_point_list;
  std::list<Arr_2::Locate_type> exp_type_list;
  
  unsigned                 expected_num_vertices,
                           expected_num_edges,
                           expected_num_faces,
                           expected_num_overlaps,
                           actual_num_overlaps;
 
  // count overlap references in arrangement,
  // that is every reference from a halfedge of an overlapped edge
  // to its overlapping curves.
  unsigned count_overlaps(Arr_2 & arr)
    {
      Arr_2::Halfedge_iterator hit;
      Arr_2::Overlap_circulator oe;
      unsigned count, counted_overlaps = 0;
      
      //std::cout << "halfedge: overlapping edges" << std::endl;
      for (hit=arr.halfedges_begin(); hit!=arr.halfedges_end(); ++hit, ++hit) 
	{
	  //std::cout << (*hit).vertex()->point();
	  //std::cout << (*hit).opposite()->vertex()->point() << ": " << std::endl;
	  oe=hit->overlap_edges();
	  // we count how many edges refer to this halfedge
	  // there is always at least one.
	  // if there is more than one, there is an overlap
	  count = 0;
	  do {
	    //std::cout << "     ";
	    //std::cout << (*oe).halfedge()->vertex()->point();
	    //std::cout << (*oe).halfedge()->opposite()->vertex()->point() << std::endl;;
	    count ++;
	  } while (++oe != hit->overlap_edges());
	  // we substract 1 from edges refering to this halfedge, see above
	  counted_overlaps = counted_overlaps + (count - 1);
	  //std::cout << std::endl;
	}
      
      return counted_overlaps; 
    }
 
  void print_vertices(Arr_2 & arr)
    {
      Arr_2::Vertex_const_iterator vit;

      std::cout << "Vertices in Arrangement:" << std::endl;
      for(vit = arr.vertices_begin(); vit != arr.vertices_end(); vit++)
	{
	  std::cout << (*vit).point() << " , ";
	}
      std::cout << std::endl;
   } 
 
  void print_kind_of_location(Arr_2::Locate_type &lt)
    {
      switch (lt) {
      case Arr_2::VERTEX:
	std::cout << "Vertex ";
	break;
      case Arr_2::EDGE:
	std::cout<< "Edge ";
	break;
      case Arr_2::FACE:
	std::cout<< "Face ";
	break;
      case Arr_2::UNBOUNDED_VERTEX:
	std::cout<< "UnBounded Vertex ";
	break;
      case Arr_2::UNBOUNDED_EDGE:
	std::cout<< "UnBounded Edge ";
	break;
      case Arr_2::UNBOUNDED_FACE:
	std::cout<< "UnBounded Face ";
	break;
      }
      std::cout << std::endl;
    }
  
  bool point_is_in_expected_place(Arr_2 & arr, Point &pnt, Arr_2::Locate_type exp_lt)
    {
      Arr_2::Locate_type    location_of_vertex;
      
      arr.locate(pnt ,location_of_vertex);
      print_kind_of_location(location_of_vertex);
      return (location_of_vertex == exp_lt);
    }
  
  void check_that_vertices_are_in_arrangement(Arr_2 & arr, Point_list & all_points_list)
    {
      Point_list::iterator pit;
      
      for (pit = all_points_list.begin(); pit != all_points_list.end(); pit++)
	{
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  std::cout << (*pit).xcoord() << " " << (*pit).ycoord() << "*** ";
#else
	  std::cout << (*pit).x() << " " << (*pit).y() << "*** ";
#endif
	  CGAL_assertion(point_is_in_expected_place(arr, *pit, Arr_2::VERTEX) ||
                         point_is_in_expected_place(arr, *pit, Arr_2::EDGE));
	}
    }
  
  void points_in_expected_place(Arr_2 &                      arr,
				  Point_list &              point_list,
				  std::list<Arr_2::Locate_type> & lt_list)
    {
      Point_list::iterator              pit;
      std::list<Arr_2::Locate_type>::iterator lt_it;
      
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
	  CGAL_assertion(point_is_in_expected_place(arr, *pit, *lt_it));
	}
    }

  void show_comparison()
    {
      std::cout << "expected # of vertices: ";
      std::cout << expected_num_vertices << std::endl;
      std::cout << "  actual # of vertices: ";
      std::cout << arr.number_of_vertices() << std::endl;

      std::cout << "expected # of edges: ";
      std::cout << expected_num_edges << std::endl;
      std::cout << "  actual # of edges: ";
      std::cout << arr.number_of_halfedges() / 2<< std::endl;      

      std::cout << "expected # of faces: ";
      std::cout << expected_num_faces << std::endl;
      std::cout << "  actual # of faces: ";
      std::cout << arr.number_of_faces() << std::endl;

      std::cout << "expected # of overlaping edges: ";
      std::cout << expected_num_overlaps << std::endl;
      std::cout << "  actual # of overlaping edges: ";
      std::cout << actual_num_overlaps << std::endl;
    }

  NT get_next_num(std::ifstream& file)
    {
      CGAL::Quotient<int> num = 0;
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
#if CGAL_ARR_TEST_TRAITS != CGAL_SEGMENT_CIRCLE_TRAITS
	      file >> num;
              result = NT(num.numerator(), num.denominator());
#else
	      num = num;
	      file >> result;
#endif
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
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
      return (int) CGAL::to_double(get_next_num(file));
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

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS

Curve read_seg_circ_curve(std::ifstream& file, bool reverse_order)
{
  Curve cv;
  
  // Get the arc type.
  char type;

  // Currently expects no comments in input file
  // Should be changed?
  file >> type;
  
  // A full circle (c) or a circular arc (a):
  if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
  {  
    // Read the circle, using the format "x0 y0 r^2"
    NT     x0, y0, r2;
    
    file >> x0 >> y0 >> r2;
//     x0 = get_next_num(file);
//     y0 = get_next_num(file);
//     r2 = get_next_num(file);
    
    Circle circle (Point (x0, y0), r2);

    if (type == 'c' || type == 'C')
    {
      // Create a full circle.
      cv = Curve(circle);  
    }
    else
    {
      // Read the end points of the circular arc.
      NT    x1, y1, x2, y2;

      file >> x1 >> y1 >> x2 >> y2;
//       x1 = get_next_num(file);
//       y1 = get_next_num(file);
//       x2 = get_next_num(file);
//       y2 = get_next_num(file);

      if ((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) != r2)
	y1 = CGAL::sqrt(r2 - (x1 - x0)*(x1 - x0)) + y0;

      if ((x2 - x0)*(x2 - x0) + (y2 - y0)*(y2 - y0) != r2)
	y2 = CGAL::sqrt(r2 - (x2 - x0)*(x2 - x0)) + y0;

      Point source (x1, y1);
      Point target (x2, y2);

      // Create the circular arc.
      cv = Curve (circle, source, target);
    }
  }
  else if (type == 's' || type == 'S')
  {
    // Read the end points of the segment.
    NT    x1, y1, x2, y2;

    file >> x1 >> y1 >> x2 >> y2;
//     x1 = get_next_num(file);
//     y1 = get_next_num(file);
//     x2 = get_next_num(file);
//     y2 = get_next_num(file);
    
    Point source (x1, y1);
    Point target (x2, y2);
   
    cv = Curve (Segment (source, target));
  }
  else
  {
    // Illegal type!
    std::cout << "Failed to read curve." << std::endl;
  }

  std::cout << "The read curve: " << cv << std::endl;
  return cv;
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

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS

        curr_curve = read_seg_circ_curve(file, reverse_order);

#else

#error No reading function defined for traits.

#endif

	arr.insert(curr_curve);
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
	exp_type_list.push_back( (Arr_2::Locate_type) exp_type);
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
      expected_num_overlaps = get_next_int(file);
      
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
      //print_vertices(arr);

      // debug
//       Arr_2::Face_handle    f  = arr->unbounded_face();
//       Arr_2::Holes_iterator it = f->holes_begin(),end=f->holes_end();
//       Arr_2::Ccb            c  = *it;
      //const X_curve& cv = curr->curve();

      // Check validity of arrangement after insertion
      CGAL_assertion(arr.is_valid());
            
      // Check that vertices read are indeed in the arrangement
      check_that_vertices_are_in_arrangement(arr, all_points_list);

      // count overlaps
      actual_num_overlaps = count_overlaps(arr); 

      show_comparison();   
      
      CGAL_assertion (arr.number_of_vertices()  == expected_num_vertices);
      // verify that test points are as located in the arrangemet as expected
      points_in_expected_place(arr, test_point_list, exp_type_list);
      CGAL_assertion (arr.number_of_halfedges() == expected_num_edges * 2);
      CGAL_assertion (arr.number_of_faces()     == expected_num_faces);
      CGAL_assertion (actual_num_overlaps       == expected_num_overlaps);

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
  
  int test_seed = rand();
  srand(test_seed); 
  std::cout << "Seed chosen for this run is " << test_seed << std::endl;

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
