#include "short_names.h"

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Pm_default_dcel.h>

#include <fstream>

#define CGAL_SEGMENT_TRAITS        1
#define CGAL_SEGMENT_LEDA_TRAITS   2
#define CGAL_POLYLINE_TRAITS      11
#define CGAL_POLYLINE_LEDA_TRAITS 12

// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_LEDA_TRAITS
#endif
   
// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && ( \
  CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
  CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS \
)

int main()
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
  #include <CGAL/Arr_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  #include <CGAL/leda_rational.h>
  #include <CGAL/Arr_leda_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  #include <CGAL/Arr_polyline_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
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
  typedef CGAL::Quotient<int>                                   NT;
  typedef CGAL::Cartesian<NT>                                   R;
  typedef CGAL::Arr_segment_traits_2<R>                         Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                                         NT;
  typedef CGAL::Arr_leda_segment_traits_2                       Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  typedef CGAL::Quotient<int>                                   NT;
  typedef CGAL::Cartesian<NT>                                   R;
  typedef CGAL::Arr_polyline_traits<R>                          Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
  typedef leda_rational                                         NT;
  typedef CGAL::Arr_leda_polyline_traits<NT>                    Traits;

#endif

typedef Traits::Point_2                                         Point_2;
typedef Traits::X_curve_2                                       X_curve_2;
typedef Traits::Curve_2                                         Curve_2;
typedef CGAL::Arr_base_node<X_curve_2>                          Base_node;
typedef CGAL::Pm_default_dcel<Traits>                           Dcel;

typedef CGAL::Planar_map_2<Dcel,Traits>                         Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>       Pmwx;
 

// we use the namespace std for compatability with MSVC
typedef std::list<Point_2>                                      Point_list;

class Arr_polyline_traits_test
{
  Pmwx subd;  

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
  
  int m_num_polylines;

  Point_list m_all_point_list;
  Point_list m_tst_point_list;
  std::list<Pmwx::Locate_type> exp_type_list;
  
  unsigned int m_exp_num_vertices;
  unsigned int m_exp_num_edges;
  unsigned int m_exp_num_faces;
  unsigned int m_exp_num_overlaps;
  unsigned int m_act_num_overlaps;
 
  void print_vertices(Pmwx & subd)
  {
    Pmwx::Vertex_const_iterator vit;

    std::cout << "Vertices in Arrangement:" << std::endl;
    for(vit = subd.vertices_begin(); vit != subd.vertices_end(); vit++) {
      std::cout << (*vit).point() << " , ";
    }
    std::cout << std::endl;
  } 
 
  void print_kind_of_location(Pmwx::Locate_type & lt)
  {
    switch (lt) {
      case Pmwx::VERTEX:
	std::cout << "Vertex ";
	break;
      case Pmwx::EDGE:
	std::cout<< "Edge ";
	break;
      case Pmwx::FACE:
	std::cout<< "Face ";
	break;
      case Pmwx::UNBOUNDED_VERTEX:
	std::cout<< "UnBounded Vertex ";
	break;
      case Pmwx::UNBOUNDED_EDGE:
	std::cout<< "UnBounded Edge ";
	break;
      case Pmwx::UNBOUNDED_FACE:
	std::cout<< "UnBounded Face ";
	break;
    }
      std::cout << std::endl;
    }
  
  bool point_is_in_expected_place(Pmwx & subd, Point_2 & pnt,
                                  Pmwx::Locate_type exp_lt)
  {
    Pmwx::Locate_type location_of_vertex;
      
    subd.locate(pnt ,location_of_vertex);
    print_kind_of_location(location_of_vertex);
    return (location_of_vertex == exp_lt);
  }
  
  void check_that_vertices_are_in_arrangement(Pmwx & subd,
                                              Point_list & m_all_point_list)
  {
    Point_list::iterator pit;
      
    for (pit = m_all_point_list.begin(); pit != m_all_point_list.end(); pit++)
    {
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
      std::cout << (*pit).xcoord() << " " << (*pit).ycoord() << "*** ";
#else
      std::cout << (*pit).x() << " " << (*pit).y() << "*** ";
#endif
      CGAL_assertion(point_is_in_expected_place(subd, *pit, Pmwx::VERTEX) ||
                     point_is_in_expected_place(subd, *pit, Pmwx::EDGE));
    }
  }
  
  void points_in_expected_place(Pmwx & subd,
				Point_list & point_list,
				std::list<Pmwx::Locate_type> & lt_list)
  {
    Point_list::iterator                          pit;
    std::list<Pmwx::Locate_type>::iterator lt_it;
      
    for (pit = point_list.begin(), lt_it = lt_list.begin();
         pit != point_list.end();
         pit++, lt_it++)
    {
#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
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
    std::cout << m_exp_num_vertices << std::endl;
    std::cout << "  actual # of vertices: ";
    std::cout << subd.number_of_vertices() << std::endl;

    std::cout << "expected # of edges: ";
    std::cout << m_exp_num_edges << std::endl;
    std::cout << "  actual # of edges: ";
    std::cout << subd.number_of_halfedges() / 2<< std::endl;      

    std::cout << "expected # of faces: ";
    std::cout << m_exp_num_faces << std::endl;
    std::cout << "  actual # of faces: ";
    std::cout << subd.number_of_faces() << std::endl;
  }

  NT get_next_num(std::ifstream & file)
  {
    CGAL::Quotient<int> num;
    NT result(INT_MAX);
    std::string s;
    char c = 0;

    //file.set_ascii_mode();
    while ( file && (result == NT(INT_MAX) )) {
      // try to convert next token to integer
      file >> c;

      // comment
      if (c == '#') std::getline(file, s);
      else {
        file.putback(c);
        file >> num;
        result = NT(num.numerator(), num.denominator());
      }
    }

    // convertion failed, data file format error
    CGAL_assertion(result != NT(INT_MAX));

    return result;
  }

  int get_next_int(std::ifstream & file)
  {

#if CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
    // The to_long precondition is that number is indeed long
    // is supplied here since input numbers are small.
    return get_next_num(file).numerator().to_long();
#else
    return get_next_num(file).numerator();
#endif
  }

#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS

  Curve_2 read_segment_curve(std::ifstream & file, bool reverse_order)
  {
    Curve_2 segment;
    NT x, y; 

    // read two segment points
    x = get_next_num(file); y = get_next_num(file);
    Point_2 p1(x,y);
    x = get_next_num(file); y = get_next_num(file);
    Point_2 p2(x,y);
    
    m_all_point_list.push_back(p1);
    m_all_point_list.push_back(p2);
    
    if (reverse_order)
      segment = Curve_2(p1,p2);
    else
      segment = Curve_2(p2,p1);

    return segment;
  }

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS

  Curve_2 read_polyline_curve(std::ifstream & file, bool reverse_order)
  {
    Curve_2 polyline;
    NT x, y; 
    int num_x_curves;
    Point_list point_list;
    Point_list::iterator plit;

    num_x_curves = get_next_int(file);
    
    while (num_x_curves--) {
      x = get_next_num(file); y = get_next_num(file);
      Point_2 s(x, y);
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
    
    m_all_point_list.splice(m_all_point_list.end(), point_list); 
    return polyline;
  }

#else
#error No curve read function defined
#endif

  void read_file_build_arrangement(std::ifstream& file, bool reverse_order)
  {
    NT x, y; 
    Curve_2 curr_curve;

    // 1. read polylines and build arrangement

    // read number of polylines
    m_num_polylines = get_next_int(file);

    // read curves (test specific)
    while (m_num_polylines--) {
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
      Point_2 s(x,y);
      m_tst_point_list.push_back(s);
	
      exp_type = get_next_int(file);
      exp_type_list.push_back( (Pmwx::Locate_type) exp_type);
    }

    // 3. read expected arrangement properties
    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    m_exp_num_vertices = get_next_int(file);

    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    m_exp_num_edges = get_next_int(file);

    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip
    m_exp_num_faces = get_next_int(file);

    //      std::getline(file, s); // skip
    //      std::getline(file, s, ':'); // skip

    // Shai: There are no overlaps in pmwx.      
    //m_exp_num_overlaps = get_next_int(file);  
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

    // DEBUG: Print polyline content
#if 0
    typedef Pmwx::Halfedge_const_iterator Halfedge_const_iterator;
    Halfedge_const_iterator it = subd.halfedges_begin();
    Halfedge_const_iterator it_end = subd.halfedges_end();
    int num = 1;
    CGAL_For_all(it, it_end) {
      std::cout << "Curve " << num << std::endl;
      std::copy(it->curve().begin(), it->curve().end(), 
                std::ostream_iterator<Point_2>(std::cout, " "));
      std::cout << std::endl << std::endl;
      ++it; ++num;
    }
      
    print_vertices(subd);

    // Pmwx::Face_handle f  = subd->unbounded_face();
    // Pmwx::Holes_iterator it = f->holes_begin(), end=f->holes_end();
    // Pmwx::Ccb c  = *it;
    // const X_curve_2_2 & cv = curr->curve();
#endif
    
    // Check validity of arrangement after insertion:
    CGAL_assertion(subd.is_valid());
            
    // Check that vertices read are indeed in the arrangement:
    check_that_vertices_are_in_arrangement(subd, m_all_point_list);

    // count overlaps
    // Shai: There are no overlaps in pmwx.
    // m_act_num_overlaps = count_overlaps(subd); 

    show_comparison();   
      
    CGAL_assertion (subd.number_of_vertices()  == m_exp_num_vertices);
    // verify that test points are as located in the arrangemet as expected
    points_in_expected_place(subd, m_tst_point_list, exp_type_list);
    CGAL_assertion (subd.number_of_halfedges() == m_exp_num_edges * 2);
    CGAL_assertion (subd.number_of_faces()     == m_exp_num_faces);

    // Shai: There are no overlaps in pmwx.
    // CGAL_assertion (m_act_num_overlaps       == m_exp_num_overlaps);
  }
};

int main(int argc, char* argv[])
{
  Arr_polyline_traits_test test;
  bool reverse_order = false;

  if (argc < 2 || argc > 3) {
    std::cout << "usage: test data_file [reverse]" << std::endl;
    exit(1);
  }

  // reverse_order = (argc == 3 && 0 == strcmp(argv[2], "reverse"));
  if (argc == 3) {
    std::string second_par(argv[2]);
    if (second_par.compare("reverse") == 0) {
      reverse_order = true;
    }
  }

  test.start(argv[1], reverse_order);
  return 0;
}

#endif
