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

#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>

#ifndef CGAL_ARR_IOSTREAM_H
#include <CGAL/IO/Arr_iostream.h>
#endif

/*#ifndef CGAL_IO_ARR_WINDOW_STREAM_H
  #include <CGAL/IO/Arr_Window_stream.h>
  #endif
  
  #ifndef CGAL_IO_ARR_POSTSCRIPT_FILE_STREAM_H
  #include <CGAL/IO/Arr_Postscript_file_stream.h>
  #endif
  
  #ifndef CGAL_LEDA_WINDOW_H
  #include <CGAL/IO/leda_window.h>
  #endif*/

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
#elif CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>
#elif CGAL_ARR_TEST_TRAITS==CGAL_POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits.h>
//#include <CGAL/IO/Arr_polyline_traits_iostream.h>
#elif CGAL_ARR_TEST_TRAITS==CGAL_POLYLINE_LEDA_TRAITS
//#error Currently not supported (July 2000)
#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_polyline_traits.h>
//#include <CGAL/IO/Arr_leda_polyline_traits_iostream.h>
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
  // Trapezoidal Decomposition
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
 
// Using my own temporary version which include arr.is_valid()
//#include "Arrangement_2_Debug.h"
//#include "Arrangement_2_Shai.h"
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

#endif

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef Traits::Curve                        Curve;
typedef CGAL::Arr_base_node<X_curve>         Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>     Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;
typedef Arr_2::Planar_map                           Planar_map;
 
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
  // Trapezoidal decomposition CGAL_ARR_TEST_POINT_LOCATION == 1
  Arr_polyline_traits_test() : arr(new CGAL::Pm_default_point_location<Planar_map>) {};
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
      
      std::cout << "halfedge: overlapping edges" << std::endl;
      for (hit=arr.halfedges_begin(); hit!=arr.halfedges_end(); ++hit, ++hit) 
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

  void compare_files(std::ifstream& file1, std::ifstream& file2)
  {
    const int STR_LEN = 80;
    
    char s1[STR_LEN], s2[STR_LEN];

    file1.seekg(0, std::ios::beg);
    file2.seekg(0, std::ios::beg);
    while (!file1.eof() && !file2.eof()){
      file1.getline(s1, STR_LEN);
      file2.getline(s2, STR_LEN);
      
      //cout<<s1<<endl;
      //cout<<s2<<endl;
      
       CGAL_assertion(std::string(s1) == std::string(s2));
    }
    
     CGAL_assertion(file1.eof() && file2.eof());

    std::cout<<"Reading and writing file  - O.K"<<std::endl;
  }
  
  void read_file_build_arrangement(std::ifstream& input_file, 
				   std::ifstream& file, bool reverse_order)
  {
      NT    x,y; 
      Curve curr_curve;
      std::ofstream  arr_file("arr.txt" , /*std::ios::in |*/ std::ios::out);
      
      arr_file.clear();

      //read arrangement
      input_file >> arr;
      arr_file << arr;
      
      //std::cout<<arr;
      //input_file.close();
      //arr_file.close();
      //arr_file.open("temp", _IO_INPUT);
      //arr_file.flush();
      
      std::ifstream  arr_input_file("arr.txt", std::ios::in);
      compare_files(input_file, arr_input_file);
      
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
      std::cout<<"expected number of  vertices ";
      std::cout<<expected_num_vertices<<std::endl;
      
      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_edges = get_next_int(file);
      std::cout<<"expected number of edges "<< expected_num_edges<<std::endl;
      
      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_faces = get_next_int(file);
      std::cout<<"expected number of faces "<<expected_num_faces <<std::endl;
      
      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_overlaps = get_next_int(file);
      std::cout<<"expected number of overlaps "<<expected_num_overlaps;
      std::cout<<std::endl;
    }

  /****************************
   * Class Interface
   ****************************/
public:
  
  void start(char * input_filename, char * filename, bool reverse_order)
    {
      // Read data from file. Build Arrangement.
      std::ifstream file(filename);
      std::ifstream input_file(input_filename);
      
      read_file_build_arrangement(input_file, file, reverse_order);
      
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
  bool        reverse_order = false;

  int test_seed = rand();
  srand(test_seed); 
  std::cout << "Seed chosen for this run is " << test_seed << std::endl;

  if (argc < 3) {
    std::cout << "usage: test data_file_io data_file" << std::endl;
    exit(1);
  }
  
  //reverse_order = (argc == 4 && 0 == strcmp(argv[3], "reverse"));

  test.start(argv[1], argv[2], reverse_order);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT
