#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Planar_map_2 Pm
#define Cartesian Cr
//#define Arr_polyline_traits APT
#define Quotient Qt
#define Planar_ma2p_traits_wrap PMTW
#define _List_iterator Lit
#define bidirectional_iterator_tag Bitt
//#define Planar_map_2 Pm2
#define Pm_default_dcel A2dd
#define Point_2 Pt2
#define allocator Altr
#define Td_X_trapezoid TdXt
#define Td_traits Tdt
#endif

#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#ifndef CGAL_PM_2_DEFAULT_DCEL_H
#include <CGAL/Pm_default_dcel.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_PM_IOSTREAM_H
#include <CGAL/IO/Pm_iostream.h>
#endif

#define CGAL_SEGMENT_TRAITS        1
#define CGAL_SEGMENT_LEDA_TRAITS   2


// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_PM_TEST_TRAITS
#define CGAL_PM_TEST_TRAITS CGAL_SEGMENT_TRAITS
//#define CGAL_PM_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS)

using namespace std;

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

#if CGAL_PM_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
#include <CGAL/Pm_segment_exact_traits.h>
#elif CGAL_PM_TEST_TRAITS==CGAL_SEGMENT_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Pm_leda_segment_exact_traits.h>
#else
  #error No traits defined for test
#endif

// Picking a default  point location strategy
// See comment above.
#ifndef CGAL_PM_TEST_POINT_LOCATION
#define CGAL_PM_TEST_POINT_LOCATION 1
  //#define CGAL_PM_TEST_POINT_LOCATION 2
  //#define CGAL_PM_TEST_POINT_LOCATION 3
#endif

#if CGAL_PM_TEST_POINT_LOCATION == 1
  // By default we use Trapezoidal Decomposition
#elif CGAL_PM_TEST_POINT_LOCATION == 2
  #include <CGAL/Pm_naive_point_location.h>  
#elif CGAL_PM_TEST_POINT_LOCATION == 3
  #include <CGAL/Pm_walk_along_line_point_location.h>
#else
  #error No point location strategy defined for test
#endif
 

#include <CGAL/Planar_map_2.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>

#if CGAL_PM_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
  typedef CGAL::Quotient<int>                  NT;
  typedef CGAL::Cartesian<NT>                  R;
  typedef CGAL::Pm_segment_exact_traits<R>    Traits;

#elif CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                        NT;
  typedef CGAL::Pm_leda_segment_exact_traits  Traits;

#endif

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef CGAL::Pm_default_dcel<Traits>        Dcel;

typedef CGAL::Planar_map_2<Dcel,Traits>      Planar_map;
 
// we use the namespace std for compatability with MSVC
typedef std::list<Point>                     Point_list;

class Pm_traits_test
{
  Planar_map   pm;  

public:
#if CGAL_PM_TEST_POINT_LOCATION == 3  
  Pm_traits_test() : pm(new CGAL::Pm_walk_along_line_point_location<Planar_map>) {};

#elif CGAL_PM_TEST_POINT_LOCATION == 2
  Pm_traits_test() : pm(new CGAL::Pm_naive_point_location<Planar_map>) {};
#else
  // None
#endif

  /****************************
   * Class Implementation
   ****************************/
private:
  
  int                      num_polylines;

  Point_list               all_points_list;
  Point_list               test_point_list;
  std::list<Planar_map::Locate_type> exp_type_list;
  
  unsigned                 expected_num_vertices,
                           expected_num_edges,
                           expected_num_faces;
 
  // count overlap references in arrangement,
  // that is every reference from a halfedge of an overlapped edge
  // to its overlapping curves.
  /*unsigned count_overlaps(Arr_2 & arr)
    {
      Planar_map::Halfedge_iterator hit;
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
      }*/
 
  void print_vertices(Planar_map &pm)
    {
      Planar_map::Vertex_const_iterator vit;

      std::cout << "Vertices in Planar_map:" << std::endl;
      for(vit = pm.vertices_begin(); vit != pm.vertices_end(); vit++)
	{
	  std::cout << (*vit).point() << " , ";
	}
      std::cout << std::endl;
   } 
 
  void print_kind_of_location(Planar_map::Locate_type &lt)
    {
      switch (lt) {
      case Planar_map::VERTEX:
	std::cout << "Vertex ";
	break;
      case Planar_map::EDGE:
	std::cout<< "Edge ";
	break;
      case Planar_map::FACE:
	std::cout<< "Face ";
	break;
      case Planar_map::UNBOUNDED_VERTEX:
	std::cout<< "UnBounded Vertex ";
	break;
      case Planar_map::UNBOUNDED_EDGE:
	std::cout<< "UnBounded Edge ";
	break;
      case Planar_map::UNBOUNDED_FACE:
	std::cout<< "UnBounded Face ";
	break;
      }
      std::cout << std::endl;
    }
  
  bool point_is_in_expected_place(Planar_map &pm, Point &pnt, 
				  Planar_map::Locate_type exp_lt)
    {
      Planar_map::Locate_type    location_of_vertex;
      
      pm.locate(pnt ,location_of_vertex);
      print_kind_of_location(location_of_vertex);
      return (location_of_vertex == exp_lt);
    }
  
  void check_that_vertices_are_in_arrangement(Planar_map &pm, 
					      Point_list & all_points_list)
    {
      Point_list::iterator pit;
      
      std::cout << "Following points should be on vertices or edges:";
      std::cout << std::endl;
      for (pit = all_points_list.begin(); pit != all_points_list.end(); pit++)
	{
#if CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
	  std::cout << (*pit).xcoord() << " " << (*pit).ycoord();
#else
	  std::cout << (*pit).x() << " " << (*pit).y();
#endif
	  std::cout << std::endl;
	  CGAL_assertion(point_is_in_expected_place(pm, *pit, 
						    Planar_map::VERTEX) ||
                         point_is_in_expected_place(pm, *pit, 
						    Planar_map::EDGE));
	}
    }
  
  void points_in_expected_place(  Planar_map &              pm,
				  Point_list &              point_list,
				  std::list<Planar_map::Locate_type> & lt_list)
    {
      Point_list::iterator              pit;
      std::list<Planar_map::Locate_type>::iterator lt_it;
      
      for (pit = point_list.begin(), lt_it = lt_list.begin();
	   pit != point_list.end();
	   pit++, lt_it++)
	{
#if CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
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

#if  CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
      // The to_long precondition is that number is indeed long
      // is supplied here since input numbers are small.
      return get_next_num(file).numerator().to_long();
#else
      return get_next_num(file).numerator();
#endif
    }

#if CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_TRAITS || \
    CGAL_PM_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS

  X_curve read_segment_curve(std::ifstream& file, bool reverse_order)
  {
    X_curve segment;
    NT    x,y; 

    // read two segment points
    x = get_next_num(file); y = get_next_num(file);
    Point p1(x,y);
    x = get_next_num(file); y = get_next_num(file);
    Point p2(x,y);
    
    all_points_list.push_back(p1);
    all_points_list.push_back(p2);
    
    if (reverse_order)
      segment = X_curve(p1,p2);
    else
      segment = X_curve(p2,p1);

    return segment;
  }

#elif CGAL_PM_TEST_TRAITS == CGAL_POLYLINE_TRAITS || \
      CGAL_PM_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS

X_curve read_polyline_curve(std::ifstream& file, bool reverse_order)
  {
    X_curve                 polyline;
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
      
      //std::cout<<s1<<endl;
      //std::cout<<s2<<endl;
      
      CGAL_assertion(std::string(s1) == std::string(s2));
    }
    
    CGAL_assertion (file1.eof() && file2.eof());

    std::cout<<"Reading and writing file  - O.K"<<std::endl;
  }
  
  void read_file_build_arrangement(std::ifstream& input_file, std::ifstream& file, bool reverse_order)
  {
      NT    x,y; 
      X_curve curr_curve;
      std::ofstream  pm_file("pm.txt", /*std::ios::in |*/ std::ios::out);
    
      pm_file.clear();

      //read planar map.
      //input_file >> pm;
      pm.read(input_file);
      pm_file << pm;
      
      // debugging!
      //std::cout<<pm;
      //input_file.close();
      //arr_file.close();
      //arr_file.open("temp", _IO_INPUT);
      //arr_file.flush();
      
      std::ifstream  pm_input_file("pm.txt", std::ios::in);
      compare_files(input_file, pm_input_file);
      
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
	exp_type_list.push_back( (Planar_map::Locate_type) exp_type);
      }

      // 3. read expected arrangement properties
      //      std::getline(file, s); // skip
      //      std::getline(file, s, ':'); // skip
      expected_num_vertices = get_next_int(file);
      std::cout<<"expected number of  vertices "<<expected_num_vertices;
      std::cout<<std::endl;
      
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
      //expected_num_overlaps = get_next_int(file);
      //std::cout<<"expected number of overlaps "<<expected_num_overlaps<<std::endl;
    }

  /****************************
   * Class Interface
   ****************************/
public:
  
  void start(char * input_filename, char * filename, bool reverse_order)
    {
      // Read data from file. Build Planar_map.
      std::ifstream input_file(input_filename);
      std::ifstream file(filename);
      
      read_file_build_arrangement(input_file, file, reverse_order);
      
      // DEBUG
      //print_vertices(arr);

      // debug
//       Arr_2::Face_handle    f  = arr->unbounded_face();
//       Arr_2::Holes_iterator it = f->holes_begin(),end=f->holes_end();
//       Arr_2::Ccb            c  = *it;
      //const X_curve& cv = curr->curve();

      // Check validity of arrangement after insertion
      CGAL_assertion(pm.is_valid());
            
      // Check that vertices read are indeed in the arrangement
      check_that_vertices_are_in_arrangement(pm, all_points_list);

      // count overlaps
      //actual_num_overlaps = count_overlaps(arr); 

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
  
  Pm_traits_test test;
  bool        reverse_order = false;

  if (argc < 3) {
    std::cout << "usage: test data_file [reverse]" << std::endl;
    exit(1);
  }
  
  //reverse_order = (argc == 3 && 0 == strcmp(argv[2], "reverse"));
  
  test.start(argv[1], argv[2], reverse_order);

  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT
