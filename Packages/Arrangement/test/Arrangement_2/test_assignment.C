#include "short_names.h"

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
#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_LEDA_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_CIRCLE_TRAITS
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
  #include <CGAL/Arr_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  #include <CGAL/leda_rational.h>
  #include <CGAL/Pm_segment_traits_leda_kernel_2.h>
  #include <CGAL/Arr_leda_segment_traits_2.h>
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
  typedef CGAL::Quotient<int>                           NT;
  typedef CGAL::Cartesian<NT>                           R;
  typedef CGAL::Arr_segment_traits_2<R>                 Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                                 NT;
  typedef CGAL::Pm_segment_traits_leda_kernel_2         Kernel;
  typedef CGAL::Arr_leda_segment_traits_2<Kernel>       Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  typedef CGAL::Quotient<int>                           NT;
  typedef CGAL::Cartesian<NT>                           R;
  typedef CGAL::Arr_polyline_traits<R>                  Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
  typedef leda_rational                                 NT;
  typedef CGAL::Arr_leda_polyline_traits<>              Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
  typedef leda_real                                     NT;
  typedef CGAL::Arr_segment_circle_traits<NT>           Traits;
  typedef Traits::Segment                               Segment;
  typedef Traits::Circle                                Circle;

#endif

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef Traits::Curve                        Curve;
typedef CGAL::Arr_base_node<X_curve>         Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>     Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;
typedef Arr_2::Planar_map                           Planar_map;
 
typedef Arr_2::Curve_iterator                  Curve_iterator;
typedef Arr_2::Subcurve_iterator               Subcurve_iterator;
typedef Arr_2::Edge_iterator                   Edge_iterator;
typedef Arr_2::Curve_const_iterator            Curve_const_iterator;
typedef Arr_2::Subcurve_const_iterator         Subcurve_const_iterator;
typedef Arr_2::Edge_const_iterator             Edge_const_iterator;
typedef Arr_2::Overlap_circulator              Overlap_circulator;
typedef Arr_2::Overlap_const_circulator        Overlap_const_circulator;

typedef Arr_2::Vertex_iterator                 Vertex_iterator;
typedef Arr_2::Halfedge_iterator               Halfedge_iterator;
typedef Arr_2::Face_iterator                   Face_iterator;
typedef Arr_2::Vertex_const_iterator           Vertex_const_iterator;
typedef Arr_2::Halfedge_const_iterator         Halfedge_const_iterator;
typedef Arr_2::Face_const_iterator             Face_const_iterator;
typedef Arr_2::Ccb_halfedge_circulator         Ccb_halfedge_circulator;
typedef Arr_2::Ccb_halfedge_const_circulator   Ccb_halfedge_const_circulator;
typedef Arr_2::Holes_iterator                  Holes_iterator;
typedef Arr_2::Holes_const_iterator            Holes_const_iterator;


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
  //std::list<Arr_2::Locate_type> exp_type_list;
 
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
    //NT    x,y; 
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

      /*
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
      expected_num_overlaps = get_next_int(file);*/
      
    }
  
  template <class Circulator>
  typename Circulator::size_type Circ_size(Circulator circ) const 
  {
    typedef typename Circulator::size_type  size_type;
    size_type n = 0;
    Circulator   d = circ;
    do {
        ++n;
        ++d;
    } while( circ != d);
    return n;
  }
    
  void compare_planar_maps(const Arr_2& arr1, const Arr_2& arr2)
  {
    Vertex_const_iterator  v_iter1, v_iter2;
    
    CGAL_assertion(arr1.number_of_vertices() == arr2.number_of_vertices());
    std::cout<<"copied arrangement and original arrangement have the same number of vertices"<<std::endl;
  
    // Comparing Vertices.
    for (v_iter1 = arr1.vertices_begin(), v_iter2 = arr2.vertices_begin(); 
         v_iter1 != arr1.vertices_end() && v_iter2 != arr2.vertices_end(); ++v_iter1, ++v_iter2)
      CGAL_assertion(v_iter1->point() == v_iter2->point());
    
    std::cout<<"copied vertices == original vertices"<< std::endl;

    
    CGAL_assertion(arr1.number_of_halfedges() == arr2.number_of_halfedges());
    std::cout<<"copied arrangement and original arrangement have the same number of halfedges"<<std::endl;
    
    // Comparing Halfedges.
    Halfedge_const_iterator  h_iter1, h_iter2;
    for (h_iter1 = arr1.halfedges_begin(), h_iter2 = arr2.halfedges_begin(); 
         h_iter1 != arr1.halfedges_end() && h_iter2 != arr2.halfedges_end(); ++h_iter1, ++h_iter2){
      CGAL_assertion(h_iter1->curve() == h_iter2->curve());
      CGAL_assertion(h_iter1->source()->point() == h_iter2->source()->point());
      CGAL_assertion(h_iter1->target()->point() == h_iter2->target()->point());
    }
    
    std::cout<<"copied halfedges == original halfedges"<< std::endl; 


    CGAL_assertion(arr1.number_of_faces() == arr2.number_of_faces());
    std::cout<<"copied arrangement and original arrangement have the same number of faces"<<std::endl;
    
    // Comparing Faces.
    Face_const_iterator  f_iter1, f_iter2;
    for (f_iter1 = arr1.faces_begin(), f_iter2 = arr2.faces_begin(); 
         f_iter1 != arr1.faces_end() && f_iter2 != arr2.faces_end(); ++f_iter1, ++f_iter2){
      
      CGAL_assertion(f_iter1->is_unbounded() == f_iter2->is_unbounded());

      if (!f_iter1->is_unbounded()){
        Ccb_halfedge_const_circulator cc1 = f_iter1->outer_ccb();
        Ccb_halfedge_const_circulator cc2 = f_iter2->outer_ccb();
        
        CGAL_assertion(Circ_size(cc1) == Circ_size(cc2));
        std::cout << "copied outer ccb and original outer ccb have the same size"<< std::endl;
          
        do {
          CGAL_assertion(cc1->curve() == cc2->curve());
        } while (++cc1 != f_iter1->outer_ccb() && ++cc2 != f_iter2->outer_ccb());

        // asserting both faces have the same size.
        //CGAL_assertion(cc1 == f_iter1->outer_ccb() && cc2 == f_iter2->outer_ccb());

        std::cout<<"copied outer ccb == original outer ccb"<<std::endl;
        
        //CGAL_assertion(std::distance(f_iter1->holes_begin(), f_iter1->holes_end()) && 
        //               std::distance(f_iter2->holes_begin(), f_iter2->holes_end()) );
        //std::cout<<"copied arrangement and original arrangement have the same number of holes"<<std::endl;

        Holes_const_iterator hole_iter1, hole_iter2; 
        for (hole_iter1 = f_iter1->holes_begin(), hole_iter2 = f_iter2->holes_begin();
             hole_iter1 != f_iter1->holes_end() && hole_iter2 != f_iter2->holes_end(); 
             ++hole_iter1, ++hole_iter2) {
          Ccb_halfedge_const_circulator cch1(*hole_iter1), cch2(*hole_iter2);
          
          CGAL_assertion(Circ_size(cch1) == Circ_size(cch2));
          std::cout<<"copied hole and original hole have the same size"<< std::endl;
          
          do{
            CGAL_assertion(cch1->curve() == cch2->curve());
          } while (++cch1 != *hole_iter1 && ++cch2 != *hole_iter2);
          std::cout<<"copied hole == original hole"<<std::endl;
        }
        
        CGAL_assertion(hole_iter1 ==  f_iter1->holes_end() && hole_iter2 == f_iter2->holes_end());
        std::cout<<"copied arrangement and original arrangement have the same number of holes"<<std::endl;
      }
    }
  }
  
  void  compare_hierarchy_tree(const Arr_2& arr1, const Arr_2& arr2)
  {
    Curve_const_iterator cn_iter1, cn_iter2;

    CGAL_assertion(arr1.number_of_curve_nodes() == arr2.number_of_curve_nodes());
    std::cout<<"copied arrangement and original arrangement have the same number of curve nodes"<< std::endl;
    
    for (cn_iter1 = arr1.curve_node_begin(), cn_iter2 = arr2.curve_node_begin();
         cn_iter1 != arr1.curve_node_end() && cn_iter2 != arr2.curve_node_end();
         ++cn_iter1, ++cn_iter2){
      
      // comparing curve for curve node.
      CGAL_assertion(cn_iter1->curve() == cn_iter2->curve());
      std::cout<<"copied curve node and original curve node have the same curve"<< std::endl;
      
      CGAL_assertion(cn_iter1->number_of_sc_levels() == cn_iter2->number_of_sc_levels());
      std::cout<<"copied curve node and original curve node have the same number of levels"<< std::endl;
      
      // comparing all subcurves nodes level.
      for (unsigned int i = 0; i < cn_iter1->number_of_sc_levels(); ++i){
        Subcurve_const_iterator scv_iter1, scv_iter2;
        // comparing all subcurves nodes data.
        for (scv_iter1 = cn_iter1->level_begin(i), scv_iter2 = cn_iter2->level_begin(i);
             scv_iter1 != cn_iter1->level_end(i), scv_iter2 != cn_iter2->level_end(i);
             ++scv_iter1, ++scv_iter2){
          
          CGAL_assertion(scv_iter1->curve() == scv_iter2->curve());
          CGAL_assertion(scv_iter1->is_edge_node() == scv_iter2->is_edge_node());
          CGAL_assertion(scv_iter1->curve_node()->curve() == scv_iter2->curve_node()->curve());
          CGAL_assertion(scv_iter1->parent()->curve() == scv_iter2->parent()->curve());
          CGAL_assertion(scv_iter1->children_begin()->curve() == scv_iter2->children_begin()->curve());

          Subcurve_const_iterator  last_child1 = scv_iter1->children_end(), 
            last_child2 = scv_iter2->children_end();

          CGAL_assertion( (--last_child1)->curve() == (--last_child2)->curve());

          CGAL_assertion(scv_iter1->edges_begin()->curve() == scv_iter2->edges_begin()->curve());
          
          Edge_const_iterator  last_edge1 = scv_iter1->edges_end(), 
            last_edge2 = scv_iter2->edges_end();

          CGAL_assertion( (--last_edge1)->curve() == (--last_edge2)->curve());
        } 
      }
      
      std::cout<<"copied hierarchy tree and original hierarchy tree are the same"<< std::endl;
      
      // comparing edge nodes.
      Edge_const_iterator edge_iter1, edge_iter2;
      for (edge_iter1 = cn_iter1->edges_begin(), edge_iter2 = cn_iter2->edges_begin();
           edge_iter1 != cn_iter1->edges_end(), edge_iter2 != cn_iter2->edges_end();
           ++edge_iter1, ++edge_iter2){
        
        CGAL_assertion(edge_iter1->is_edge_node() == edge_iter2->is_edge_node());
        CGAL_assertion(edge_iter1->halfedge()->curve() == edge_iter2->halfedge()->curve());
      }
      
      std::cout<<"copied edge nodes list and original edge nodes list are the same"<< std::endl;
    }

    // Finally checking for equality with the overlapping edges.
    Halfedge_const_iterator h_iter1, h_iter2;
    for (h_iter1 = arr1.halfedges_begin(), h_iter2 = arr2.halfedges_begin(); 
         h_iter1 != arr1.halfedges_end() && h_iter2 != arr2.halfedges_end(); ++h_iter1, ++h_iter2){
      CGAL_assertion(h_iter1->edge_node()->curve() == h_iter2->edge_node()->curve());
      
      Overlap_const_circulator  ovlp_edges1 = h_iter1->overlap_edges(), ovlp_edges2 = h_iter2->overlap_edges();
      CGAL_assertion(Circ_size(ovlp_edges1) == Circ_size(ovlp_edges2));
      std::cout<<"copied hole and original hole have the same size"<< std::endl;

      do{
        CGAL_assertion(ovlp_edges1->curve() == ovlp_edges2->curve());
      } while (++ovlp_edges1 != h_iter1->overlap_edges() && ++ovlp_edges2 != h_iter2->overlap_edges());
    }
  }

  void  compare_arr(const Arr_2& arr1, const Arr_2& arr2)
  {
    compare_planar_maps(arr1, arr2);
    compare_hierarchy_tree(arr1, arr2);
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
      
      // Check validity of arrangement after insertion
      CGAL_assertion(arr.is_valid());

      Arr_2  copy_arr(arr);
      
      // Check validity of arrangement after copying
      CGAL_assertion(copy_arr.is_valid());
      
      compare_arr(copy_arr,arr);
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
