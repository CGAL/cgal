/***********************************************************************



************************************************************************/
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>


typedef double coord_type;
typedef CGAL::Simple_cartesian<coord_type>  SC;
typedef CGAL::Filtered_kernel<SC> K;

typedef K::Point_3  Point;
typedef K::Segment_3  Segment;


typedef CGAL::Alpha_shape_vertex_base_3<K> Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>   Fb;

typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds> Triangulation_3;

typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

typedef Alpha_shape_3::Cell  Cell;
typedef Alpha_shape_3::Vertex Vertex;
typedef Alpha_shape_3::Edge Edge;
typedef Alpha_shape_3::Cell_handle  Cell_handle;
typedef Alpha_shape_3::Vertex_handle Vertex_handle;

typedef Alpha_shape_3::Cell_circulator  Cell_circulator;

typedef Alpha_shape_3::Locate_type Locate_type;

typedef Alpha_shape_3::Finite_cells_iterator  Finite_cells_iterator;
typedef Alpha_shape_3::Finite_facets_iterator Finite_facets_iterator;
typedef Alpha_shape_3::Finite_edges_iterator  Finite_edges_iterator;
typedef Alpha_shape_3::Finite_vertices_iterator 
                                              Finite_vertices_iterator;
typedef Alpha_shape_3::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_3::Edge_iterator  Edge_iterator;


typedef Alpha_shape_3::Coord_type Coord_type;
typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;
typedef Alpha_shape_3::Alpha_shape_cells_iterator 
                                      Alpha_shape_cells_iterator;
typedef Alpha_shape_3::Alpha_shape_vertices_iterator 
                                      Alpha_shape_vertices_iterator;
typedef Alpha_shape_3::Alpha_shape_facets_iterator 
                                      Alpha_shape_facets_iterator;

//---------------------------------------------------------------------

void
construct_alpha_shape(const std::list<Point> &V_p,
		      Alpha_shape_3::Mode mode,
		      Alpha_shape_3& A)
  // Generate Alpha Shape
{ 
  //  std::vector<K::Segment_3> V_seg;
  
  int  n = A.make_alpha_shape(V_p.begin(), V_p.end());
  std::cout << "Inserted " << n  << " points" << std::endl;
  
  A.set_mode(mode);
}

bool
file_input(std::list<Point>& L)
{

  std::ifstream is("./data/fin", std::ios::in);

  if(is.fail())
    {
      std::cerr << "unable to open file for input" << std::endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  Point p;
  for( ; n>0 ; n--)
    {
      is >> p;
      L.push_back(p);
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}

void count_faces(const Alpha_shape_3 &A)
{
  Alpha_shape_cells_iterator cit=A.alpha_shape_cells_begin();
  int count_cells=0;
  for ( ; cit != A.alpha_shape_cells_end() ; ++cit) ++count_cells;
  Finite_cells_iterator fcit=A.finite_cells_begin();
  int count_cint=0;
  for (; fcit != A.finite_cells_end(); ++fcit) 
    if ( A.classify(fcit) ==  Alpha_shape_3::INTERIOR) ++count_cint;
  assert(count_cells == count_cint);

  int count_exterior_facets = 0;
  int count_singular_facets = 0;
  int count_regular_facets  = 0;
  int count_interior_facets = 0;
  int count_facets = 0;
  Alpha_shape_facets_iterator face_iterator = A.alpha_shape_facets_begin();
  for (;face_iterator!=A.alpha_shape_facets_end();face_iterator++) 
    count_facets++;
  Finite_facets_iterator fit = A.finite_facets_begin();
  for ( ;  fit != A.finite_facets_end(); ++fit) {
    switch(A.classify(*fit) ) {
    case Alpha_shape_3::EXTERIOR : ++count_exterior_facets; break;
    case Alpha_shape_3::SINGULAR : ++count_singular_facets;  break;
    case Alpha_shape_3::REGULAR  : ++count_regular_facets;  break;
    case Alpha_shape_3::INTERIOR : ++count_interior_facets;  break;
    }
  }
  
  if (A.get_mode() == Alpha_shape_3::REGULARIZED)
    assert( count_facets == count_regular_facets );
  if( A.get_mode() == Alpha_shape_3::GENERAL)
    assert( count_facets == count_regular_facets + count_singular_facets);
   
  int count_exterior_edges = 0;
  int count_singular_edges = 0;
  int count_regular_edges = 0;
  int count_interior_edges = 0;
  Triangulation_3::Finite_edges_iterator e_iterator;
  for (e_iterator = A.finite_edges_begin();
       e_iterator!=A.finite_edges_end();e_iterator++){
       switch(A.classify(*e_iterator)) {
       case Alpha_shape_3::EXTERIOR : ++count_exterior_edges; break;
       case Alpha_shape_3::SINGULAR : ++count_singular_edges;  break;
       case Alpha_shape_3::REGULAR  : ++count_regular_edges;  break;
       case Alpha_shape_3::INTERIOR : ++count_interior_edges;  break;
    }
  }

  int count_regular_vertices = 0;
  int count_singular_vertices = 0;
  int count_exterior_vertices = 0;
  int count_interior_vertices = 0;
  int count_vertices = 0;
  Alpha_shape_vertices_iterator vit=A.alpha_shape_vertices_begin();
  for ( ; vit != A.alpha_shape_vertices_end() ; ++vit) ++count_vertices;
  Finite_vertices_iterator fvit = A.finite_vertices_begin();
  for ( ; fvit != A.finite_vertices_end(); ++fvit) {
    switch(A.classify(fvit)) {
    case Alpha_shape_3::EXTERIOR : ++count_exterior_vertices; break;
    case Alpha_shape_3::SINGULAR : ++count_singular_vertices;  break;
    case Alpha_shape_3::REGULAR  : ++count_regular_vertices;  break;
    case Alpha_shape_3::INTERIOR : ++count_interior_vertices;  break;
    }
  }
  if (A.get_mode() == Alpha_shape_3::REGULARIZED)
    assert(count_vertices == count_regular_vertices );
  if( A.get_mode() == Alpha_shape_3::GENERAL)
    assert( count_facets == count_regular_facets + count_singular_facets);
	                  
  int ncc = A.number_of_solid_components();

 //  std::cout << "alpha " << A.get_alpha() << std::endl;
//   std::cout << "The current alpha shape has " << std::endl;
//   std::cout << "\t" << count_cells << " cells" << std::endl
// 	    << std::endl;

//   std::cout << "\t" << count_singular_facets << " singular_facets " 
// 	    << std::endl;
//   std::cout << "\t" << count_regular_facets  << " regular_facets " 
// 	    << std::endl << std::endl;

//   std::cout << "\t" << count_exterior_edges << " exterior_edges" << std::endl;
//   std::cout << "\t" << count_singular_edges << " singular edges" <<  std::endl;
//   std::cout << "\t" << count_regular_edges  << " regular edges" << std::endl;
//   std::cout << "\t" << count_interior_edges  << " interior edges" << std::endl;  
//   std::cout << "\t" << count_regular_vertices << " regular vertices" 
// 	    << std::endl;
//   std::cout << "\t" << count_singular_vertices << " singular vertices" 
// 	    << std::endl;
//   std::cout << "number_of_components " 
// 	    << A.number_of_solid_components() << std::endl;

 if (count_cells >= 1) {
     assert(count_regular_facets == 2*count_regular_vertices - 4*ncc);
    //this relation might not be valid for any alpha_shape
    assert(count_regular_edges == 3*count_regular_vertices - 6*ncc);
  }
}


//TO DEBUG
// void show_triangulation(Alpha_shape_3& A)
// {
//   std::cout << "number of finite cells " << A.number_of_finite_cells() 
// 	    << std::endl;
//   std::cout << "number of finite faces " << A.number_of_finite_facets() 
// 	    << std::endl;
//   std::cout << "number of finite edges " << A.number_of_finite_edges() 
// 	    << std::endl;

//   Triangulation_3::Finite_cells_iterator 
//     cit = A.finite_cells_begin();
//   for ( ; cit != A.finite_cells_end(); ++cit) {
//     std::cout << "cell " <<  "alpha " << cit->get_alpha() << std::endl;
//     std::cout << cit->vertex(0)->point() << std::endl
// 	      << cit->vertex(1)->point() << std::endl
// 	      << cit->vertex(2)->point() << std::endl
// 	      << cit->vertex(3)->point() << std::endl;
//   }
// }

// void explore_alpha_shape(Alpha_shape_3& A )
// {
//   A.show_interval_facet_map();

//   Alpha_iterator alpha_it = A.alpha_begin();
//   std::cout<<"The values of alpha are:\n\n";
//   for (;alpha_it!=A.alpha_end();alpha_it++){
//     std::cout<<"\t"<<*alpha_it<<std::endl;
//   }
// }


//------------------ main -------------------------------------------

int main()
{
  std::list<Point> L;
  
  // first a known small case
  // a cube with one corner less and two small pyramids
  // on back and front face
  L.push_back(Point(0.,0.,0.));
  L.push_back(Point(2.,0.,0.));
  L.push_back(Point(0.,2.,0.));
  L.push_back(Point(0.,0.,2.));
  L.push_back(Point(2.,0.,2.));
  L.push_back(Point(0.,2.,2.));
  L.push_back(Point(0.,2.,2.));
  L.push_back(Point(2.,2.,0));
  L.push_back(Point(1.,-1.5 ,1.));
  L.push_back(Point(1.,3.5 ,1.));
  
  Alpha_shape_3 A( L.begin(), L.end(), 0, Alpha_shape_3::REGULARIZED);
  assert(A.number_of_alphas() == 6 );
  std::cout << "test_classify_and_iterators in regularised mode" 
	    << std::endl;
  Alpha_iterator alpha_it = A.alpha_begin();
  for (;alpha_it!=A.alpha_end();alpha_it++){
    A.set_alpha(*alpha_it);
    count_faces(A);
  }

  A.set_mode(Alpha_shape_3::GENERAL);
  assert(A.number_of_alphas() == 6) ;
	 std::cout << "test_classify_and_iterators in general mode" 
	    << std::endl;
 
  for(alpha_it = A.alpha_begin();alpha_it!=A.alpha_end();alpha_it++){
    A.set_alpha(*alpha_it);
    count_faces(A);
  }
// alpha values 1.38942
//              2
//              2.00694
//              2.66667
//              3
//              4.35417
  std::cout << "test number_of_components - find_optimal_alpha "<< std::endl;
  assert( *(A.find_optimal_alpha(1)) == A.get_nth_alpha(5));
  assert( *(A.find_optimal_alpha(2)) == A.get_nth_alpha(3));
  assert (A.number_of_solid_components(*(A.find_optimal_alpha(2))) == 2);
  assert (A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);
  assert (A.get_nth_alpha(1) == *(A.alpha_lower_bound(2)));
  assert (A.get_nth_alpha(2) == *(A.alpha_upper_bound(2)));

  // test a bigger alpha_shapes
  A.clear();
  L.clear();
  file_input(L);  
  construct_alpha_shape(L,Alpha_shape_3::GENERAL,A);
  std::cout << "Alpha Shape computed" << std::endl;
  A.set_alpha(*A.find_optimal_alpha(2));
  std::cout << " test number_of_components - find_optimal_alpha "<< std::endl;
  assert( A.number_of_solid_components() == 2);

 

  return 0;
}
