// standard includes

//#define CGAL_SDG_VERBOSE 

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

// define the kernel
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

typedef CGAL::Simple_cartesian<double>    CK;
typedef CGAL::Filtered_kernel<CK>         Kernel;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>             SDG2;

using namespace std;

int main( int argc, char *argv[] ) {
  if ( not (( argc == 1 ) or (argc == 2)) ) {
    std::cout <<"usage: "<< argv[0] <<" [filename]\n";
  }

  ifstream ifs( (argc == 1) ? "data/sites2.cin" : argv[1] );
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site;

  // read the sites from the stream and insert them in the diagram
  while ( ifs >> site ) { sdg.insert( site ); }

  ifs.close();


  std::cout << "About to validate diagram ..." << std::endl;

  // validate the diagram
  //assert( sdg.is_valid(true, 1) );
  cout << endl << endl;

  std::cout << "Diagram validated." << std::endl;

  string inf_vertex("infinite vertex");
  char vid[] = {'A', 'B', 'C'};

  cout << "Sandeep: the number of vertices in sdg = " << sdg.number_of_vertices() << endl;
  cout << "Sandeep: the number of faces in sdg = " << sdg.number_of_faces() << endl;
  cout << "Sandeep: the number of input sites in sdg = " << sdg.number_of_input_sites() << endl;
  cout << "Sandeep: the number of output sites in sdg = " << sdg.number_of_output_sites() << endl;
  // print the vertices of the segment Delaunay graph
  cout << std::endl;
  cout << "Output sites/ vertices of sdg:" << std::endl;
  cout << "------------------------------" << std::endl;
  SDG2::Finite_vertices_iterator     fvit;
  SDG2::All_vertices_iterator avit;
  SDG2::Vertex_handle vh;
  
  int v_count=0;
  for (avit = sdg.all_vertices_begin();
       avit != sdg.all_vertices_end(); ++avit) {
    vh = avit;
    if (sdg.is_infinite(vh)) {
      cout << "vertex " << ++v_count << " : " << inf_vertex << endl;
    } else {
      cout << "vertex " << ++v_count << " : " << avit->site() << endl;
    }
  }
  
  cout << std::endl;
  cout << "Edges of sdg:" << std::endl;
  cout << "-------------" << std::endl;
  
  SDG2::All_edges_iterator eit = sdg.all_edges_begin();
  for (int k = 1; eit != sdg.all_edges_end(); ++eit, ++k) {
    SDG2::Edge e = *eit;
    // get the vertices defining the edge
    SDG2::Vertex_handle v[] = { e.first->vertex(sdg.ccw(e.second)),
      e.first->vertex(sdg.cw(e.second))
    };
    
    cout << "--- Edge " << k << " ---" << endl;
    for (int i = 0; i < 2; i++) {
      // check if the vertex is the vertex at infinity; if yes, print
      // the corresponding string, otherwise print the site
      if ( sdg.is_infinite(v[i]) ) {
        cout << vid[i] << ": " << inf_vertex << endl;
      } else {
        cout << vid[i] << ": " << v[i]->site() << endl;
      }
    }
    cout << endl;
  }
  
  cout << std::endl;
  cout << "Faces of sdg:" << std::endl;
  cout << "-------------" << std::endl;
  
  SDG2::All_faces_iterator fit = sdg.all_faces_begin();
  for (int k = 1; fit != sdg.all_faces_end(); ++fit, ++k) {
    SDG2::Face f = *fit;
    // get the vertices defining the face
    SDG2::Vertex_handle v[] = { f.vertex(0),
      f.vertex(1),
      f.vertex(2)
    };
    
    cout << "--- Face " << k << " ---" << endl;
    for (int i = 0; i < 3; i++) {
      // check if the vertex is the vertex at infinity; if yes, print
      // the corresponding string, otherwise print the site
      if ( sdg.is_infinite(v[i]) ) {
        cout << vid[i] << ": " << inf_vertex << endl;
      } else {
        cout << vid[i] << ": " << v[i]->site() << endl;
      }
    }
    cout << endl;
  }
  //Sandeep: A counterclockwise traversal of the vertices adjacent to the infinite_vertex
  //is a clockwise traversal of the convex hull.
  cout << std::endl;
  cout << "Convex-hull of sdg:" << std::endl;
  cout << "-------------------" << std::endl;
  SDG2::Vertex_handle v_inf = sdg.infinite_vertex();
  SDG2::Vertex_circulator	 vc1 = sdg.incident_vertices(v_inf);
  SDG2::Vertex_circulator	 vc2 = vc1;
  if (sdg.is_infinite(v_inf)){
    cout << "vertex 0 : " << inf_vertex << endl;
  }
  int cnt = 0;
  if (vc1 != 0) {
    do {
      vh = vc1;
      if (sdg.is_infinite(vh)) {
        cout << "vertex " << ++cnt << " : " << inf_vertex << endl;
      } else {
        cout << "vertex " << ++cnt << " : " << vc1->site() << endl;
      }
    } while (++vc1 != vc2);
  }
  
  //Write the current state of the segment Delaunay graph to an output stream.
  //ostream os;
  //sdg.file_output (ostream& os);
  cout << sdg ;
  
  return 0;
}
