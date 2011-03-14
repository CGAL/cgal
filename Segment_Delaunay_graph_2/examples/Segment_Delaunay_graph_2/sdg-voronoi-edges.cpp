// standard includes
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
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>

typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>             SDG2;

using namespace std;

int main()
{
  ifstream ifs("data/sites2.cin");
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site;

  // read the sites from the stream and insert them in the diagram
  while ( ifs >> site ) { sdg.insert( site ); }

  ifs.close();

  // validate the diagram
  assert( sdg.is_valid(true, 1) );
  cout << endl << endl;

  /*
  // now walk through the non-infinite edges of the segment Delaunay
  // graphs (which are dual to the edges in the Voronoi diagram) and
  // print the sites defining each Voronoi edge.
  //
  // Each oriented Voronoi edge (horizontal segment in the figure
  // below) is defined by four sites A, B, C and D.
  //
  //     \                     /
  //      \         B         /
  //       \                 /
  //     C  -----------------  D
  //       /                 \
  //      /         A         \
  //     /                     \
  //
  // The sites A and B define the (oriented) bisector on which the
  // edge lies whereas the sites C and D, along with A and B define
  // the two endpoints of the edge. These endpoints are the Voronoi
  // vertices of the triples A, B, C and B, A, D.
  // If one of these vertices is the vertex at infinity the string
  // "infinite vertex" is printed; the corresponding Voronoi edge is
  // actually a stright-line or parabolic ray.
  // The sites below are printed in the order A, B, C, D.
  */

  string inf_vertex("infinite vertex");
  char vid[] = {'A', 'B', 'C', 'D'};

  SDG2::Finite_edges_iterator eit = sdg.finite_edges_begin();
  for (int k = 1; eit != sdg.finite_edges_end(); ++eit, ++k) {
    SDG2::Edge e = *eit;
    // get the vertices defining the Voronoi edge
    SDG2::Vertex_handle v[] = { e.first->vertex( sdg.ccw(e.second) ),
                                e.first->vertex( sdg.cw(e.second) ),
                                e.first->vertex( e.second ),
                                sdg.tds().mirror_vertex(e.first, e.second) };

    cout << "--- Edge " << k << " ---" << endl;
    for (int i = 0; i < 4; i++) {
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

  return 0;
}
