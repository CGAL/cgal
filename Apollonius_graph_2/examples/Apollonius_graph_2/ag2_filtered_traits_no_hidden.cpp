// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits

// choose the representation
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Rep;

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>

// typedef for the traits; the filtered traits class is used
typedef CGAL::Apollonius_graph_filtered_traits_2<Rep> Traits;

// typedefs for the algorithm

// With the second template argument in the vertex base class being
// false, we indicate that there is no need to store the hidden sites.
// One case where this is indeed not needed is when we only do
// insertions, like in the main program below.
typedef CGAL::Apollonius_graph_vertex_base_2<Traits,false>   Vb;
typedef CGAL::Triangulation_face_base_2<Traits>              Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       Agds;
typedef CGAL::Apollonius_graph_2<Traits,Agds>    Apollonius_graph;


int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Site_2 site;

  // read the sites and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  // now remove all sites
  std::cout << "Removing all sites... " << std::flush;
  while ( ag.number_of_vertices() > 0 ) {
    ag.remove( ag.finite_vertex() );
  }
  std::cout << "done!" << std::endl << std::endl;

  return 0;
}
