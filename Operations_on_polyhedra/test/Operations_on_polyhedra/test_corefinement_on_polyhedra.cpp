#define CGAL_COREFINEMENT_POLYHEDRA_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#define CGAL_TODO_WARNINGS

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/internal/corefinement/Polyhedra_output_builder.h>
#include <CGAL/iterator.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;
typedef std::map<Polyhedron::Facet_const_handle,std::size_t>       Facet_id_map;
typedef boost::associative_property_map<Facet_id_map>             Facet_id_pmap;
typedef CGAL::Corefinement
            ::Polyhedra_output_builder< Polyhedron,
                                        Facet_id_pmap>           Output_builder;
typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron,
                                            Output_builder>       Split_visitor;

int main(int argc,char** argv) {
  if (argc!=3){
    std::cerr << "Usage "<< argv[0] << " file1.off file2.off\n";
    return 1;
  }

  Polyhedron P, Q;
  std::ifstream file(argv[1]);
  file >> P;
  file.close();
  file.open(argv[2]);
  file >> Q;
  file.close();

  CGAL::set_pretty_mode(std::cerr);

  CGAL::Emptyset_iterator output_it;
  Facet_id_map P_facet_id_map, Q_facet_id_map;

  CGAL::cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output;
  Polyhedron inter, union_;
  desired_output[Output_builder::P_MINUS_Q]=boost::make_optional( &P );
  desired_output[Output_builder::Q_MINUS_P]=boost::make_optional( &Q );
  desired_output[Output_builder::P_INTER_Q]=boost::make_optional( &inter );
  desired_output[Output_builder::P_UNION_Q]=boost::make_optional( &union_ );
  
  Output_builder output_builder(P, Q,
                                desired_output,
                                Facet_id_pmap(P_facet_id_map),
                                Facet_id_pmap(Q_facet_id_map) );
  Split_visitor visitor(output_builder);

  CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections(visitor);
  std::cout << "Vertices before " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;
  polyline_intersections(P, Q, output_it);
  std::cout << "Vertices after " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;

  if ( output_builder.union_valid() ) std::cout << "Union is valid\n";
  else std::cout << "Union is invalid\n";
  if ( output_builder.intersection_valid() ) std::cout << "Intersection is valid\n";
  else std::cout << "Intersection is invalid\n";
  if ( output_builder.P_minus_Q_valid() ) std::cout << "P-Q is valid\n";
  else std::cout << "P-Q is invalid\n";
  if ( output_builder.Q_minus_P_valid() ) std::cout << "Q-P is valid\n";
  else std::cout << "Q-P is invalid\n";

  CGAL_assertion(P.is_valid());
  CGAL_assertion(Q.is_valid());
  CGAL_assertion(inter.is_valid());
  CGAL_assertion(union_.is_valid());

  std::ofstream output("P_minus_Q.off");
  output << P;
  output.close();
  
  output.open("P_inter_Q.off");
  output << inter;
  output.close();
  
  output.open("P_union_Q.off");
  output << union_;
  output.close();
  
  output.open("Q_minus_P.off");
  output << Q;
  output.close();

  return 0;
}
