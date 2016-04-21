#define CGAL_COREFINEMENT_POLYHEDRA_DEBUG
// #define CGAL_COREFINEMENT_DEBUG

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
#include <CGAL/array.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;
typedef std::map<Polyhedron::Facet_const_handle,std::size_t>       Facet_id_map;
typedef boost::associative_property_map<Facet_id_map>             Facet_id_pmap;
typedef CGAL::Corefinement
            ::Polyhedra_output_builder< Polyhedron,
                                        Facet_id_pmap>           Output_builder;
typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron,
                                            Output_builder>       Split_visitor;


void run_boolean_operations(
  Polyhedron& P,
  Polyhedron& Q,
  Polyhedron& union_,
  Polyhedron& inter,
  Polyhedron& P_minus_Q,
  Polyhedron& Q_minus_P,
  std::string scenario,
  int id)
{
  std::cout << "Scenario #" << id << " - " << scenario << "\n";

  CGAL::Emptyset_iterator output_it;
  Facet_id_map P_facet_id_map, Q_facet_id_map;

  CGAL::cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output;

  desired_output[Output_builder::P_UNION_Q]=boost::make_optional( &union_ );
  desired_output[Output_builder::P_INTER_Q]=boost::make_optional( &inter );
  desired_output[Output_builder::P_MINUS_Q]=boost::make_optional( &P_minus_Q );
  desired_output[Output_builder::Q_MINUS_P]=boost::make_optional( &Q_minus_P );

  Output_builder output_builder(P, Q,
                                desired_output,
                                Facet_id_pmap(P_facet_id_map),
                                Facet_id_pmap(Q_facet_id_map) );
  Split_visitor visitor(output_builder);

  CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections(visitor);
  std::cout << "  Vertices before " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;
  polyline_intersections(P, Q, output_it);
  std::cout << "  Vertices after " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;


  if ( output_builder.union_valid() ){
   std::cout << "  Union is a valid operation\n";
   CGAL_assertion(union_.is_valid());
   std::stringstream fname;
   fname << scenario << "_P_union_Q.off";
   std::ofstream output(fname.str().c_str());
   output << union_;
  }
  else
    std::cout << "  Union is invalid\n";

  if ( output_builder.intersection_valid() ){
   std::cout << "  Intersection is a valid operation\n";
   CGAL_assertion(inter.is_valid());
   std::stringstream fname;
   fname << scenario << "_P_inter_Q.off";
   std::ofstream output(fname.str().c_str());
   output << inter;
  }
  else
    std::cout << "  Intersection is invalid\n";

  if ( output_builder.P_minus_Q_valid() ){
   std::cout << "  P-Q is a valid operation\n";
   CGAL_assertion(P_minus_Q.is_valid());
   std::stringstream fname;
   fname << scenario << "_P_minus_Q.off";
   std::ofstream output(fname.str().c_str());
   output << P_minus_Q;
  }
  else
    std::cout << "  P-Q is invalid\n";

  if ( output_builder.Q_minus_P_valid() ){
   std::cout << "  Q-P is a valid operation\n";
   CGAL_assertion(Q_minus_P.is_valid());
   std::stringstream fname;
   fname << scenario << "_Q_minus_P.off";
   std::ofstream output(fname.str().c_str());
   output << Q_minus_P;
  }
  else
    std::cout << "  Q-P is invalid\n";
}

void run(char* P_fname, char* Q_fname, int k=-1)
{
  Polyhedron P, Q, W, X, Y, Z;
  std::vector< CGAL::cpp11::array<Polyhedron*,4> > scenarios;
  std::vector< std::string > scenarios_str;
  scenarios.reserve(21); // 21 = An2 + 2* An1 + An0
  scenarios_str.reserve(21);

  // no P nor Q
  scenarios.push_back( CGAL::make_array(&W, &X, &Y, &Z) );
  scenarios_str.push_back( "NNNN" );
  // P for union
  scenarios.push_back( CGAL::make_array(&P, &Q, &X, &Y) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Q, &Y) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Y, &Q) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Y, &Z) );
  scenarios.push_back( CGAL::make_array(&Q, &X, &Y, &Z) );
  scenarios_str.push_back( "PQNN" );
  scenarios_str.push_back( "PNQN" );
  scenarios_str.push_back( "PNNQ" );
  scenarios_str.push_back( "PNNN" );
  scenarios_str.push_back( "QNNN" );
  // P for intersection
  scenarios.push_back( CGAL::make_array(&Q, &P, &X, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Q, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Y, &Q) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Y, &Z) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &Y, &Z) );
  scenarios_str.push_back( "QPNN" );
  scenarios_str.push_back( "NPQN" );
  scenarios_str.push_back( "NPNQ" );
  scenarios_str.push_back( "NPNN" );
  scenarios_str.push_back( "NQNN" );
  // P for P-Q
  scenarios.push_back( CGAL::make_array(&Q, &X, &P, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &P, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &P, &Q) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &P, &Z) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Q, &Z) );
  scenarios_str.push_back( "QNPN" );
  scenarios_str.push_back( "NQPN" );
  scenarios_str.push_back( "NNPQ" );
  scenarios_str.push_back( "NNPN" );
  scenarios_str.push_back( "NNQN" );
  // P for Q-P
  scenarios.push_back( CGAL::make_array(&Q, &X, &Y, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &Y, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Q, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Z, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Z, &Q) );
  scenarios_str.push_back( "QNNP" );
  scenarios_str.push_back( "NQNP" );
  scenarios_str.push_back( "NNQP" );
  scenarios_str.push_back( "NNNP" );
  scenarios_str.push_back( "NNNQ" );

  if (k!=-1)
  {
    std::ifstream file(P_fname);
    file >> P;
    file.close();
    file.open(Q_fname);
    file >> Q;
    file.close();

    run_boolean_operations(
      P, Q,
      *scenarios[k][0], *scenarios[k][1], *scenarios[k][2], *scenarios[k][3],
      scenarios_str[k], k );
  }
  else
    for (std::size_t i=0; i<scenarios.size();++i)
    {
      P.clear(); Q.clear(); W.clear(); X.clear(); Y.clear(); Z.clear();
      std::ifstream file(P_fname);
      file >> P;
      file.close();
      file.open(Q_fname);
      file >> Q;
      file.close();

      run_boolean_operations(
        P, Q,
        *scenarios[i][0], *scenarios[i][1], *scenarios[i][2], *scenarios[i][3],
        scenarios_str[i], i );
    }
}

int main(int argc,char** argv)
{
  if (argc<3){
    std::cerr << "Usage "<< argv[0] << " file1.off file2.off [scenario_id]\n";
    return 1;
  }
  if (argc==3)
    run(argv[1], argv[2]);
  else
    run(argv[1], argv[2], atoi(argv[3]));
  return 0;
}
