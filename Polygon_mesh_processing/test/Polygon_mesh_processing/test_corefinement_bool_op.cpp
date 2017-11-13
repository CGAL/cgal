// #define CGAL_COREFINEMENT_POLYHEDRA_DEBUG
// #define CGAL_COREFINEMENT_DEBUG
// #define CGAL_TODO_WARNINGS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

// includes for Operations on polyhedra
#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/internal/corefinement/Polyhedra_output_builder.h>

#include <CGAL/iterator.h>
#include <CGAL/array.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace CFR = CGAL::Corefinement;

//includes typedefs for Operations on polyhedra
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;
typedef std::map<Polyhedron::Facet_const_handle,std::size_t>       Facet_id_map;
typedef boost::associative_property_map<Facet_id_map>             Facet_id_pmap;
typedef CGAL::Corefinement
            ::Polyhedra_output_builder< Polyhedron,
                                        Facet_id_pmap>           Output_builder;
typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron,
                                            Output_builder&>       Split_visitor;

struct Result_checking
{
  bool check;
  bool union_res;
  bool inter_res;
  bool P_minus_Q_res;
  bool Q_minus_P_res;

  Result_checking() : check(false) {}
};

void run_boolean_operations(
  Polyhedron& P,
  Polyhedron& Q,
  Polyhedron& union_,
  Polyhedron& inter,
  Polyhedron& P_minus_Q,
  Polyhedron& Q_minus_P,
  std::string scenario,
  std::size_t id,
  const Result_checking& rc)
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
   assert(union_.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_P_union_Q.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << union_;
   #endif
  }
  else
    std::cout << "  Union is invalid\n";

  if ( output_builder.intersection_valid() ){
   std::cout << "  Intersection is a valid operation\n";
   assert(inter.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_P_inter_Q.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << inter;
   #endif
  }
  else
    std::cout << "  Intersection is invalid\n";

  if ( output_builder.P_minus_Q_valid() ){
   std::cout << "  P-Q is a valid operation\n";
   assert(P_minus_Q.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_P_minus_Q.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << P_minus_Q;
   #endif
  }
  else
    std::cout << "  P-Q is invalid\n";

  if ( output_builder.Q_minus_P_valid() ){
   std::cout << "  Q-P is a valid operation\n";
   assert(Q_minus_P.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_Q_minus_P.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << Q_minus_P;
   #endif
  }
  else
    std::cout << "  Q-P is invalid\n";

  if (rc.check)
  {
    if(output_builder.union_valid()!=rc.union_res ||
       output_builder.intersection_valid()!=rc.inter_res ||
       output_builder.P_minus_Q_valid()!=rc.P_minus_Q_res ||
       output_builder.Q_minus_P_valid()!=rc.Q_minus_P_res)
    {
      std::cerr << "  ERROR: at least one operation is not as expected!\n";
      exit(EXIT_FAILURE);
    }
    else
      std::cout << "  All operations are as expected.\n";
  }
}


void run_boolean_operations(
  Surface_mesh& tm1,
  Surface_mesh& tm2,
  Surface_mesh& union_,
  Surface_mesh& inter,
  Surface_mesh& tm1_minus_tm2,
  Surface_mesh& tm2_minus_tm1,
  std::string scenario,
  std::size_t id,
  const Result_checking& rc)
{
  std::cout << "Scenario #" << id << " - " << scenario << "\n";

  typedef boost::optional<Surface_mesh*> OSM;

  CGAL::cpp11::array<OSM,4> desired_output;
  
  desired_output[CFR::UNION]=OSM( &union_ );
  desired_output[CFR::INTER]=OSM( &inter );
  desired_output[CFR::TM1_MINUS_TM2]=OSM( &tm1_minus_tm2 );
  desired_output[CFR::TM2_MINUS_TM1]=OSM( &tm2_minus_tm1 );

  std::cout << "  Vertices before " <<  tm1.number_of_vertices()
            << " " << tm2.number_of_vertices() << std::endl;

  CGAL::cpp11::array<bool,4> res = PMP::boolean_operation(tm1, tm2, desired_output);

  std::cout << "  Vertices after " <<  tm1.number_of_vertices()
            << " " << tm2.number_of_vertices() << std::endl;

  if ( res[CFR::UNION] ){
   std::cout << "  Union is a valid operation\n";
   assert(union_.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_tm1_union_tm2.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << union_;
   #endif
  }
  else
    std::cout << "  Union is invalid\n";

  if ( res[CFR::INTER] ){
   std::cout << "  Intersection is a valid operation\n";
   assert(inter.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_tm1_inter_tm2.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << inter;
   #endif
  }
  else
    std::cout << "  Intersection is invalid\n";

  if ( res[CFR::TM1_MINUS_TM2] ){
   std::cout << "  tm1-tm2 is a valid operation\n";
   assert(tm1_minus_tm2.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_tm1_minus_tm2.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << tm1_minus_tm2;
   #endif
  }
  else
    std::cout << "  tm1-tm2 is invalid\n";

  if ( res[CFR::TM2_MINUS_TM1] ){
   std::cout << "  tm2-tm1 is a valid operation\n";
   assert(tm2_minus_tm1.is_valid());
   #ifdef CGAL_COREFINEMENT_DEBUG
   std::stringstream fname;
   fname << scenario << "_tm2_minus_tm1.off";
   std::ofstream output(fname.str().c_str());
   output << std::setprecision(17) << tm2_minus_tm1;
   #endif
  }
  else
    std::cout << "  tm2-tm1 is invalid\n";

  if ( rc.check )
  {
    if (res[CFR::UNION]!=rc.union_res ||
        res[CFR::INTER]!=rc.inter_res ||
        res[CFR::TM1_MINUS_TM2]!=rc.P_minus_Q_res ||
        res[CFR::TM2_MINUS_TM1]!=rc.Q_minus_P_res)
    {
      std::cerr << "  ERROR: at least one operation is not as expected!\n";
      exit(EXIT_FAILURE);
    }
    else
      std::cout << "  All operations are as expected.\n";
  }
}

template <class TriangleMesh>
void run(char* P_fname, char* Q_fname, int k,
         const Result_checking& rc)
{
  TriangleMesh P, Q, W, X, Y, Z;
  std::vector< CGAL::cpp11::array<TriangleMesh*,4> > scenarios;
  std::vector< std::string > scenarios_str;
  scenarios.reserve(21); // 21 = An2 + 2* An1 + An0
  scenarios_str.reserve(21);

  // no P nor Q
  scenarios.push_back( CGAL::make_array(&W, &X, &Y, &Z) );
  scenarios_str.push_back( "NNNN" ); // #0
  // P for union
  scenarios.push_back( CGAL::make_array(&P, &Q, &X, &Y) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Q, &Y) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Y, &Q) );
  scenarios.push_back( CGAL::make_array(&P, &X, &Y, &Z) );
  scenarios.push_back( CGAL::make_array(&Q, &X, &Y, &Z) );
  scenarios_str.push_back( "PQNN" ); // #1
  scenarios_str.push_back( "PNQN" ); // #2
  scenarios_str.push_back( "PNNQ" ); // #3
  scenarios_str.push_back( "PNNN" ); // #4
  scenarios_str.push_back( "QNNN" ); // #5
  // P for intersection
  scenarios.push_back( CGAL::make_array(&Q, &P, &X, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Q, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Y, &Q) );
  scenarios.push_back( CGAL::make_array(&X, &P, &Y, &Z) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &Y, &Z) );
  scenarios_str.push_back( "QPNN" ); // #6
  scenarios_str.push_back( "NPQN" ); // #7
  scenarios_str.push_back( "NPNQ" ); // #8
  scenarios_str.push_back( "NPNN" ); // #9
  scenarios_str.push_back( "NQNN" ); // #10
  // P for P-Q
  scenarios.push_back( CGAL::make_array(&Q, &X, &P, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &P, &Y) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &P, &Q) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &P, &Z) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Q, &Z) );
  scenarios_str.push_back( "QNPN" ); // #11
  scenarios_str.push_back( "NQPN" ); // #12
  scenarios_str.push_back( "NNPQ" ); // #13
  scenarios_str.push_back( "NNPN" ); // #14
  scenarios_str.push_back( "NNQN" ); // #15
  // P for Q-P
  scenarios.push_back( CGAL::make_array(&Q, &X, &Y, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Q, &Y, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Q, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Z, &P) );
  scenarios.push_back( CGAL::make_array(&X, &Y, &Z, &Q) );
  scenarios_str.push_back( "QNNP" ); // #16
  scenarios_str.push_back( "NQNP" ); // #17
  scenarios_str.push_back( "NNQP" ); // #18
  scenarios_str.push_back( "NNNP" ); // #19
  scenarios_str.push_back( "NNNQ" ); // #20

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
      scenarios_str[k], k, rc);
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
        scenarios_str[i], i, rc );
    }
}

int main(int argc,char** argv)
{
  if (argc<3){
    std::cerr << "Usage "<< argv[0] << " file1.off file2.off [SM/POLY] [scenario_id/ALL] [0/1 0/1 0/1 0/1 (expected valid operations U I P-Q Q-P)]\n";
    return 1;
  }

  Result_checking rc;

  if (argc==9)
  {
    rc.check=true;
    rc.union_res = atoi(argv[5])!=0;
    rc.inter_res = atoi(argv[6])!=0;
    rc.P_minus_Q_res = atoi(argv[7])!=0;
    rc.Q_minus_P_res = atoi(argv[8])!=0;
  }

  int scenario = -1;
  if (argc>=5 && std::string(argv[4])!=std::string("ALL"))
    scenario = atoi(argv[4]);

  if ( argc==3 || std::string(argv[3])==std::string("SM") )
    run<Surface_mesh>(argv[1], argv[2], scenario, rc);
  else
    if ( std::string(argv[3])==std::string("POLY") )
      run<Polyhedron>(argv[1], argv[2], scenario, rc);
    else
    {
      std::cout <<"Invalid value, only SM or POLY can be given\n";
      return 1;
    }
  return 0;
}
