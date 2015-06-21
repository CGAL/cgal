#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include "PolyhedralSurf_rings.h"
#include "compute_normals.h"
#include <CGAL/Ridges.h>
#include <CGAL/Umbilics.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <fstream>
#include <cassert>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT                      FT;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;

typedef CGAL::Surface_mesh<Point_3> PolyhedralSurf;

typedef boost::graph_traits<PolyhedralSurf>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<PolyhedralSurf>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<PolyhedralSurf>::face_descriptor   face_descriptor;

typedef T_PolyhedralSurf_rings<PolyhedralSurf> Poly_rings;
typedef CGAL::Monge_via_jet_fitting<Kernel>    Monge_via_jet_fitting;
typedef Monge_via_jet_fitting::Monge_form      Monge_form;

typedef PolyhedralSurf::Property_map<vertex_descriptor,FT> VertexFT_property_map;
typedef PolyhedralSurf::Property_map<vertex_descriptor,Vector_3> VertexVector_property_map;
//RIDGES
typedef CGAL::Ridge_line<PolyhedralSurf> Ridge_line;
typedef CGAL::Ridge_approximation < PolyhedralSurf,
				    VertexFT_property_map,
				    VertexVector_property_map > Ridge_approximation;
//UMBILICS
typedef CGAL::Umbilic<PolyhedralSurf> Umbilic;
typedef CGAL::Umbilic_approximation < PolyhedralSurf,
				      VertexFT_property_map,
				      VertexVector_property_map > Umbilic_approximation;

//create property maps

PolyhedralSurf::Property_map<vertex_descriptor,FT> 
vertex_k1_pm, vertex_k2_pm,
  vertex_b0_pm, vertex_b3_pm,
  vertex_P1_pm, vertex_P2_pm;

PolyhedralSurf::Property_map<vertex_descriptor,Vector_3> vertex_d1_pm, vertex_d2_pm;

PolyhedralSurf::Property_map<face_descriptor,Vector_3> face2normal_pm;

// default fct parameter values and global variables
unsigned int d_fitting = 3;
unsigned int d_monge = 3;
unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
unsigned int nb_points_to_use = 0;//
CGAL::Ridge_order tag_order = CGAL::Ridge_order_3;
double umb_size = 2;
bool verbose = false;
unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

/* gather points around the vertex v using rings on the
   polyhedralsurf. the collection of points resorts to 3 alternatives:
   1. the exact number of points to be used
   2. the exact number of rings to be used
   3. nothing is specified
*/
template <typename VertexPointMap>
void gather_fitting_points(vertex_descriptor v,
			   std::vector<Point_3> &in_points,
			   Poly_rings& poly_rings,
                           VertexPointMap vpm)
{
  //container to collect vertices of v on the PolyhedralSurf
  std::vector<vertex_descriptor> gathered;
  //initialize
  in_points.clear();

  //OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
  //enough rings and discard some points of the last collected ring to
  //get the exact "nb_points_to_use"
  if ( nb_points_to_use != 0 ) {
    poly_rings.collect_enough_rings(v, nb_points_to_use, gathered);
    if ( gathered.size() > nb_points_to_use ) gathered.resize(nb_points_to_use);
  }
  else { // nb_points_to_use=0, this is the default and the option -p is not considered;
    // then option -a nb_rings is checked. If nb_rings=0, collect
    // enough rings to get the min_nb_points required for the fitting
    // else collect the nb_rings required
    if ( nb_rings == 0 )
      poly_rings.collect_enough_rings(v, min_nb_points, gathered);
    else poly_rings.collect_i_rings(v, nb_rings, gathered);
  }

  //store the gathered points
  std::vector<vertex_descriptor>::const_iterator
    itb = gathered.begin(), ite = gathered.end();
  CGAL_For_all(itb,ite) in_points.push_back(get(vpm,*itb));
}

/* Use the jet_fitting package and the class Poly_rings to compute
   diff quantities.
*/
void compute_differential_quantities(PolyhedralSurf& P, Poly_rings& poly_rings)
{
  //container for approximation points
  std::vector<Point_3> in_points;

  typedef boost::property_map<PolyhedralSurf,CGAL::vertex_point_t>::type VPM;
  VPM vpm = get(CGAL::vertex_point,P);

  //MAIN LOOP
  vertex_iterator vitb = P.vertices_begin(), vite = P.vertices_end();
  for (; vitb != vite; vitb++) {
    //initialize
    vertex_descriptor v = * vitb;
    in_points.clear();
    Monge_form monge_form;
    Monge_via_jet_fitting monge_fit;

    //gather points around the vertex using rings
    gather_fitting_points(v, in_points, poly_rings, vpm);

    //exit if the nb of points is too small
    if ( in_points.size() < min_nb_points )
      {std::cerr << "Too few points to perform the fitting" << std::endl; exit(1);}

    //For Ridges we need at least 3rd order info
    assert( d_monge >= 3);
    // run the main fct : perform the fitting
    monge_form = monge_fit(in_points.begin(), in_points.end(),
			   d_fitting, d_monge);

    //switch min-max ppal curv/dir wrt the mesh orientation
    const Vector_3 normal_mesh = computeFacetsAverageUnitNormal(P,v, face2normal_pm, Kernel());
    monge_form.comply_wrt_given_normal(normal_mesh);

    //Store monge data needed for ridge computations in property maps
    vertex_d1_pm[v] = monge_form.maximal_principal_direction();
    vertex_d2_pm[v] = monge_form.minimal_principal_direction();
    vertex_k1_pm[v] = monge_form.coefficients()[0];
    vertex_k2_pm[v] = monge_form.coefficients()[1];
    vertex_b0_pm[v] = monge_form.coefficients()[2];
    vertex_b3_pm[v] = monge_form.coefficients()[5];
    if ( d_monge >= 4) {
      //= 3*b1^2+(k1-k2)(c0-3k1^3)
      vertex_P1_pm[v] =
	3*monge_form.coefficients()[3]*monge_form.coefficients()[3]
	+(monge_form.coefficients()[0]-monge_form.coefficients()[1])
	*(monge_form.coefficients()[6]
	  -3*monge_form.coefficients()[0]*monge_form.coefficients()[0]
	  *monge_form.coefficients()[0]);
      //= 3*b2^2+(k2-k1)(c4-3k2^3)
      vertex_P2_pm[v] =
	3*monge_form.coefficients()[4]*monge_form.coefficients()[4]
	+(-monge_form.coefficients()[0]+monge_form.coefficients()[1])
	*(monge_form.coefficients()[10]
	  -3*monge_form.coefficients()[1]*monge_form.coefficients()[1]
	  *monge_form.coefficients()[1]);
    }
  } //END FOR LOOP
}


///////////////MAIN///////////////////////////////////////////////////////
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
int main(int argc, char *argv[])
#else
int main()
#endif
{
  std::string if_name, of_name;// of_name same as if_name with '/' -> '_'

  try {
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
   unsigned int int_tag;
   po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message.")
      ("input-file,f", po::value<std::string>(&if_name)->default_value("data/poly2x^2+y^2-0.062500.off"),
       "name of the input off file")
      ("degree-jet,d", po::value<unsigned int>(&d_fitting)->default_value(3),
       "degree of the jet,  3 <= degre-jet <= 4")
      ("degree-monge,m", po::value<unsigned int>(&d_monge)->default_value(3),
       "degree of the Monge rep, 3<= degree-monge <= degree-jet")
      ("nb-rings,a", po::value<unsigned int>(&nb_rings)->default_value(0),
       "number of rings to collect neighbors. 0 means collect enough rings to make appro possible a>=1 fixes the nb of rings to be collected")
      ("nb-points,p", po::value<unsigned int>(&nb_points_to_use)->default_value(0),
       "number of neighbors to use.  0 means this option is not considered, this is the default p>=1 fixes the nb of points to be used")
      ("ridge_order,t", po::value<unsigned int>(&int_tag)->default_value(3),
       "Order of differential quantities used, must be 3 or 4")
      ("umbilic-patch-size,u", po::value<double>(&umb_size)->default_value(2),
       "size of umbilic patches (as multiple of 1ring size)")
      ("verbose,v", po::value<bool>(&verbose)->default_value(false),
       "verbose output on text file")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cerr << desc << "\n";
      return 1;
    }


    if (vm.count("ridge_order")){
      if ( int_tag == 3 ) tag_order = CGAL::Ridge_order_3;
      if ( int_tag == 4 ) tag_order = CGAL::Ridge_order_4;
      if ( int_tag != 3 && int_tag != 4 )
	{std::cerr << "ridge_order must be CGAL::Ridge_order_3 or CGAL::Ridge_order_4";
	  return 1;}
    }
#else 
    std::cerr << "Command-line options require Boost.ProgramOptions" << std::endl;
    if_name = "data/poly2x^2+y^2-0.062500.off";
    d_fitting = 3;
    d_monge = 3;
    nb_rings = 0;
    nb_points_to_use = 0;
    umb_size = 2;
    verbose = false;
#endif
  }

  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
    }
    catch(...) {
      std::cerr << "Exception of unknown type!\n";
    }

  //modify global variables
  min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

  //prepare output file names
  assert(!if_name.empty());
  of_name = if_name;
  for(unsigned int i=0; i<of_name.size(); i++)
    if (of_name[i] == '/') of_name[i]='_';
  std::ostringstream str_4ogl;
  str_4ogl << "data/"
	   << of_name << "RIDGES"
	   << "-d" << d_fitting
	   << "-m" << d_monge
	   << "-t" << tag_order
	   << "-a" << nb_rings
	   << "-p" << nb_points_to_use
	   << ".4ogl.txt";
  std::cout << str_4ogl.str() << std::endl ;
  std::ofstream out_4ogl(str_4ogl.str().c_str() , std::ios::out);

  //if verbose only...
  std::ostringstream str_verb;
  str_verb << "data/"
	   << of_name << "RIDGES"
	   << "-d" << d_fitting
	   << "-m" << d_monge
	   << "-t" << tag_order
	   << "-a" << nb_rings
	   << "-p" << nb_points_to_use
	   << ".verb.txt";
  std::cout << str_verb.str() << std::endl ;
  std::ofstream out_verb(str_verb.str().c_str() , std::ios::out);

  //load the model from <mesh.off>
  PolyhedralSurf P;
  std::ifstream stream(if_name.c_str());
  stream >> P;
  fprintf(stderr, "loadMesh %d Ves %d Facets\n",
	  (int)num_vertices(P), (int)num_faces(P));
  if(verbose)
    out_verb << "Polysurf with " << num_vertices(P)
	     << " vertices and " << num_faces(P)
	     << " facets. " << std::endl;


vertex_k1_pm = P.add_property_map<vertex_descriptor,FT>("v:k1",0).first;
vertex_k2_pm = P.add_property_map<vertex_descriptor,FT>("v:k2",0).first;
vertex_b0_pm = P.add_property_map<vertex_descriptor,FT>("v:b0",0).first; 
vertex_b3_pm = P.add_property_map<vertex_descriptor,FT>("v:b3",0).first;
vertex_P1_pm = P.add_property_map<vertex_descriptor,FT>("v:P1",0).first; 
vertex_P2_pm = P.add_property_map<vertex_descriptor,FT>("v:P2",0).first;

vertex_d1_pm = P.add_property_map<vertex_descriptor,Vector_3>("v:d1",Vector_3(0,0,0)).first;
vertex_d2_pm = P.add_property_map<vertex_descriptor,Vector_3>("v:d2",Vector_3(0,0,0)).first;

face2normal_pm = P.add_property_map<face_descriptor,Vector_3>("f:n",Vector_3(0,0,0)).first;

  //exit if not enough points in the model
  if (min_nb_points > num_vertices(P))
    {std::cerr << "not enough points in the model" << std::endl;   exit(1);}

  //initialize Polyhedral data : normal of facets
  compute_facets_normals(P,face2normal_pm, Kernel());

  //create a Poly_rings object
  Poly_rings poly_rings(P);

  std::cout << "Compute differential quantities via jet fitting..." << std::endl;
  //initialize the diff quantities property maps
  compute_differential_quantities(P, poly_rings);

  //---------------------------------------------------------------------------
  //Ridges
  //--------------------------------------------------------------------------
  std::cout << "Compute ridges..." << std::endl;
  Ridge_approximation ridge_approximation(P,
					  vertex_k1_pm, vertex_k2_pm,
					  vertex_b0_pm, vertex_b3_pm,
					  vertex_d1_pm, vertex_d2_pm,
					  vertex_P1_pm, vertex_P2_pm );
  std::vector<Ridge_line*> ridge_lines;
  std::back_insert_iterator<std::vector<Ridge_line*> > ii(ridge_lines);

  //Find MAX_RIDGE, MIN_RIDGE, CREST_RIDGES
  //   ridge_approximation.compute_max_ridges(ii, tag_order);
  //   ridge_approximation.compute_min_ridges(ii, tag_order);
  ridge_approximation.compute_crest_ridges(ii, tag_order);

  // or with the global function
  CGAL::compute_max_ridges(P,
			   vertex_k1_pm, vertex_k2_pm,
			   vertex_b0_pm, vertex_b3_pm,
			   vertex_d1_pm, vertex_d2_pm,
			   vertex_P1_pm, vertex_P2_pm,
			   ii, tag_order);

  std::vector<Ridge_line*>::iterator iter_lines = ridge_lines.begin(),
    iter_end = ridge_lines.end();
  //OpenGL output

  typedef boost::property_map<PolyhedralSurf,CGAL::vertex_point_t>::type VPM;
  VPM vpm = get(CGAL::vertex_point,P);

  for (;iter_lines!=iter_end;iter_lines++) (*iter_lines)->dump_4ogl(out_4ogl, vpm);


  for (iter_lines = ridge_lines.begin();iter_lines!=iter_end;iter_lines++){
    //verbose txt output
    if (verbose){
      (*iter_lines)->dump_verbose(out_verb,vpm);
    }
    delete *iter_lines;
    }

  //---------------------------------------------------------------------------
  // UMBILICS
  //--------------------------------------------------------------------------
  std::cout << "Compute umbilics..." << std::endl;
  std::vector<Umbilic*> umbilics;
  std::back_insert_iterator<std::vector<Umbilic*> > umb_it(umbilics);

  //explicit construction of the class
 //  Umbilic_approximation umbilic_approximation(P,
// 					      vertex_k1_pm, vertex_k2_pm,
// 					      vertex_d1_pm, vertex_d2_pm);
//   umbilic_approximation.compute(umb_it, umb_size);
  //or global function call
  CGAL::compute_umbilics(P,
			 vertex_k1_pm, vertex_k2_pm,
			 vertex_d1_pm, vertex_d2_pm,
			 umb_it, umb_size);

  std::vector<Umbilic*>::iterator iter_umb = umbilics.begin(),
    iter_umb_end = umbilics.end();
  // output
  std::cout << "nb of umbilics " << umbilics.size() << std::endl;
  for (;iter_umb!=iter_umb_end;iter_umb++) std::cout << **iter_umb;

  //verbose txt output
  if (verbose) {
    out_verb << "nb of umbilics " << umbilics.size() << std::endl;
  }
  for ( iter_umb = umbilics.begin();iter_umb!=iter_umb_end;iter_umb++){
    if (verbose) {
      out_verb << **iter_umb;
    }
    delete *iter_umb;
  }
  return 0;
}
