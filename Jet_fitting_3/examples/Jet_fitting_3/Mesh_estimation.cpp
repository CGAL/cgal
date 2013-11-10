#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>

#include <fstream>
#include <cassert>

#include <CGAL/property_map.h>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

using namespace std;

#include "PolyhedralSurf.h"
#include "PolyhedralSurf_operations.h"
#include "PolyhedralSurf_rings.h"

//Kernel of the PolyhedralSurf
typedef double                DFT;
typedef CGAL::Simple_cartesian<DFT>  Data_Kernel;
typedef Data_Kernel::Point_3  DPoint;
typedef Data_Kernel::Vector_3 DVector;

//HDS
typedef PolyhedralSurf::Vertex_handle Vertex_handle;
typedef PolyhedralSurf::Vertex Vertex;
typedef PolyhedralSurf::Halfedge_handle Halfedge_handle;
typedef PolyhedralSurf::Halfedge Halfedge;
typedef PolyhedralSurf::Vertex_iterator Vertex_iterator;
typedef PolyhedralSurf::Facet_handle Facet_handle;
typedef PolyhedralSurf::Facet Facet;

struct Hedge_cmp{
  bool operator()(Halfedge_handle a,  Halfedge_handle b) const{
    return &*a < &*b;
  }
};

struct Facet_cmp{
  bool operator()(Facet_handle a, Facet_handle b) const{
    return &*a < &*b;
  }
};

//Vertex property map, with std::map
typedef std::map<Vertex*, int> Vertex2int_map_type;
typedef boost::associative_property_map< Vertex2int_map_type > Vertex_PM_type;
typedef T_PolyhedralSurf_rings<PolyhedralSurf, Vertex_PM_type > Poly_rings;

//Hedge property map, with enriched Halfedge with its length
// typedef HEdge_PM<PolyhedralSurf> Hedge_PM_type;
// typedef T_PolyhedralSurf_hedge_ops<PolyhedralSurf, Hedge_PM_type> Poly_hedge_ops;
//Hedge property map, with std::map
typedef std::map<Halfedge_handle, double, Hedge_cmp> Hedge2double_map_type;
typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
typedef T_PolyhedralSurf_hedge_ops<PolyhedralSurf, Hedge_PM_type> Poly_hedge_ops;

// //Facet property map with enriched Facet with its normal
// typedef Facet_PM<PolyhedralSurf> Facet_PM_type;
// typedef T_PolyhedralSurf_facet_ops<PolyhedralSurf, Facet_PM_type> Poly_facet_ops;
//Facet property map, with std::map
typedef std::map<Facet_handle, Vector_3, Facet_cmp> Facet2normal_map_type;
typedef boost::associative_property_map<Facet2normal_map_type> Facet_PM_type;
typedef T_PolyhedralSurf_facet_ops<PolyhedralSurf, Facet_PM_type> Poly_facet_ops;

typedef double                   LFT;
typedef CGAL::Simple_cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;


// default parameter values and global variables
unsigned int d_fitting = 2;
unsigned int d_monge = 2;
unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
unsigned int nb_points_to_use = 0;//
bool verbose = false;
unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;


//gather points around the vertex v using rings on the
//polyhedralsurf. the collection of points resorts to 3 alternatives:
// 1. the exact number of points to be used
// 2. the exact number of rings to be used
// 3. nothing is specified
void gather_fitting_points(Vertex* v,
			   std::vector<DPoint> &in_points,
			   Vertex_PM_type& vpm)
{
  //container to collect vertices of v on the PolyhedralSurf
  std::vector<Vertex*> gathered;
  //initialize
  in_points.clear();

  //OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
  //enough rings and discard some points of the last collected ring to
  //get the exact "nb_points_to_use"
  if ( nb_points_to_use != 0 ) {
    Poly_rings::collect_enough_rings(v, nb_points_to_use, gathered, vpm);
    if ( gathered.size() > nb_points_to_use ) gathered.resize(nb_points_to_use);
  }
  else { // nb_points_to_use=0, this is the default and the option -p is not considered;
    // then option -a nb_rings is checked. If nb_rings=0, collect
    // enough rings to get the min_nb_points required for the fitting
    // else collect the nb_rings required
    if ( nb_rings == 0 )
      Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);
    else Poly_rings::collect_i_rings(v, nb_rings, gathered, vpm);
  }

  //store the gathered points
  std::vector<Vertex*>::iterator
    itb = gathered.begin(), ite = gathered.end();
  CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
}

///////////////MAIN///////////////////////////////////////////////////////
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
int main(int argc, char *argv[])
#else
int main()
#endif
{
  string if_name_string;
  string if_name; //input file name
  string w_if_name;  //as above, but / replaced by _
  string res4openGL_fname;
  string verbose_fname;
  std::ofstream out_4ogl, out_verbose;

  try {
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message.")
      ("input-file,f", po::value<string>(&if_name_string)->default_value("data/ellipe0.003.off"),
       "name of the input off file")
      ("degree-jet,d", po::value<unsigned int>(&d_fitting)->default_value(2),
       "degree of the jet, 1 <= degre-jet <= 4")
      ("degree-monge,m", po::value<unsigned int>(&d_monge)->default_value(2),
       "degree of the Monge rep, 1 <= degree-monge <= degree-jet")
      ("nb-rings,a", po::value<unsigned int>(&nb_rings)->default_value(0),
       "number of rings to collect neighbors. 0 means collect enough rings to make appro possible a>=1 fixes the nb of rings to be collected")
      ("nb-points,p", po::value<unsigned int>(&nb_points_to_use)->default_value(0),
       "number of neighbors to use.  0 means this option is not considered, this is the default p>=1 fixes the nb of points to be used")
      ("verbose,v", po::value<bool>(&verbose)->default_value(false),
       "verbose output on text file")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }
#else
    std::cerr << "Command-line options require Boost.ProgramOptions" << std::endl;
    if_name_string = "data/ellipe0.003.off";
    d_fitting = 2;
    d_monge = 2;
    nb_rings = 0;
    nb_points_to_use = 0;
    verbose = false;
#endif
  }
  catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }

  //modify global variables which are fct of options:
  min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;
  if (nb_points_to_use < min_nb_points && nb_points_to_use != 0)
    {std::cerr << "the nb of points asked is not enough to perform the fitting" << std::endl; exit(0);}

  //prepare output file names
  //--------------------------
  std::cerr << "if_name_string" << if_name_string  << std::endl;
  if_name = if_name_string;

  w_if_name = if_name;
  for(unsigned int i=0; i<w_if_name.size(); i++)
    if (w_if_name[i] == '/') w_if_name[i]='_';
  cerr << if_name << '\n';
  cerr << w_if_name << '\n';

  res4openGL_fname = w_if_name + ".4ogl.txt";
std::cerr << "res4openGL_fname" << res4openGL_fname  << std::endl;
  out_4ogl.open(res4openGL_fname.c_str(), std::ios::out);
  assert(out_4ogl.good());
  //if verbose only...
  if(verbose){
    verbose_fname  = w_if_name + ".verb.txt";
    out_verbose.open(verbose_fname.c_str(), std::ios::out);
    assert(out_verbose.good());
    CGAL::set_pretty_mode(out_verbose);
  }
  unsigned int nb_vertices_considered = 0;//count vertices for verbose

  //load the model from <mesh.off>
  //------------------------------
  PolyhedralSurf P;
  std::ifstream stream(if_name.c_str());
  stream >> P;
  std::cout << "loadMesh...  "<< "Polysurf with " << P.size_of_vertices()
	    << " vertices and " << P.size_of_facets()
	    << " facets. " << std::endl;

  if(verbose)
    out_verbose << "Polysurf with " << P.size_of_vertices()
                << " vertices and " << P.size_of_facets()
                << " facets. " << std::endl;
  //exit if not enough points in the model
  if (min_nb_points > P.size_of_vertices())    exit(0);

  //create property maps
  //-----------------------------
  //Vertex, using a std::map
  Vertex2int_map_type vertex2props;
  Vertex_PM_type vpm(vertex2props);

  //Hedge, with enriched hedge
  //HEdgePM_type hepm = get_hepm(boost::edge_weight_t(), P);
  //Hedge, using a std::map
  Hedge2double_map_type hedge2props;
  Hedge_PM_type hepm(hedge2props);

  //Facet PM, with enriched Facet
  //FacetPM_type fpm = get_fpm(boost::vertex_attribute_t(), P);
  //Facet PM, with std::map
  Facet2normal_map_type facet2props;
  Facet_PM_type fpm(facet2props);

  //initialize Polyhedral data : length of edges, normal of facets
  Poly_hedge_ops::compute_edges_length(P, hepm);
  Poly_facet_ops::compute_facets_normals(P, fpm);

  //MAIN LOOP: perform calculation for each vertex
  //----------------------------------------------
  std::vector<DPoint> in_points;  //container for data points
  Vertex_iterator vitb, vite;

  //initialize the tag of all vertices to -1
  vitb = P.vertices_begin(); vite = P.vertices_end();
  CGAL_For_all(vitb,vite) put(vpm, &(*vitb), -1);

  vitb = P.vertices_begin(); vite = P.vertices_end();
  for (; vitb != vite; vitb++) {
    //initialize
    Vertex* v = &(*vitb);
    in_points.clear();
    My_Monge_form monge_form;

    //gather points around the vertex using rings
    gather_fitting_points(v, in_points, vpm);

    //skip if the nb of points is to small
    if ( in_points.size() < min_nb_points )
      {std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
	continue;}

    // perform the fitting
    My_Monge_via_jet_fitting monge_fit;
    monge_form = monge_fit(in_points.begin(), in_points.end(),
			   d_fitting, d_monge);
    //switch min-max ppal curv/dir wrt the mesh orientation
    const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
    monge_form.comply_wrt_given_normal(normal_mesh);

    //OpenGL output. Scaling for ppal dir, may be optimized with a
    //global mean edges length computed only once on all edges of P
    DFT scale_ppal_dir = Poly_hedge_ops::compute_mean_edges_length_around_vertex(v, hepm)/2;

    out_4ogl << v->point()  << " ";
    monge_form.dump_4ogl(out_4ogl, scale_ppal_dir);

    //verbose txt output
    if (verbose) {
      std::vector<DPoint>::iterator itbp = in_points.begin(), itep = in_points.end();
      out_verbose << "in_points list : " << std::endl ;
      for (;itbp!=itep;itbp++) out_verbose << *itbp << std::endl ;

      out_verbose << "--- vertex " <<  ++nb_vertices_considered
                  <<	" : " << v->point() << std::endl
                  << "number of points used : " << in_points.size() << std::endl
	;// << monge_form;
    }
  } //all vertices processed

  //cleanup filenames
  //------------------
  out_4ogl.close();
  if(verbose) {
    out_verbose.close();
  }  
  return 0;
}
