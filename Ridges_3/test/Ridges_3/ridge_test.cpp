#include <cassert>
#include <fstream>
#include <vector>
//This Is an enriched Polyhedron with facets' normal
#include "PolyhedralSurf.h"
#include "PolyhedralSurf_rings.h"
#include "compute_normals.h"
#include <CGAL/Ridges.h>
#include <CGAL/Umbilics.h>
#include <CGAL/Monge_via_jet_fitting.h>

// Functions declared in PolyhedralSurf.h
// They were previously defined in a separate file PolyhedralSurf.cpp,
// but I prefer to avoid custom CMakeLists.txt files in the testsuite.
// -- Laurent Rineau, 2008/11/10

typedef PolyhedralSurf::Traits          Kernel;
typedef Kernel::FT                      FT;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;

typedef boost::graph_traits<PolyhedralSurf>::vertex_descriptor  vertex_descriptor;
typedef boost::graph_traits<PolyhedralSurf>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<PolyhedralSurf>::face_descriptor face_descriptor;

typedef T_PolyhedralSurf_rings<PolyhedralSurf> Poly_rings;
typedef CGAL::Monge_via_jet_fitting<Kernel>    Monge_via_jet_fitting;
typedef Monge_via_jet_fitting::Monge_form      Monge_form;


typedef std::map<vertex_descriptor, FT> VertexFT_map;
typedef boost::associative_property_map< VertexFT_map > VertexFT_property_map;

typedef std::map<vertex_descriptor, Vector_3> VertexVector_map;
typedef boost::associative_property_map< VertexVector_map > VertexVector_property_map;

typedef std::map<face_descriptor, Vector_3> Face2Vector_map;
typedef boost::associative_property_map< Face2Vector_map > Face2Vector_property_map;
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

//create property maps, to be moved in main?
VertexFT_map vertex_k1_map, vertex_k2_map,
  vertex_b0_map, vertex_b3_map,
  vertex_P1_map, vertex_P2_map;
VertexVector_map vertex_d1_map, vertex_d2_map;
Face2Vector_map face2normal_map;

VertexFT_property_map vertex_k1_pm(vertex_k1_map), vertex_k2_pm(vertex_k2_map),
  vertex_b0_pm(vertex_b0_map), vertex_b3_pm(vertex_b3_map),
  vertex_P1_pm(vertex_P1_map), vertex_P2_pm(vertex_P2_map);
VertexVector_property_map vertex_d1_pm(vertex_d1_map), vertex_d2_pm(vertex_d2_map);
Face2Vector_property_map  face2normal_pm(face2normal_map);

  // default fct parameter values and global variables
  unsigned int d_fitting = 4;
  unsigned int d_monge = 4;
  unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
  unsigned int nb_points_to_use = 0;//
  CGAL::Ridge_order tag_order = CGAL::Ridge_order_3;
  double umb_size = 1;
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
    vertex_d1_map[v] = monge_form.maximal_principal_direction();
    vertex_d2_map[v] = monge_form.minimal_principal_direction();
    vertex_k1_map[v] = monge_form.coefficients()[0];
    vertex_k2_map[v] = monge_form.coefficients()[1];
    vertex_b0_map[v] = monge_form.coefficients()[2];
    vertex_b3_map[v] = monge_form.coefficients()[5];
    if ( d_monge >= 4) {
      //= 3*b1^2+(k1-k2)(c0-3k1^3)
      vertex_P1_map[v] =
        3*monge_form.coefficients()[3]*monge_form.coefficients()[3]
        +(monge_form.coefficients()[0]-monge_form.coefficients()[1])
        *(monge_form.coefficients()[6]
          -3*monge_form.coefficients()[0]*monge_form.coefficients()[0]
          *monge_form.coefficients()[0]);
      //= 3*b2^2+(k2-k1)(c4-3k2^3)
      vertex_P2_map[v] =
        3*monge_form.coefficients()[4]*monge_form.coefficients()[4]
        +(-monge_form.coefficients()[0]+monge_form.coefficients()[1])
        *(monge_form.coefficients()[10]
          -3*monge_form.coefficients()[1]*monge_form.coefficients()[1]
          *monge_form.coefficients()[1]);
    }
  } //END FOR LOOP
}

int main()
{
  //load the model from <mesh.off>
  PolyhedralSurf P;
  std::ifstream stream("data/ellipsoid.off");
  stream >> P;
  fprintf(stderr, "loadMesh %d Ves %d Facets\n",
          (int)P.size_of_vertices(), (int)P.size_of_facets());

  //exit if not enough points in the model
  if (min_nb_points > P.size_of_vertices())
    {std::cerr << "not enough points in the model" << std::endl;   return 1;}

  //initialize Polyhedral data : normal of facets
  compute_facets_normals(P,face2normal_pm, Kernel());

  //create a Poly_rings object
  Poly_rings poly_rings(P);

  std::cout << "Compute differential quantities via jet fitting..." << std::endl;
  //initialize the diff quantities property maps
  compute_differential_quantities(P, poly_rings);

  std::cout << "Compute ridges with tag_3" << std::endl;
  //---------------------------------------------------------------------------
  //Ridges
  //--------------------------------------------------------------------------
   Ridge_approximation ridge_approximation(P,
                                          vertex_k1_pm, vertex_k2_pm,
                                          vertex_b0_pm, vertex_b3_pm,
                                          vertex_d1_pm, vertex_d2_pm,
                                          vertex_P1_pm, vertex_P2_pm );
  std::vector<Ridge_line*> ridge_lines;
  std::back_insert_iterator<std::vector<Ridge_line*> > ii(ridge_lines);

  //Find MAX_RIDGE, RED_RIDGE, CREST or all ridges
  ridge_approximation.compute_max_ridges(ii, tag_order);
  ridge_approximation.compute_min_ridges(ii, tag_order);
  ridge_approximation.compute_crest_ridges(ii, tag_order);

  std::cout << "Compute ridges with tag_4" << std::endl;
  tag_order =  CGAL::Ridge_order_4;
   //Find MAX_RIDGE, RED_RIDGE, CREST or all ridges

  ridge_approximation.compute_max_ridges(ii, tag_order);
  ridge_approximation.compute_min_ridges(ii, tag_order);
  ridge_approximation.compute_crest_ridges(ii, tag_order);
  std::vector<Ridge_line*>::iterator iter_lines = ridge_lines.begin(),
    iter_end = ridge_lines.end();
  for (;iter_lines!=iter_end;iter_lines++){
    delete *iter_lines;
  }
  //---------------------------------------------------------------------------
  // UMBILICS
  //--------------------------------------------------------------------------
  Umbilic_approximation umbilic_approximation(P,
                                              vertex_k1_pm, vertex_k2_pm,
                                              vertex_d1_pm, vertex_d2_pm);
  std::vector<Umbilic*> umbilics;
  std::back_insert_iterator<std::vector<Umbilic*> > umb_it(umbilics);
  std::cout << "compute umbilics u=1" << std::endl;
  umbilic_approximation.compute(umb_it, umb_size);
  umb_size=2;
  std::cout << "compute umbilics u=2" << std::endl;
  umb_size=3.5;
  std::cout << "compute umbilics u=3.5" << std::endl;

  assert(umbilics.size() == 4);
  std::vector<Umbilic*>::iterator iter_umb = umbilics.begin(),
   iter_umb_end = umbilics.end();
  // output
  std::cout << "nb of umbilics " << umbilics.size() << std::endl;
  for (;iter_umb!=iter_umb_end;iter_umb++){
    std::cout << **iter_umb;
    delete *iter_umb;
  }
  std::cout << "success\n";
  return 0;
}

