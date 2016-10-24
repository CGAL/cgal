#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/ARAP_parameterizer_3.h>
#include <CGAL/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <CGAL/parameterize.h>

#include <CGAL/Simple_cartesian.h>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>            Kernel;
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Point_3                           Point_3;

#define POLY_MESH
#ifdef POLY_MESH
#define MVC_POLY_MESH
#define BARY_POLY_MESH
#define ARAP_POLY_MESH
#endif

#define SURF_MESH
#ifdef SURF_MESH
#define ARAP_SURF_MESH
#endif

//#define SEAM_MESH
#ifdef SEAM_MESH
#define SURF_MESH
#define ARAP_SEMESH
#endif

// #define REDIRECT_OUTPUT

#ifdef POLY_MESH
  typedef CGAL::Polyhedron_3<Kernel> PMesh;

  typedef boost::graph_traits<PMesh>::vertex_descriptor             vertex_descriptor;
  typedef boost::graph_traits<PMesh>::halfedge_descriptor           halfedge_descriptor;
  typedef boost::graph_traits<PMesh>::face_descriptor               face_descriptor;

  typedef CGAL::Parameterizer_traits_3<PMesh>::Error_code           Error_code;
#endif

#ifdef SURF_MESH
  typedef CGAL::Surface_mesh<Point_3> SMesh;

typedef boost::graph_traits<SMesh>::vertex_descriptor     SM_vertex_descriptor;
typedef boost::graph_traits<SMesh>::halfedge_descriptor   SM_halfedge_descriptor;
typedef boost::graph_traits<SMesh>::face_descriptor       SM_face_descriptor;
#endif

#ifdef SEAM_MESH
typedef boost::graph_traits<SMesh>::edge_descriptor               SM_edge_descriptor;

typedef SMesh::Property_map<SM_halfedge_descriptor, Point_2>      UV_pmap;
typedef SMesh::Property_map<SM_edge_descriptor, bool>             Seam_edge_pmap;
typedef SMesh::Property_map<SM_vertex_descriptor, bool>           Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SMesh, Seam_edge_pmap, Seam_vertex_pmap>  TMesh;

typedef boost::graph_traits<TMesh>::vertex_descriptor      TM_vertex_descriptor;
typedef boost::graph_traits<TMesh>::halfedge_descriptor    TM_halfedge_descriptor;
typedef boost::graph_traits<TMesh>::face_descriptor        TM_face_descriptor;
#endif

int main(int argc, char * argv[])
{
  std::cout.precision(17);
  CGAL::set_pretty_mode(std::cout);

#ifdef REDIRECT_OUTPUT
  std::freopen("/home/mrouxell/POLY_MESH.txt", "w", stdout);
#endif

//  std::ifstream in((argc>1)?argv[1]:"../data/tiny_nef.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/nefertiti.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/lion.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/blob.off");
  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/mushroom.off");
//  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/lion-head.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/three_peaks.off");
//  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/three_peaks_dense.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_with_hole.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_dense.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_densified.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/mushroom_hole_2.off");

  if(!in){
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }

  // ***************************************************************************
  // Default case
  // ***************************************************************************

#ifdef MVC_POLY_MESH
  {
    PMesh pm;
    in >> pm;

    halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                          boost::hash<vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                            boost::hash<vertex_descriptor> > > uvpm(uvhm);
    // Go to default (aka Mean values)
    CGAL::parameterize(pm, hd, uvpm);

    std::cout << "Parameterized with Default (Mean Values)!" << std::endl;
  }
#endif

  // ***************************************************************************
  // Barycentric mapping
  // ***************************************************************************

#ifdef BARY_POLY_MESH
  {
    PMesh pm;
    in >> pm;

    halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    // UV map
    CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                          boost::hash<vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                            boost::hash<vertex_descriptor> > > uvpm(uvhm);
    // Indices map
    typedef boost::unordered_map<vertex_descriptor,int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(face(opposite(hd, pm), pm),
                                                       pm,
      boost::make_function_output_iterator(
          CGAL::Parameterization::Vertices<PMesh, Indices>(pm, indices)));

    // Vertex parameterized map
    boost::unordered_set<vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);
    typename CGAL::Barycentric_mapping_parameterizer_3<PMesh> parameterizer;

    parameterizer.parameterize(pm, hd, uvpm, boost::make_assoc_property_map(indices), vpm);

    std::cout << "Parameterized with Barycentric!" << std::endl;
  }
#endif

  // ***************************************************************************
  // ARAP WITH POLY_MESH
  // ***************************************************************************

#ifdef ARAP_POLY_MESH
  {
    PMesh pm;
    in >> pm;

    halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    // UV map
    CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                          boost::hash<vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<vertex_descriptor, Point_2,
                            boost::hash<vertex_descriptor> > > uvpm(uvhm);

    // Indices map
    typedef boost::unordered_map<vertex_descriptor, int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(hd, pm), pm),
      pm,
      boost::make_function_output_iterator(
        CGAL::Parameterization::Vertices<PMesh,
                                         Indices>(pm, indices)));
    boost::associative_property_map<Indices> vipm(indices);

    // Vertex parameterized map
    boost::unordered_set<vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);

    // Parameterizer
    typename CGAL::ARAP_parameterizer_3<PMesh> parameterizer;
    Error_code status =
      parameterizer.parameterize(pm, hd, uvpm, vipm, vpm);

    if(status != CGAL::Parameterizer_traits_3<PMesh>::OK)
      std::cout << "Encountered a problem: " << status << std::endl;
    else
      std::cout << "Parameterized with ARAP!" << std::endl;

    std::ofstream out("ARAP_final_param.txt");
    BOOST_FOREACH(face_descriptor fd, faces(pm)){
      halfedge_descriptor hd = halfedge(fd, pm);
      out << "4 " << uvhm[target(hd, pm)] << " 0 ";
      hd = next(hd, pm);
      BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, pm)){
        out << uvhm[vd] << " 0 ";
      }
      out << std::endl;
    }
  }
#endif

  // ***************************************************************************
  // ARAP WITH SURF_MESH
  // ***************************************************************************

#ifdef REDIRECT_OUTPUT
  std::freopen("/home/mrouxell/SURF_MESH.txt", "w", stdout);
#endif

#ifdef ARAP_SURF_MESH
  {
    SMesh sm;
    in.clear();
    in.seekg(0, std::ios::beg);
    in >> sm;

    SM_halfedge_descriptor bhd =
                        CGAL::Polygon_mesh_processing::longest_border(sm).first;

    CGAL::Unique_hash_map<SM_vertex_descriptor, Point_2,
                          boost::hash<SM_vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<SM_vertex_descriptor,
                            Point_2,
                            boost::hash<SM_vertex_descriptor> > > uv_pm(uvhm);

    // Indices map
    typedef boost::unordered_map<SM_vertex_descriptor, int> Vertex_index_map;
    Vertex_index_map indices;
    CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(bhd, sm), sm),
      sm,
      boost::make_function_output_iterator(
        CGAL::Parameterization::Vertices<SMesh, Vertex_index_map>(sm, indices)));
    boost::associative_property_map<Vertex_index_map> vipm(indices);

    boost::unordered_set<SM_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map< boost::unordered_set<SM_vertex_descriptor> > vpm(vs);

    typename CGAL::ARAP_parameterizer_3<SMesh> parameterizer;
    parameterizer.parameterize(sm, bhd, uv_pm, vipm, vpm);
  }
#endif

  // ***************************************************************************
  // ARAP WITH SEAM_MESH
  // ***************************************************************************

#ifdef ARAP_SEMESH
  {
    SMesh sm;
    in.clear();
    in.seekg(0, std::ios::beg);
    in >> sm;

    Seam_edge_pmap seam_edge_pm =
        sm.add_property_map<SM_edge_descriptor,bool>("e:on_seam", false).first;
    Seam_vertex_pmap seam_vertex_pm =
        sm.add_property_map<SM_vertex_descriptor,bool>("v:on_seam",false).first;

    std::ifstream in((argc>2) ? argv[2] : "../data/lion.selection.txt");
    std::string vertices;
    std::getline(in,vertices);
    std::istringstream iss(vertices);
    int p1, p2;
    bool two_vertices_given = false;
    if(iss >> p1 >> p2){
      two_vertices_given = true;
    }

    int s, t;
    SM_halfedge_descriptor smhd;
    while(in >> s >> t){
      SM_vertex_descriptor svd(s), tvd(t);
      SM_edge_descriptor ed = edge(svd, tvd, sm).first;
      if(! is_border(ed, sm)){
        put(seam_edge_pm, ed, true);
        put(seam_vertex_pm, svd, true);
        put(seam_vertex_pm, tvd, true);
        if(smhd == boost::graph_traits<SMesh>::null_halfedge()){
          smhd = halfedge(edge(svd, tvd, sm).first, sm);
        }
      }
    }

    TMesh mesh(sm, seam_edge_pm, seam_vertex_pm);

    // The 2D points of the uv parametrisation will be written into this map
    // Note that this is a halfedge property map, and that the uv
    // is only stored for the canonical halfedges representing a vertex
    UV_pmap uv_pm = sm.add_property_map<SM_halfedge_descriptor,
                                        Point_2>("h:uv").first;

    TM_halfedge_descriptor bhd(smhd);
    bhd = opposite(bhd, mesh); // a halfedge on the virtual border


    typedef boost::unordered_map<TM_vertex_descriptor, int> Vertex_index_map;
    Vertex_index_map vim;
    boost::associative_property_map<Vertex_index_map> vipm(vim);
    mesh.initialize_vertex_index_map(bhd, vipm);

    boost::unordered_set<TM_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map< boost::unordered_set<TM_vertex_descriptor> > vpm(vs);

    typename CGAL::ARAP_parameterizer_3<TMesh> parameterizer;
    parameterizer.parameterize(mesh, bhd, uv_pm, vipm, vpm);
  }
#endif

  return 0;
}
