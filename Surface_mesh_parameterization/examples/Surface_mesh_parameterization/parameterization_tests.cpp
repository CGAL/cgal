#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

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

#define POLY_MESH // Polyhedron mesh
#ifdef POLY_MESH
#define MVC_POLY_MESH
#define BARY_POLY_MESH
#define ARAP_POLY_MESH
#endif

#define SURF_MESH // Surface_mesh
#ifdef SURF_MESH
#define ARAP_SURF_MESH
#endif

#define PM_SEAM_MESH // Polyhedron-based seam mesh
#ifdef PM_SEAM_MESH
#define POLY_MESH
#define ARAP_PM_SEAM_MESH
#endif

#define SM_SEAM_MESH // Surface_mesh-based seam mesh
#ifdef SM_SEAM_MESH
#define SURF_MESH
#define ARAP_SM_SEAM_MESH
#endif

// #define REDIRECT_OUTPUT

#ifdef POLY_MESH
typedef CGAL::Polyhedron_3<Kernel>                                PMesh;

typedef boost::graph_traits<PMesh>::vertex_descriptor             PM_vertex_descriptor;
typedef boost::graph_traits<PMesh>::halfedge_descriptor           PM_halfedge_descriptor;

typedef CGAL::Unique_hash_map<PM_halfedge_descriptor, Point_2>    PM_UV_hmap;
typedef boost::associative_property_map<PM_UV_hmap>               PM_UV_pmap;

typedef CGAL::Parameterizer_traits_3<PMesh>::Error_code           Error_code;
#endif

#ifdef SURF_MESH
typedef CGAL::Surface_mesh<Point_3>                               SMesh;

typedef boost::graph_traits<SMesh>::vertex_descriptor             SM_vertex_descriptor;
typedef boost::graph_traits<SMesh>::halfedge_descriptor           SM_halfedge_descriptor;

typedef SMesh::Property_map<SM_halfedge_descriptor, Point_2>      SM_UV_pmap;
#endif

#ifdef SM_SEAM_MESH
typedef boost::graph_traits<SMesh>::edge_descriptor               SM_edge_descriptor;

typedef SMesh::Property_map<SM_edge_descriptor, bool>             SM_seam_edge_pmap;
typedef SMesh::Property_map<SM_vertex_descriptor, bool>           SM_seam_vertex_pmap;

typedef CGAL::Seam_mesh<SMesh, SM_seam_edge_pmap, SM_seam_vertex_pmap>
                                                                  SM_Seam_mesh;

typedef boost::graph_traits<SM_Seam_mesh>::vertex_descriptor      SM_SE_vertex_descriptor;
typedef boost::graph_traits<SM_Seam_mesh>::halfedge_descriptor    SM_SE_halfedge_descriptor;
#endif

#ifdef PM_SEAM_MESH
typedef boost::graph_traits<PMesh>::edge_descriptor               PM_edge_descriptor;

typedef CGAL::Unique_hash_map<PM_edge_descriptor, bool>           PM_seam_edge_hmap;
typedef boost::associative_property_map<PM_seam_edge_hmap>        PM_seam_edge_pmap;
typedef CGAL::Unique_hash_map<PM_vertex_descriptor, bool>         PM_seam_vertex_hmap;
typedef boost::associative_property_map<PM_seam_vertex_hmap>      PM_seam_vertex_pmap;

typedef CGAL::Seam_mesh<PMesh, PM_seam_edge_pmap, PM_seam_vertex_pmap>
                                                                  PM_Seam_mesh;

typedef boost::graph_traits<PM_Seam_mesh>::vertex_descriptor      PM_SE_vertex_descriptor;
typedef boost::graph_traits<PM_Seam_mesh>::halfedge_descriptor    PM_SE_halfedge_descriptor;
#endif

int main(int argc, char * argv[])
{
  std::cout.precision(17);
  CGAL::set_pretty_mode(std::cout);

#ifdef REDIRECT_OUTPUT
  std::freopen("/home/mrouxell/POLY_MESH.txt", "w", stdout);
#endif

  std::ifstream in((argc>1)?argv[1]:"../data/tiny_nef.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/nefertiti.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/lion.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/blob.off");
//  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/mushroom.off");
//  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/lion-head.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/three_peaks.off");
//  std::ifstream in((argc>1)?argv[1]:"/home/mrouxell/Data/OFF/three_peaks_dense.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_with_hole.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_dense.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/cow_densified.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/mushroom_big_hole.off");
//  std::ifstream in((argc>1)?argv[1]:"../data/mushroom_hole_1.off");

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

    PM_halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                          boost::hash<PM_vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                            boost::hash<PM_vertex_descriptor> > > uvpm(uvhm);
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

    PM_halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    // UV map
    CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                          boost::hash<PM_vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                            boost::hash<PM_vertex_descriptor> > > uvpm(uvhm);
    // Indices map
    typedef boost::unordered_map<PM_vertex_descriptor, int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(face(opposite(hd, pm), pm),
                                                       pm,
      boost::make_function_output_iterator(
          CGAL::Parameterization::Vertices<PMesh, Indices>(pm, indices)));

    // Vertex parameterized map
    boost::unordered_set<PM_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<PM_vertex_descriptor> > vpm(vs);
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

    PM_halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(pm).first;

    // UV map
    CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                          boost::hash<PM_vertex_descriptor> > uvhm;
    boost::associative_property_map<
      CGAL::Unique_hash_map<PM_vertex_descriptor, Point_2,
                            boost::hash<PM_vertex_descriptor> > > uvpm(uvhm);

    // Indices map
    typedef boost::unordered_map<PM_vertex_descriptor, int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(hd, pm), pm),
      pm,
      boost::make_function_output_iterator(
        CGAL::Parameterization::Vertices<PMesh,
                                         Indices>(pm, indices)));
    boost::associative_property_map<Indices> vipm(indices);

    // Vertex parameterized map
    boost::unordered_set<PM_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<PM_vertex_descriptor> > vpm(vs);

    // Parameterizer
    typename CGAL::ARAP_parameterizer_3<PMesh> parameterizer;
    Error_code status =
      parameterizer.parameterize(pm, hd, uvpm, vipm, vpm);

    if(status != CGAL::Parameterizer_traits_3<PMesh>::OK)
      std::cout << "Encountered a problem: " << status << std::endl;
    else
      std::cout << "Parameterized with ARAP!" << std::endl;
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
  // ARAP WITH SEAM_MESH (SM)
  // ***************************************************************************

#ifdef ARAP_SM_SEAM_MESH
  {
    SMesh sm;
    in.clear();
    in.seekg(0, std::ios::beg);
    in >> sm;

    SM_seam_edge_pmap seam_edge_pm =
        sm.add_property_map<SM_edge_descriptor,bool>("e:on_seam", false).first;
    SM_seam_vertex_pmap seam_vertex_pm =
        sm.add_property_map<SM_vertex_descriptor,bool>("v:on_seam",false).first;

    SM_Seam_mesh mesh(sm, seam_edge_pm, seam_vertex_pm);
    SM_halfedge_descriptor smhd = mesh.add_seams("../data/lion.selection.txt");

    // The 2D points of the uv parametrisation will be written into this map
    // Note that this is a halfedge property map, and that the uv
    // is only stored for the canonical halfedges representing a vertex
    SM_UV_pmap uv_pm = sm.add_property_map<SM_halfedge_descriptor,
                                        Point_2>("h:uv").first;

    SM_SE_halfedge_descriptor bhd(smhd);
    bhd = opposite(bhd, mesh); // a halfedge on the virtual border

    // Indices
    typedef boost::unordered_map<SM_SE_vertex_descriptor, int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(
              face(opposite(bhd, mesh), mesh),
              mesh,
              boost::make_function_output_iterator(
                CGAL::Parameterization::Vertices<SM_Seam_mesh, Indices>(mesh, indices)));
    boost::associative_property_map<Indices> vipm(indices);

    // Parameterized
    boost::unordered_set<SM_SE_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<SM_SE_vertex_descriptor> > vpm(vs);

    typename CGAL::ARAP_parameterizer_3<SM_Seam_mesh> parameterizer;
    parameterizer.parameterize(mesh, bhd, uv_pm, vipm, vpm);

  }
#endif

  // ***************************************************************************
  // ARAP WITH SEAM_MESH (POLY)
  // ***************************************************************************

#ifdef ARAP_PM_SEAM_MESH
  {
    PMesh pm;
    in.clear();
    in.seekg(0, std::ios::beg);
    in >> pm;

    PM_seam_edge_hmap seam_edge_hm(false);
    PM_seam_edge_pmap seam_edge_pm(seam_edge_hm);
    PM_seam_vertex_hmap seam_vertex_hm(false);
    PM_seam_vertex_pmap seam_vertex_pm(seam_vertex_hm);

    PM_Seam_mesh mesh(pm, seam_edge_pm, seam_vertex_pm);
    PM_halfedge_descriptor pmhd = mesh.add_seams("../data/lion.selection.txt");

    // The 2D points of the uv parametrisation will be written into this map
    // Note that this is a halfedge property map, and that the uv
    // is only stored for the canonical halfedges representing a vertex
    PM_UV_hmap uv_hm;
    PM_UV_pmap uv_pm(uv_hm);

    PM_SE_halfedge_descriptor bhd(pmhd);
    bhd = opposite(bhd, mesh); // a halfedge on the virtual border

    // Indices
    typedef boost::unordered_map<PM_SE_vertex_descriptor, int> Indices;
    Indices indices;
    CGAL::Polygon_mesh_processing::connected_component(
              face(opposite(bhd, mesh), mesh),
              mesh,
              boost::make_function_output_iterator(
                CGAL::Parameterization::Vertices<PM_Seam_mesh, Indices>(mesh, indices)));
    boost::associative_property_map<Indices> vipm(indices);

    // Parameterized
    boost::unordered_set<PM_SE_vertex_descriptor> vs;
    CGAL::internal::Bool_property_map<boost::unordered_set<PM_SE_vertex_descriptor> > vpm(vs);

    typename CGAL::ARAP_parameterizer_3<PM_Seam_mesh> parameterizer;
    parameterizer.parameterize(mesh, bhd, uv_pm, vipm, vpm);
  }
#endif

  return 0;
}
