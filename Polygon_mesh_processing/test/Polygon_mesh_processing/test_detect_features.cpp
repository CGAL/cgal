#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>               Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor   face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
    const char* filename = (argc > 1) ? argv[1] : "data/cube_quad.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh))
    {
        std::cerr << "Not a valid input file." << std::endl;
        return 1;
    }

    typedef boost::property_map<Mesh,CGAL::face_patch_id_t<int> >::type PatchID;
    typedef boost::property_map<Mesh, CGAL::vertex_incident_patches_t<int> >::type VIP;
    typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIF_map;

    PatchID pid = get(CGAL::face_patch_id_t<int>(), mesh);
    VIP vip     = get(CGAL::vertex_incident_patches_t<int>(), mesh);
    EIF_map eif = get(CGAL::edge_is_feature, mesh);
    std::size_t number_of_patches
      = PMP::sharp_edges_segmentation(mesh, 90, eif, pid);

    std::size_t nb_sharp_edges = 0;
    for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
    {
      if(get(eif, e))
        ++nb_sharp_edges;
    }

    CGAL_assertion(nb_sharp_edges == 12);
    CGAL_assertion(number_of_patches == 6);
    CGAL_USE(number_of_patches);

    number_of_patches
      = PMP::sharp_edges_segmentation(mesh, 90, eif, pid,
                                      PMP::parameters::first_index(1)
                                     .vertex_incident_patches_map(vip));

    CGAL_assertion(number_of_patches == 6);

    PMP::detect_sharp_edges(mesh, 90, eif);
    number_of_patches = PMP::internal::detect_surface_patches(mesh, pid, eif);
    PMP::detect_vertex_incident_patches(mesh, pid, vip, eif);

    nb_sharp_edges = 0;
    for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
    {
      if(get(eif, e))
        ++nb_sharp_edges;
    }


    CGAL_assertion(nb_sharp_edges == 12);
    CGAL_assertion(number_of_patches == 6);

    Mesh::Property_map<face_descriptor,std::pair<int, int> > patch_id_map
            = mesh.add_property_map<face_descriptor,std::pair<int, int> >("f:pid",std::pair<int,int>()).first;
    Mesh::Property_map<vertex_descriptor,std::set<std::pair<int, int> > > vertex_incident_patch_map
            = mesh.add_property_map<vertex_descriptor,std::set<std::pair<int, int> > >("f:vip",std::set<std::pair<int, int> >()).first;
    PMP::detect_sharp_edges(mesh, 90, eif);
    number_of_patches
      = PMP::internal::detect_surface_patches(mesh, patch_id_map, eif,
                                              PMP::parameters::first_index(1));
    PMP::detect_vertex_incident_patches(mesh, patch_id_map, vertex_incident_patch_map, eif);

    nb_sharp_edges =0;
    for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
    {
      if(get(eif, e))
        ++nb_sharp_edges;
    }

    CGAL_assertion(nb_sharp_edges == 12);
    CGAL_assertion(number_of_patches == 6);

    return 0;
}


