#include <CGAL/mesh_segmentation.h>
#include <CGAL/approx_decomposition.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

typedef CGAL::Face_filtered_graph<Polyhedron> Filtered_graph;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template <class ValueType>
struct Facet_with_id_pmap : public boost::put_get_helper<ValueType&, Facet_with_id_pmap<ValueType> >
{
    typedef face_descriptor key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;
    
    Facet_with_id_pmap(std::vector<ValueType>& internal_vector)
    : internal_vector(internal_vector) {}

    reference operator[] (key_type key) const { return internal_vector[key->id()]; }

private:
    std::vector<ValueType>& internal_vector;
};

int main()
{
    std::ifstream input("data/sword.off");
    
    Polyhedron mesh;

    if (!input || !(input >> mesh))
    {
        std::cout << "Failed to read mesh" << std::endl;
        return EXIT_FAILURE;
    }

    if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cout << "Input mesh is invalid" << std::endl;
        return EXIT_FAILURE;
    }  

    // init the polyhedron simplex indices
    CGAL::set_halfedgeds_items_id(mesh);
    
    //for each input vertex compute its distance to the convex hull
    typedef std::map<vertex_descriptor, double> Vertex_double_map;
    Vertex_double_map distances_map;
    boost::associative_property_map<Vertex_double_map> distances_property_map(distances_map);

	std::cout << CGAL::concavity_value<Concurrency_tag>(mesh, distances_property_map) << std::endl;
    
    // create a property-map for sdf values
    std::vector<double> sdf_values(num_faces(mesh));
    Facet_with_id_pmap<double> sdf_property_map(sdf_values);

    // compute sdf values with skeleton
    BOOST_FOREACH(face_descriptor f, faces(mesh))
    {
        double dist = 0;
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, mesh), mesh))
        {
            dist += distances_map[target(hd, mesh)];
        }
        sdf_property_map[f] = dist / 3.;
    }

    // post-process the sdf values
    CGAL::sdf_values_postprocessing(mesh, sdf_property_map);
    
    // create a property-map for segment-ids (it is an adaptor for this case)
    std::vector<std::size_t> segment_ids(num_faces(mesh));
    Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
    
    // segment the mesh using default parameters
    std::size_t clusters_num = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);
    std::cout << "Number of segments: " << clusters_num << std::endl;

    // write clusters to disk
    for (std::size_t i = 0; i < clusters_num; ++i)
    {
        Filtered_graph filtered_mesh(mesh, i, segment_property_map);
        Polyhedron cluster;
        CGAL::copy_face_graph(filtered_mesh, cluster);

        std::ofstream out("cluster_" + std::to_string(i) + ".off");
        out << cluster;
        out.close();
	}

    return EXIT_SUCCESS;
}
