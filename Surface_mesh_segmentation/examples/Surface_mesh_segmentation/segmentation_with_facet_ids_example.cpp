#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

// Property map associating a face with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Face_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Face_with_id_pmap<ValueType> >
{
    typedef face_descriptor key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    Face_with_id_pmap(
      std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }

    reference operator[](key_type key) const
    { return internal_vector[key->id()]; }
private:
    std::vector<ValueType>& internal_vector;
};

int main()
{
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input("data/cactus.off");
    if ( !input || !(input >> mesh) || mesh.empty() || ( !CGAL::is_triangle_mesh(mesh)) ) {
      std::cerr << "Input is not a triangle mesh" << std::endl;
      return EXIT_FAILURE;
    }

    // assign id field for each face
    std::size_t face_id = 0;
    for(face_descriptor f : faces( mesh) ) {
        f->id() = face_id++;
    }

    // create a property-map for SDF values
    std::vector<double> sdf_values(num_faces(mesh));
    Face_with_id_pmap<double> sdf_property_map(sdf_values);

    CGAL::sdf_values(mesh, sdf_property_map);

    // access SDF values (with constant-complexity)
    for(face_descriptor f : faces(mesh)) {
        std::cout << sdf_property_map[f] << " ";
    }
    std::cout << std::endl;

    // create a property-map for segment-ids
    std::vector<std::size_t> segment_ids(num_faces(mesh));
    Face_with_id_pmap<std::size_t> segment_property_map(segment_ids);

    CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

    // access segment-ids (with constant-complexity)
    for(face_descriptor f : faces(mesh)) {
        std::cout << segment_property_map[f] << " ";
    }
    std::cout << std::endl;
    return EXIT_SUCCESS;
}
