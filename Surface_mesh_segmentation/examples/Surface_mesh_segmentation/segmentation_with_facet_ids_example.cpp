#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Facet_with_id_pmap<ValueType> >
{
    typedef Polyhedron::Facet_const_handle key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    Facet_with_id_pmap(
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
    if ( !input || !(input >> mesh) || mesh.empty() ) {
      std::cerr << "Not a valid off file." << std::endl;
      return EXIT_FAILURE;
    }

    // assign id field for each facet
    std::size_t facet_id = 0;
    for(Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
        facet_it->id() = facet_id;
    }

    // create a property-map for SDF values
    std::vector<double> sdf_values(mesh.size_of_facets());
    Facet_with_id_pmap<double> sdf_property_map(sdf_values);

    CGAL::sdf_values(mesh, sdf_property_map);

    // access SDF values (with constant-complexity)
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
        std::cout << sdf_property_map[facet_it] << " ";
    }
    std::cout << std::endl;

    // create a property-map for segment-ids
    std::vector<std::size_t> segment_ids(mesh.size_of_facets());
    Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

    CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

    // access segment-ids (with constant-complexity)
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
        std::cout << segment_property_map[facet_it] << " ";
    }
    std::cout << std::endl;
}
