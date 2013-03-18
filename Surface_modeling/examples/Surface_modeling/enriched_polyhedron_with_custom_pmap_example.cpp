#include <CGAL/Deform_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Eigen_solver_traits.h>

#include <fstream>
#include <boost/property_map/property_map.hpp>

#include <Eigen/SuperLUSupport>

#include <CGAL/Polyhedron_items_with_id_3.h>

// Property map using extra field in KeyType
template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;
        
    reference operator[](key_type key) const { return key->id(); }
};

typedef CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> > DefaultSolver;
  
typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron; //enriched polyhedron

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor  	  edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Deform_mesh<Polyhedron, DefaultSolver, Vertex_index_map, Edge_index_map> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream("models/plane.off") >> mesh;

  Deform_mesh deform_mesh(mesh, Vertex_index_map(), Edge_index_map()); 

  // Deform mesh as you wish
}

