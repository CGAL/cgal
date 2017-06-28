#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_iterator edge_iterator;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

class Constrained_edge_map
{
public:
  typedef boost::read_write_property_map_tag    category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(Surface_mesh& sm)
    : sm_(sm)
  {
    sm_.add_property(constraint);
  }

  inline friend reference get(const Constrained_edge_map& em, key_type e)
  {
    bool b = em.sm_.property(em.constraint,em.sm_.edge_handle(e.idx())); 
    return b;
  }
  
  inline friend void put(const Constrained_edge_map& em, key_type e, value_type b)
  {
    em.sm_.property(em.constraint,em.sm_.edge_handle(e.idx())) = b;
  }

private:
  Surface_mesh& sm_;
  OpenMesh::EPropHandleT<bool> constraint;
};


class Face_partition_map
{
public:
  typedef boost::read_write_property_map_tag    category;
  typedef int                                  value_type;
  typedef int                                  reference;
  typedef face_descriptor                       key_type;

  Face_partition_map(Surface_mesh& sm)
    : sm_(sm)
  {
    sm_.add_property(constraint);
  }

  inline friend reference get(const Face_partition_map& em, key_type e)
  {
    int i = em.sm_.property(em.constraint,em.sm_.face_handle(e.idx())); 
    return i;
  }
  
  inline friend void put(const Face_partition_map& em, key_type f, value_type i)
  {
    em.sm_.property(em.constraint,em.sm_.face_handle(f.idx())) = i;
  }

private:
  Surface_mesh& sm_;
  OpenMesh::FPropHandleT<int> constraint;
};


namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface_mesh surface_mesh;
  Constrained_edge_map constraints_map(surface_mesh);
  Face_partition_map face_partition_map(surface_mesh);
  if (argc==2)
    OpenMesh::IO::read_mesh(surface_mesh, argv[1]);
  else
    OpenMesh::IO::read_mesh(surface_mesh, "cube.off");

  if (!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // For the pupose of the example we mark 10 edges as constrained edges
  edge_iterator b,e;
  int count=0;
  for(boost::tie(b,e) = edges(surface_mesh); b!= e; ++b){
      put(constraints_map,*b,(count++ <100));
  }
  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface_mesh> stop(0);
     
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.

  int r = SMS::parallel_edge_collapse
            (surface_mesh
             ,stop
             , face_partition_map
             ,1
             ,CGAL::parameters::halfedge_index_map  (get(CGAL::halfedge_index  ,surface_mesh)) 
                               .vertex_point_map(get(boost::vertex_point, surface_mesh))
                               .edge_is_constrained_map(constraints_map) 
             );
  
  surface_mesh.garbage_collection();
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << num_edges(surface_mesh) << " final edges.\n" ;
        
   OpenMesh::IO::write_mesh(surface_mesh, "out.off");
  
  return EXIT_SUCCESS;
}
