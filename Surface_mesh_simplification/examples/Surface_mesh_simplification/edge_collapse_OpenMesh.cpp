#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Surface_mesh;

typedef CGAL::Simple_cartesian<double> Kern;
typedef Kern::Point_3 Point_3;
typedef Kern::Vector_3 Vector_3;

struct OM_Squared_distance_3 {
  double operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    return CGAL::squared_distance(Point_3(p[0],p[1],p[2]),Point_3(q[0],q[1],q[2]));
  }
};

struct OM_Midpoint_3 {
  Surface_mesh::Point operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    Point_3 m = CGAL::midpoint(Point_3(p[0],p[1],p[2]),Point_3(q[0],q[1],q[2]));
    std::cerr << m << std::endl;
    return Surface_mesh::Point(m.x(), m.y(), m.z());
  }
};

struct OM_Construct_cross_product_vector_3 {
  Surface_mesh::Point operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    Vector_3 v = CGAL::cross_product(Point_3(p[0],p[1],p[2]) - CGAL::ORIGIN,Point_3(q[0],q[1],q[2])-CGAL::ORIGIN);
    return Surface_mesh::Point(v.x(), v.y(), v.z());
  }
};

struct OM_Compute_scalar_product_3 {
  double operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    return (Point_3(p[0],p[1],p[2]) - CGAL::ORIGIN) * (Point_3(q[0],q[1],q[2])-CGAL::ORIGIN);
  }
};

struct OM_Construct_vector_3 {
  Surface_mesh::Point operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    return Surface_mesh::Point(p[0]-q[0],p[1]-q[1],p[2]-q[2]);
  }
};

struct OM_Equal_3 {
  bool operator()(const Surface_mesh::Point& p, const Surface_mesh::Point& q) const
  {
    return p[0]==q[0] && p[1]==q[1] && p[2]==q[2];
  }
};

class OM_kernel {
  
public:
  typedef Surface_mesh::Point Point_3;
  typedef Surface_mesh::Point Vector_3;
  typedef double  FT;

  typedef OM_Squared_distance_3 Squared_distance_3;
  typedef OM_Midpoint_3 Midpoint_3;
  typedef OM_Construct_vector_3 Construct_vector_3;
  typedef OM_Construct_cross_product_vector_3 Construct_cross_product_vector_3;
  typedef OM_Compute_scalar_product_3 Compute_scalar_product_3;
  typedef OM_Equal_3 Equal_3;

  Squared_distance_3
  compute_squared_distance_3_object() const
  {
    return Squared_distance_3();
  }

  Midpoint_3
  construct_midpoint_3_object() const
  {
    return Midpoint_3();
  }

  Construct_vector_3
  construct_vector_3_object() const
  {
    return Construct_vector_3();
  }

  Construct_cross_product_vector_3
  construct_cross_product_vector_3_object() const
  {
    return Construct_cross_product_vector_3();
  }


  Compute_scalar_product_3
  compute_scalar_product_3_object() const
  {
    return Compute_scalar_product_3();
  }

Equal_3
  equal_3_object() const
  {
    return Equal_3();
  }
};




namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface_mesh surface_mesh;
  
 OpenMesh::IO::read_mesh(surface_mesh, argv[1]);

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface_mesh> stop(100);
     
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.

  int r = SMS::edge_collapse
            (surface_mesh
            ,stop
             , OM_kernel()
                          , CGAL::get_cost (SMS::Edge_length_cost <Surface_mesh, OM_kernel>())
                          .get_placement(SMS::Midpoint_placement<Surface_mesh,OM_kernel>())
             );
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << num_edges(surface_mesh) << " final edges.\n" ;
        
  //  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;
  
  return 0 ;      
}

// EOF //
