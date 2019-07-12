#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
  Mesh tm1, tm2;
  // TODO Read real meshes and try method with it.
  //      std::ofstream(argv[1]) >> tm1;
  //      http://www.pmp-book.org/download/meshes/iphi_fullres.zip
  //      https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga028a80dc84395650f67714fa7618ec53
/*
  CGAL::make_tetrahedron(Point_3(.0,.0,.0),
                         Point_3(2,.0,.0),
                         Point_3(1,1,1),
                         Point_3(1,.0,2),
                         tm1);
  tm2=tm1;
  CGAL::Polygon_mesh_processing::isotropic_remeshing(tm2.faces(),.05, tm2);
*/
  tm1 = Mesh();
  tm1.add_vertex( Point_3(-1.,1.,1.) );
  tm1.add_vertex( Point_3(0.,-1.,1.) );
  tm1.add_vertex( Point_3(1.,1.,1.) );
  tm1.add_face( tm1.vertices() );

  std::cout << "TM1 is valid: " << (tm1.is_valid() ? "true" : "false") << std::endl;

  tm2 = Mesh();
  Vertex_index w0 = tm2.add_vertex( Point_3(-1.,1.,0.) );
  Vertex_index w1 = tm2.add_vertex( Point_3(0.,-1.,0.) );
  Vertex_index w2 = tm2.add_vertex( Point_3(1.,1.,0.) );
  Vertex_index w3 = tm2.add_vertex( Point_3(0.,1.,-1.) );
  Vertex_index w4 = tm2.add_vertex( Point_3(-0.5,0.,-1.) );
  Vertex_index w5 = tm2.add_vertex( Point_3(0.5,0.,-1.) );
  tm2.add_face( w0, w3, w4 );
  tm2.add_face( w1, w4, w5 );
  tm2.add_face( w2, w5, w3 );

  std::cout << "TM2 is valid: " << (tm2.is_valid() ? "true" : "false") << std::endl;

  std::cout << "Approximated Hausdorff distance: "
            << CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance
                  <TAG>(tm1, tm2, 0.001)
            << std::endl;
}
