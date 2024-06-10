#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

typedef CGAL::Surface_mesh<CGAL::Epeck::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<CGAL::Epeck> Polyhedron_3;
typedef CGAL::Epeck::Point_3 Point_3;

int main()
{
  Point_3 p(0.3, 0.3454, 0), q(0.6694, 0,0), r(0,0,0);
  Surface_mesh sm;
  Polyhedron_3 po;
  CGAL::make_triangle(p,q,r,sm);
  CGAL::make_triangle(p,q,r,po);

  // from <CGAL/boost/graph/IO/polygon_mesh_io.h>
  // https://doc.cgal.org/latest/BGL/group__PkgBGLIoFuncsPLY.html#ga959dcd88ca979d3b6b0806d883a0247f
  CGAL::IO::write_polygon_mesh("write_polygon_mesh_sm.ply", sm,
                               CGAL::parameters::stream_precision(17).use_binary_mode(true));

  Surface_mesh osm;
  if(!CGAL::IO::read_polygon_mesh("write_polygon_mesh_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'write_polygon_mesh_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }

  // ERROR this produces an invalid plyfile
  CGAL::IO::write_polygon_mesh("write_polygon_mesh_po.ply", po,
                               CGAL::parameters::stream_precision(17).use_binary_mode(true));
  if(!CGAL::IO::read_polygon_mesh("write_polygon_mesh_po.ply", osm))
  {
    std::cerr << "Error: failed to read 'write_polygon_mesh_po.ply'" << std::endl;
    return EXIT_FAILURE;
  }

  //  OK
  // from #include <CGAL/boost/graph/IO/PLY.h>
  // https://doc.cgal.org/latest/BGL/group__PkgBGLIoFuncsPLY.html#ga959dcd88ca979d3b6b0806d883a0247f
  CGAL::IO::write_PLY("generic_write_PLY_sm.ply", sm, "generic write_PLY(Surface_mesh)",
                      CGAL::parameters::stream_precision(17).use_binary_mode(true));
  if(!CGAL::IO::read_polygon_mesh("generic_write_PLY_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'generic_write_PLY_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }

 // ERROR this produces an invalid plyfile
  CGAL::IO::write_PLY("generic_write_PLY_po.ply", po, "generic write_PLY(Polyhedron)",
                      CGAL::parameters::stream_precision(17).use_binary_mode(true));
  if(!CGAL::IO::read_polygon_mesh("generic_write_PLY_po.ply", osm))
  {
    std::cerr << "Error: failed to read 'generic_write_PLY_po.ply'" << std::endl;
    return EXIT_FAILURE;
  }

  // OK
  // from #include <CGAL/Surface_mesh/IO/PLY.h>
  // https://doc.cgal.org/latest/Surface_mesh/group__PkgSurfaceMeshIOFuncPLY.html#ga50f0e9f2b293855d2c7f1a62939cbe8d
  std::ofstream out("overloaded_write_PLY_sm.ply", std::ios::binary);
  CGAL::IO::write_PLY(out, sm, "overloaded_write_PLY(Surface_mesh)");
  if(!CGAL::IO::read_polygon_mesh("overloaded_write_PLY_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'overloaded_write_PLY_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }

  return 0;
}
