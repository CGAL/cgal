#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/generators.h>

typedef CGAL::Surface_mesh<CGAL::Epeck::Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Polyhedron_3<CGAL::Epeck> Polyhedron_3;
typedef boost::graph_traits<Polyhedron_3>::vertex_descriptor pvertex_descriptor;
typedef CGAL::Epeck::Point_3 Point_3;
typedef CGAL::Epeck::Vector_3 Vector_3;
typedef CGAL::Simple_cartesian<double> SC;
typedef Surface_mesh::template Property_map<vertex_descriptor, Vector_3> VNMap;
typedef CGAL::Cartesian_converter_property_map< SC::Vector_3, VNMap> ScSmVertexNormalMap;

typedef std::map<pvertex_descriptor,Vector_3> PolyhedronVN;
typedef boost::associative_property_map<PolyhedronVN> PolyhedronVNMap;
typedef CGAL::Cartesian_converter_property_map< SC::Vector_3, PolyhedronVNMap> ScPolyhedronVertexNormalMap;

typedef boost::property_map<Surface_mesh, CGAL::vertex_point_t >::type SmVertexPointMap;
typedef CGAL::Cartesian_converter_property_map< SC::Point_3, SmVertexPointMap> ScSmVertexPointMap;

typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t >::type PolyhedronVertexPointMap;
typedef CGAL::Cartesian_converter_property_map< SC::Point_3, PolyhedronVertexPointMap> ScPolyhedronVertexPointMap;

int main()
{
  Point_3 p(0.3, 0.3454, 0), q(0.6694, 0,0), r(0,0,0);
  Surface_mesh sm;
  Polyhedron_3 po;
  CGAL::make_triangle(p,q,r,sm);
  CGAL::make_triangle(p,q,r,po);

  VNMap smvnpm = sm.template add_property_map<vertex_descriptor, Vector_3>("v:normal", CGAL::NULL_VECTOR).first;

  CGAL::Polygon_mesh_processing::compute_vertex_normals(sm, smvnpm);

  ScSmVertexNormalMap scsmvnpm(smvnpm);

  SmVertexPointMap smvpm  = get(CGAL::vertex_point_t(), sm);
  ScSmVertexPointMap scsmvpm(smvpm);
  PolyhedronVertexPointMap povpm  = get(CGAL::vertex_point_t(), po);
  ScPolyhedronVertexPointMap scpovpm(povpm);

  // from <CGAL/boost/graph/IO/polygon_mesh_io.h>
  // https://doc.cgal.org/latest/BGL/group__PkgBGLIoFuncsPLY.html#ga959dcd88ca979d3b6b0806d883a0247f
  CGAL::IO::write_polygon_mesh("write_polygon_mesh_sm.ply", sm,
                               CGAL::parameters::stream_precision(17).use_binary_mode(true).vertex_point_map(scsmvpm).vertex_normal_map(scsmvnpm) );

{
  Surface_mesh osm;

  if(!CGAL::IO::read_polygon_mesh("write_polygon_mesh_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'write_polygon_mesh_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }else{
    for(const std::string& s : osm.properties<vertex_descriptor>()){
      std::cout << s << std::endl;
    }
  }
}

  PolyhedronVN vn;
  PolyhedronVNMap povnpm(vn);
  ScPolyhedronVertexNormalMap scpovnpm(povnpm);


  CGAL::Polygon_mesh_processing::compute_vertex_normals(po, povnpm);

  // ERROR this produces an invalid plyfile
  CGAL::IO::write_polygon_mesh("write_polygon_mesh_po.ply", po,
                               CGAL::parameters::stream_precision(17).use_binary_mode(true).vertex_point_map(scpovpm).vertex_normal_map(scpovnpm)) ;

{
  Surface_mesh osm;

  if(!CGAL::IO::read_polygon_mesh("write_polygon_mesh_po.ply", osm))
  {
    std::cerr << "Error: failed to read 'write_polygon_mesh_po.ply'" << std::endl;
    return EXIT_FAILURE;
  }else{
    for(const std::string& s : osm.properties<vertex_descriptor>()){
      std::cout << s << std::endl;
    }
  }
}
  //  OK
  // from #include <CGAL/boost/graph/IO/PLY.h>
  // https://doc.cgal.org/latest/BGL/group__PkgBGLIoFuncsPLY.html#ga959dcd88ca979d3b6b0806d883a0247f
  CGAL::IO::write_PLY("generic_write_PLY_sm.ply", sm, "generic write_PLY(Surface_mesh)",
                      CGAL::parameters::stream_precision(17).use_binary_mode(true).vertex_point_map(scsmvpm).vertex_normal_map(scsmvnpm));

{
  Surface_mesh osm;

  if(!CGAL::IO::read_polygon_mesh("generic_write_PLY_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'generic_write_PLY_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }else{
    for(const std::string& s : osm.properties<vertex_descriptor>()){
      std::cout << s << std::endl;
    }
  }
}

 // ERROR this produces an invalid plyfile
  CGAL::IO::write_PLY("generic_write_PLY_po.ply", po, "generic write_PLY(Polyhedron)",
                      CGAL::parameters::stream_precision(17).use_binary_mode(true).vertex_point_map(scpovpm));

{
  Surface_mesh osm;

  if(!CGAL::IO::read_polygon_mesh("generic_write_PLY_po.ply", osm))
  {
    std::cerr << "Error: failed to read 'generic_write_PLY_po.ply'" << std::endl;
    return EXIT_FAILURE;
  }else{
    for(const std::string& s : osm.properties<vertex_descriptor>()){
      std::cout << s << std::endl;
    }
  }
}
{
  // OK
  // from #include <CGAL/Surface_mesh/IO/PLY.h>
  // https://doc.cgal.org/latest/Surface_mesh/group__PkgSurfaceMeshIOFuncPLY.html#ga50f0e9f2b293855d2c7f1a62939cbe8d
  std::ofstream out("overloaded_write_PLY_sm.ply", std::ios::binary);
  CGAL::IO::set_binary_mode(out);
  CGAL::IO::write_PLY(out, sm, "overloaded_write_PLY(Surface_mesh)",CGAL::parameters::stream_precision(17).use_binary_mode(true)
    .vertex_point_map(scsmvpm).vertex_normal_map(scsmvnpm)  );
}
{
  Surface_mesh osm;

  if(!CGAL::IO::read_polygon_mesh("overloaded_write_PLY_sm.ply", osm))
  {
    std::cerr << "Error: failed to read 'overloaded_write_PLY_sm.ply'" << std::endl;
    return EXIT_FAILURE;
  }else{
    for(const std::string& s : osm.properties<vertex_descriptor>()){
      std::cout << s << std::endl;
    }
  }
}

  return 0;
}
