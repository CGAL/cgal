#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/helpers.h>

#if defined(CGAL_USE_OPENMESH)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#endif

#include <CGAL/boost/graph/io.h>
#include <CGAL/IO/STL/STL_writer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;

typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef CGAL::Surface_mesh<Point>                            SM;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, MyTraits>::type LCC;

#if defined(CGAL_USE_OPENMESH)

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OMesh;

#endif

template<typename Mesh>
void test_bgl_read_write(const char* filename)
{
  Mesh sm;
  std::ifstream in(filename);
  CGAL::read_off(in,sm);
  CGAL::write_off(std::cout, sm);
}

#ifdef CGAL_USE_VTK
template<typename Mesh>
bool test_bgl_vtp()
{
  Mesh fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);

  std::ofstream os("tetrahedron.vtp");
  CGAL::write_vtp(os, fg);
  if(!os)
  {
    std::cerr<<"vtp writing failed."<<std::endl;
    return false;
  }
  os.close();
  Mesh fg2;
  if(!CGAL::read_vtp("tetrahedron.vtp", fg2))
  {
    std::cerr<<"vtp reading failed."<<std::endl;
    return false;
  }
  if(num_vertices(fg) != num_vertices(fg2)
     || num_faces(fg) != num_faces(fg2))
  {
    std::cerr<<"Coherence problem. Wrong number of vertices or faces."<<std::endl;
    std::cerr<<num_faces(fg)<<" != "<<num_faces(fg2)<<std::endl;
    return false;
  }
  return true;
}

template<>
bool test_bgl_vtp<Polyhedron>()
{
  Polyhedron fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);

  typedef boost::property_map<Polyhedron, CGAL::dynamic_vertex_property_t<std::size_t> >::type VertexIdMap;
   VertexIdMap vid  = get(CGAL::dynamic_vertex_property_t<std::size_t>(), fg);
   std::size_t id = 0;
   for(auto v : vertices(fg))
     put(vid,v, id++);
  std::ofstream os("tetrahedron.vtp");
  CGAL::write_vtp(os, fg, CGAL::parameters::vertex_index_map(vid));
  if(!os)
  {
    std::cerr<<"vtp writing failed."<<std::endl;
    return false;
  }
  os.close();
  Polyhedron fg2;
  if(!CGAL::read_vtp("tetrahedron.vtp", fg2))
  {
    std::cerr<<"vtp reading failed."<<std::endl;
    return false;
  }
  if(num_vertices(fg) != num_vertices(fg2)
     || num_faces(fg) != num_faces(fg2))
  {
    std::cerr<<"Coherence problem. Wrong number of vertices or faces."<<std::endl;
    return false;
  }
  return true;
}
#endif

template<class FaceGraph>
bool test_gocad()
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  CGAL::write_gocad(fg, out, "tetrahedron");
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  std::istringstream in( out.str());
  std::string name, color;
  CGAL::read_gocad(fg2, in, name, color);
  if(name != "tetrahedron"){
    std::cerr<<"reading error: tetrahedron != "<<name<<std::endl;
    return 1;
  }
  if( !color.empty()){
    std::cerr<<"reading error: there should be no color."<<std::endl;
    return 1;
  }

  if(in.fail()){
    std::cerr<<"Tetrahedron reading failed."<<std::endl;
    return false;
  }

  if(num_vertices(fg2) != 4){
    std::cerr<<"Wrong number of vertices: 4 != "<<num_vertices(fg2)<<std::endl;
    return false;
  }

   if(num_faces(fg2) != 4)
  {
    std::cerr<<"Wrong number of faces: 4 != "<<num_faces(fg2)<<std::endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  const char* filename=(argc>1)?argv[1]:"data/prim.off";

  test_bgl_read_write<Polyhedron>(filename);
  test_bgl_read_write<SM>(filename);
  test_bgl_read_write<LCC>(filename);
  if(!test_gocad<Polyhedron>())
    return 1;
  if(!test_gocad<SM>())
    return 1;
  if(!test_gocad<LCC>())
    return 1;
#ifdef CGAL_USE_OPENMESH
  test_bgl_read_write<OMesh>(filename);
#endif
#ifdef CGAL_USE_VTK
  if(!test_bgl_vtp<Polyhedron>())
    return 1;
  if(!test_bgl_vtp<SM>())
    return 1;
  if(!test_bgl_vtp<LCC>())
    return 1;
#endif

  return EXIT_SUCCESS;
}
