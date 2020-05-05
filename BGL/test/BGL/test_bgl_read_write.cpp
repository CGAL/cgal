#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

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

#include <CGAL/IO/VTK.h>

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

template<typename PointRange, typename PolygonRange>
void fill_soup(PointRange& points, PolygonRange& polygons)
{
  points.resize(4);
  polygons.resize(4);
  points[0] = Point(0, 0, 0);
  points[1] = Point(1, 1, 0);
  points[2] = Point(2, 0, 1);
  points[3] = Point(3, 0, 0);
  std::vector<std::size_t> poly(3);
  poly[0] = 0;
  poly[1] = 1;
  poly[2] = 2;
  polygons[0] = poly;
  poly[0] = 3;
  poly[1] = 1;
  poly[2] = 0;
  polygons[1] = poly;
  poly[0] = 3;
  poly[1] = 2;
  poly[2] = 1;
  polygons[2] = poly;
  poly[0] = 3;
  poly[1] = 0;
  poly[2] = 2;
  polygons[3] = poly;
}

template<typename Mesh>
void test_bgl_OFF(const char* filename)
{
  Mesh sm;
  std::ifstream in(filename);
  CGAL::read_polygon_mesh(in,sm);

  CGAL::write_OFF(std::cout, sm);
}

//todo check the result.
template<typename Mesh>
void test_bgl_OFF_with_np()
{
  Mesh fg;
  std::ifstream in("data/full.off");
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Kernel::Vector_3> >::type VertexNormalMap;
  VertexNormalMap vnm  = get(CGAL::dynamic_vertex_property_t<Kernel::Vector_3>(), fg);

  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<CGAL::Color> >::type VertexColorMap;
  VertexColorMap vcm  = get(CGAL::dynamic_vertex_property_t<CGAL::Color>(), fg);

  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Kernel::Point_2> >::type VertexTextureMap;
  VertexTextureMap vtm  = get(CGAL::dynamic_vertex_property_t<Kernel::Point_2>(), fg);

  typedef typename boost::property_map<Mesh, CGAL::dynamic_face_property_t<CGAL::Color> >::type FaceColorMap;
  FaceColorMap fcm  = get(CGAL::dynamic_face_property_t<CGAL::Color>(), fg);


  bool ok = CGAL::read_polygon_mesh(in, fg, CGAL::parameters::vertex_normal_map(vnm)
                 .vertex_color_map(vcm)
                 .vertex_texture_map(vtm)
                 .face_color_map(fcm));
  CGAL_assertion(ok);

  ok = CGAL::write_OFF(std::cout, fg, CGAL::parameters::vertex_normal_map(vnm)
                  .vertex_color_map(vcm)
                  .vertex_texture_map(vtm)
                  .face_color_map(fcm));
  CGAL_assertion(ok);
}

void test_soup_off(const char* filename)
{
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > polygons;
  std::ifstream in(filename);
  CGAL::read_OFF(in,points, polygons);
  CGAL::write_OFF(std::cout, points, polygons);
}

template<typename Mesh>
bool test_bgl_OBJ()
{
  Mesh fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  CGAL::write_OBJ(out, fg);
  std::istringstream in( out.str());
  fg.clear();
  CGAL::read_polygon_mesh(in, fg);
  CGAL_assertion(num_vertices(fg) == 4);
  CGAL_assertion(num_faces(fg) == 4);

  return true;
}

template<typename Mesh>
bool test_bgl_OBJ_with_np()
{
  Mesh fg, fg2;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Kernel::Vector_3> >::type VertexNormalMap;
  VertexNormalMap vnm  = get(CGAL::dynamic_vertex_property_t<Kernel::Vector_3>(), fg);
  VertexNormalMap vnm2  = get(CGAL::dynamic_vertex_property_t<Kernel::Vector_3>(), fg2);


  for(const auto& v : vertices(fg))
    put(vnm, v, Kernel::Vector_3(0.0,1.0,0.0));
  std::ostringstream out;
  CGAL::write_OBJ(out, fg, CGAL::parameters::vertex_normal_map(vnm));
  std::istringstream in( out.str());

  CGAL::read_polygon_mesh(in, fg2, CGAL::parameters::vertex_normal_map(vnm2));

  CGAL_assertion(num_vertices(fg2) == 4);
  CGAL_assertion(num_faces(fg2) == 4);
  typename boost::graph_traits<Mesh>::vertex_iterator vit, vit2;
  for( vit = vertices(fg).begin(), vit2 = vertices(fg2).begin();
       vit != vertices(fg).end(), vit2 != vertices(fg2).end();
       ++vit, ++vit2)
  {
    CGAL_assertion(get(vnm, *vit) == get(vnm2, *vit2));
  }

  return true;
}

void test_bgl_soup_obj()
{
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > polygons;
  std::ostringstream out;
  fill_soup(points, polygons);

  CGAL::write_OBJ(out, points, polygons);
  points.clear();
  polygons.clear();
  std::istringstream in( out.str());
  CGAL::read_OBJ(in,points, polygons);
  CGAL_assertion(points.size() == 4);
  CGAL_assertion(polygons.size() == 4);

}
#ifdef CGAL_USE_VTK
template<typename Mesh>
bool test_bgl_vtp(bool binary = false)
{
  Mesh fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);

  std::ofstream os("tetrahedron.vtp");
  CGAL::write_VTP(os, fg, CGAL::parameters::use_binary_mode(binary));
  if(!os)
  {
    std::cerr<<"vtp writing failed."<<std::endl;
    return false;
  }
  os.close();
  Mesh fg2;
  if(!CGAL::read_VTP("tetrahedron.vtp", fg2))
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

template<>
bool test_bgl_vtp<Polyhedron>(bool binary)
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
  CGAL::write_VTP(os, fg, CGAL::parameters::vertex_index_map(vid).use_binary_mode(binary));
  if(!os)
  {
    std::cerr<<"vtp writing failed."<<std::endl;
    return false;
  }
  os.close();
  Polyhedron fg2;
  if(!CGAL::read_VTP("tetrahedron.vtp", fg2))
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

bool test_soup_vtp(bool binary = false)
{
  //generate a test file
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > polys;
  fill_soup(points, polys);

  std::ofstream os("tetrahedron_soup.vtp");
  CGAL::write_VTP(os, points, polys, CGAL::parameters::use_binary_mode(binary));
  if(!os)
  {
    std::cerr<<"vtp writing failed."<<std::endl;
    return false;
  }
  os.close();


  std::vector<Point> soup_points;
  std::vector<std::vector<std::size_t> > soup_polygons;
  CGAL::read_VTP("tetrahedron_soup.vtp", soup_points, soup_polygons);
  if(4 != soup_points.size()
     || 4 != soup_polygons.size())
  {
    std::cerr<<"Coherence problem. Wrong number of vertices or faces."<<std::endl;
    return false;
  }

  return true;
}

#endif

template<class FaceGraph>
bool test_bgl_gocad()
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  CGAL::write_GOCAD(out, "tetrahedron", fg);
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  std::istringstream in( out.str());
  std::pair<std::string, std::string> cnn;
  CGAL::read_GOCAD(in, cnn, fg2);
  if(cnn.first != "tetrahedron"){
    std::cerr<<"reading error: tetrahedron != "<<cnn.first<<std::endl;
    return 1;
  }
  if( !cnn.second.empty()){
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
template<class FaceGraph>
bool test_bgl_gocad_with_np()
{

  typedef typename boost::property_map<FaceGraph,CGAL::vertex_point_t>::type VertexPointMap;
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  VertexPointMap vpm  = get(CGAL::vertex_point, fg);

  std::ostringstream out;
  CGAL::write_GOCAD(out, "tetrahedron", fg, CGAL::parameters::vertex_point_map(vpm));
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  VertexPointMap vpm2  = get(CGAL::vertex_point, fg2);
  std::istringstream in( out.str());
  std::pair<std::string, std::string> cnn;
  CGAL::read_GOCAD(in, cnn, fg2, CGAL::parameters::vertex_point_map(vpm2));
  if(cnn.first != "tetrahedron"){
    std::cerr<<"reading error: tetrahedron != "<<cnn.first<<std::endl;
    return 1;
  }
  if( !cnn.second.empty()){
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

bool test_soup_gocad()
{
  //generate a test file
  std::vector<Point> points(4);
  std::vector<std::vector<std::size_t> > polys(4);
  points[0] = Point(0, 0, 0);
  points[1] = Point(1, 1, 0);
  points[2] = Point(2, 0, 1);
  points[3] = Point(3, 0, 0);
  std::vector<std::size_t> poly(3);
  poly[0] = 0;
  poly[1] = 1;
  poly[2] = 2;
  polys[0] = poly;
  poly[0] = 3;
  poly[1] = 1;
  poly[2] = 0;
  polys[1] = poly;
  poly[0] = 3;
  poly[1] = 2;
  poly[2] = 1;
  polys[2] = poly;
  poly[0] = 3;
  poly[1] = 0;
  poly[2] = 2;
  polys[3] = poly;

  std::ofstream os("tetrahedron.TS");
  CGAL::write_GOCAD(os, points, polys);
  if(!os)
  {
    std::cerr<<"gocad writing failed."<<std::endl;
    return false;
  }
  os.close();


  std::vector<Point> soup_points;
  std::vector<std::vector<std::size_t> > soup_polygons;
  CGAL::read_GOCAD("tetrahedron.TS", soup_points, soup_polygons);
  if(4 != soup_points.size()
     || 4 != soup_polygons.size())
  {
    std::cerr<<"Coherence problem. Wrong number of vertices or faces."<<std::endl;
    return false;
  }

  return true;
}

template<class FaceGraph>
bool test_bgl_stl()
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;
  CGAL::write_STL(out, fg);
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  FaceGraph fg2;
  std::istringstream in(out.str());
  if(!CGAL::read_polygon_mesh(in, fg2)){
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

template<class FaceGraph>
bool test_bgl_PLY(bool binary = false)
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  std::ostringstream out;

  if(binary)
    CGAL::set_mode(out, CGAL::IO::BINARY);

  CGAL::write_PLY(out, fg, "hello");
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  std::istringstream in(out.str());

  if(binary)
    CGAL::set_mode(in, CGAL::IO::BINARY);

  fg.clear();
  if(!CGAL::read_polygon_mesh(in, fg)){
    std::cerr<<"Tetrahedron reading failed."<<std::endl;
    return false;
  }
  CGAL_assertion(num_vertices(fg) == 4);
  CGAL_assertion(num_faces(fg) == 4);
  return true;
}

template<typename FaceGraph>
bool test_bgl_PLY_with_np(bool binary)
{
  FaceGraph fg;
  CGAL::make_tetrahedron(Point(0, 0, 0), Point(1, 1, 0),
                         Point(2, 0, 1), Point(3, 0, 0), fg);
  typedef typename boost::property_map<FaceGraph, CGAL::dynamic_vertex_property_t<CGAL::Color> >::type VertexColorMap;
  typedef typename boost::property_map<FaceGraph, CGAL::dynamic_face_property_t<CGAL::Color> >::type FaceColorMap;
  FaceColorMap fcm  = get(CGAL::dynamic_face_property_t<CGAL::Color>(), fg);
  VertexColorMap vcm  = get(CGAL::dynamic_vertex_property_t<CGAL::Color>(), fg);

  auto fit = faces(fg).begin();
  put(fcm,*fit++,CGAL::Color(155,0,0));
  put(fcm,*fit++,CGAL::Color(0,155,0));
  put(fcm,*fit++,CGAL::Color(0,0,155));
  put(fcm,*fit++,CGAL::Color(155,0,155));

  auto vit = vertices(fg).begin();
  put(vcm,*vit++,CGAL::Color(255,0,0));
  put(vcm,*vit++,CGAL::Color(0,255,0));
  put(vcm,*vit++,CGAL::Color(0,0,255));
  put(vcm,*vit++,CGAL::Color(255,0,255));
  std::ostringstream out;

  if(binary)
    CGAL::set_mode(out, CGAL::IO::BINARY);

  CGAL::write_PLY(out, fg, "hello", CGAL::parameters::vertex_color_map(vcm)
                  .face_color_map(fcm));
  if(out.fail())
  {
    std::cerr<<"Tetrahedron writing failed."<<std::endl;
    return false;
  }
  std::istringstream in(out.str());

  if(binary)
    CGAL::set_mode(in, CGAL::IO::BINARY);

  fg.clear();

  VertexColorMap vcm2  = get(CGAL::dynamic_vertex_property_t<CGAL::Color>(), fg);
  FaceColorMap fcm2  = get(CGAL::dynamic_face_property_t<CGAL::Color>(), fg);
  if(!CGAL::read_polygon_mesh(in, fg, CGAL::parameters::vertex_color_map(vcm2).face_color_map(fcm2))){
    std::cerr<<"Tetrahedron reading failed."<<std::endl;
    return false;
  }
  CGAL_assertion(num_vertices(fg) == 4);
  CGAL_assertion(num_faces(fg) == 4);
  vit = vertices(fg).begin();
  CGAL_assertion(get(vcm2, *vit++) == CGAL::Color(255,0,0));
  CGAL_assertion(get(vcm2, *vit++) == CGAL::Color(0,255,0));
  CGAL_assertion(get(vcm2, *vit++) == CGAL::Color(0,0,255));
  CGAL_assertion(get(vcm2, *vit++) == CGAL::Color(255,0,255));

  fit = faces(fg).begin();
  CGAL_assertion(get(fcm2,*fit++)==CGAL::Color(155,0,0));
  CGAL_assertion(get(fcm2,*fit++)==CGAL::Color(0,155,0));
  CGAL_assertion(get(fcm2,*fit++)==CGAL::Color(0,0,155));
  CGAL_assertion(get(fcm2,*fit++)==CGAL::Color(155,0,155));


  return true;
}


int main(int argc, char** argv)
{
  const char* filename=(argc>1) ? argv[1] : "data/prim.off";
  // OFF
  test_bgl_OFF<Polyhedron>(filename);
  test_bgl_OFF<SM>(filename);
  test_bgl_OFF<LCC>(filename);
#ifdef CGAL_USE_OPENMESH
  test_bgl_OFF<OMesh>(filename);
#endif
  //polyhedron's overload doesn't care for any np that is not vpm
  test_bgl_OFF_with_np<Polyhedron>();
  test_bgl_OFF_with_np<SM>();
  test_bgl_OFF_with_np<LCC>();
#ifdef CGAL_USE_OPENMESH
  test_bgl_OFF_with_np<OMesh>();
#endif
    test_soup_off(filename);

  // OBJ
  test_bgl_OBJ<Polyhedron>();
  test_bgl_OBJ<SM>();
  test_bgl_OBJ<LCC>();
  test_bgl_soup_obj();
#ifdef CGAL_USE_OPENMESH
  test_bgl_OBJ<OMesh>();
#endif

  test_bgl_OBJ_with_np<Polyhedron>();
  test_bgl_OBJ_with_np<SM>();
  test_bgl_OBJ_with_np<LCC>();
#ifdef CGAL_USE_OPENMESH
  test_bgl_OBJ_with_np<OMesh>();
#endif


  //PLY
  if(!test_bgl_PLY<Polyhedron>())
    return 1;
  if(!test_bgl_PLY<Polyhedron>(true))
    return 1;
  if(!test_bgl_PLY<SM>())
    return 1;
  if(!test_bgl_PLY<SM>(true))
    return 1;

  if(!test_bgl_PLY_with_np<Polyhedron>(false))
    return 1;
  if(!test_bgl_PLY_with_np<Polyhedron>(true))
    return 1;
  if(!test_bgl_PLY_with_np<SM>(false))
    return 1;
  if(!test_bgl_PLY_with_np<SM>(true))
    return 1;


  // GOCAD
  if(!test_bgl_gocad<Polyhedron>())
    return 1;
  if(!test_bgl_gocad<SM>())
    return 1;
  if(!test_bgl_gocad<LCC>())
    return 1;
  if(!test_bgl_gocad_with_np<Polyhedron>())
    return 1;
  if(!test_bgl_gocad_with_np<SM>())
    return 1;
  if(!test_bgl_gocad_with_np<LCC>())
    return 1;
  if(!test_soup_gocad())
    return 1;

  // STL
  if(!test_bgl_stl<Polyhedron>())
    return 1;
  if(!test_bgl_stl<SM>())
    return 1;
  if(!test_bgl_stl<LCC>())
    return 1;

  // VTP
#ifdef CGAL_USE_VTK
  if(!test_bgl_vtp<Polyhedron>(false))
    return 1;
  if(!test_bgl_vtp<SM>(false))
    return 1;
  if(!test_bgl_vtp<LCC>(false))
    return 1;
  if(!test_soup_vtp(false))
    return 1;

  if(!test_bgl_vtp<Polyhedron>(true))
    return 1;
  if(!test_bgl_vtp<SM>(true))
    return 1;
  if(!test_bgl_vtp<LCC>(true))
    return 1;
  if(!test_soup_vtp(true))
    return 1;
#endif

  return EXIT_SUCCESS;
}
