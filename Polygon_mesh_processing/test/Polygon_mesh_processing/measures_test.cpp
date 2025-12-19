#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Bbox_3.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;


template<typename Mesh, typename K>
void test_pmesh(const Mesh& pmesh)
{
  typedef typename K::FT FT;

  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

  bool has_border = false;
  halfedge_descriptor border_he;
  for(halfedge_descriptor h : halfedges(pmesh))
  {
    if (is_border(h, pmesh))
    {
      border_he = h;
      has_border = true;
      break;
    }
  }
  FT border_l = PMP::face_border_length(border_he, pmesh);
  std::cout << "length of hole border = " << border_l << std::endl;
  if (has_border)
    assert(border_l > 0);

  face_descriptor valid_patch_face;
  unsigned int count = 0;
  for(halfedge_descriptor h : halfedges(pmesh))
  {
    if (is_border(h, pmesh) || is_border(opposite(h, pmesh), pmesh))
      continue;
    else
    {
       FT edge_length = PMP::edge_length(h, pmesh);
       FT squared_edge_length = PMP::squared_edge_length(h, pmesh);
       std::cout << "squared edge length = " << squared_edge_length << std::endl;
       std::cout << "edge length = " << edge_length << std::endl;

      FT  squared_face_area = PMP::squared_face_area(face(h, pmesh), pmesh);
      std::cout << "squared face area = " << squared_face_area << std::endl;
      FT face_area = PMP::face_area(face(h, pmesh), pmesh);
      std::cout << "face area = " << face_area << std::endl;
      assert(face_area > 0);

      if(++count == 20)
      {
        valid_patch_face = face(h, pmesh);
        break;
      }
    }
  }

  std::list<face_descriptor> patch;
  patch.push_back(valid_patch_face);
  while (patch.size() < 5)
  {
    face_descriptor f = patch.front();
    patch.pop_front();
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      if (boost::graph_traits<Mesh>::null_halfedge() != opposite(h, pmesh))
        patch.push_back(face(opposite(h, pmesh), pmesh));
      patch.push_back(f);
    }
    if (patch.front() == valid_patch_face)
      break;//back to starting point
  }

  FT patch_area = PMP::area(patch, pmesh);
  std::cout << "patch area = " << patch_area << std::endl;
  assert(patch_area > 0);

  FT mesh_area = PMP::area(pmesh);
  std::cout << "mesh area = " << mesh_area << std::endl;
  assert(mesh_area >= patch_area);

  FT mesh_area_np = PMP::area(pmesh,
    CGAL::parameters::geom_traits(K()));
  std::cout << "mesh area (NP) = " << mesh_area_np << std::endl;
  assert(mesh_area_np > 0);

  std::pair<halfedge_descriptor, FT> res = PMP::longest_border(pmesh);
  if(res.first == boost::graph_traits<Mesh>::null_halfedge()){
    std::cout << "mesh has no border" << std::endl;
  } else {
    std::cout << "longest border has length = " << res.second <<std::endl;
  }

  CGAL::Bbox_3 bb = PMP::bbox(pmesh);
  std::cout << "bbox x[" << bb.xmin() << "; " << bb.xmax() << "]" << std::endl;
  std::cout << "     y[" << bb.ymin() << "; " << bb.ymax() << "]" << std::endl;
  std::cout << "     z[" << bb.zmin() << "; " << bb.zmax() << "]" << std::endl;

  CGAL::Bbox_3 bb_v;
  for(vertex_descriptor vd : vertices(pmesh))
    bb_v+=PMP::vertex_bbox(vd, pmesh);

  CGAL::Bbox_3 bb_f;
  for(face_descriptor fd : faces(pmesh))
    bb_f+=PMP::face_bbox(fd, pmesh);

  CGAL::Bbox_3 bb_e;
  for(edge_descriptor ed : edges(pmesh))
    bb_e+=PMP::edge_bbox(ed, pmesh);

  assert(bb==bb_v);
  assert(bb==bb_f);
  assert(bb==bb_e);
}

template <typename Polyhedron, typename K>
void test_polyhedron(const std::string filename)
{
  std::cout << "Test Polyhedron " << filename
    << " with Kernel " << typeid(K).name() << std::endl;

  //run test for a Polyhedron
  Polyhedron poly; // file should contain oriented polyhedron
  std::ifstream input(filename);

  if (!input || !(input >> poly))
  {
    std::cerr << "Error: cannot read Polyhedron : " << filename << "\n";
    assert(!poly.empty());
    assert(false);
    return;
  }

  test_pmesh<Polyhedron, K>(poly);
}

template <typename Surface_mesh, typename K>
void test_closed_surface_mesh(const std::string filename)
{
  std::cout << "Test Surface_mesh " << filename
    << " with Kernel " << typeid(K).name() << std::endl;
  Surface_mesh sm;
  std::ifstream input(filename);

  if (!input || !(input >> sm))
  {
    std::cerr << "Error: cannot read Surface mesh : " << filename << "\n";
    assert(sm.number_of_vertices() > 0);
    assert(false);
    return;
  }

  test_pmesh<Surface_mesh, K>(sm);

  typename K::FT vol = PMP::volume(sm);
  std::cout << "volume = " << vol << std::endl;
  assert(vol > 0);

}


template <typename Surface_mesh, typename K>
void test_centroid(const std::string filename)
{
  std::cout << "Test Surface_mesh " << filename
    << " with Kernel " << typeid(K).name() << std::endl;
  Surface_mesh sm;
  std::ifstream input(filename);
  input >> sm;

  typename K::Point_3 p = PMP::centroid(sm);

  // For data/elephant.off
  // compare with centroid of 1.000.000 points inside the mesh:
  //  0.00772887 -0.134923 0.011703
  assert (p.x() > 0.007 && p.x() < 0.008);
  assert (p.y() > -0.14 && p.y() < -0.13);
  assert (p.z() > 0.01 && p.z() < 0.02);

  typename K::Vector_3 v(10,20,30);
  for(typename boost::graph_traits<Surface_mesh>::vertex_descriptor vd : vertices(sm)){
    sm.point(vd) = sm.point(vd) + v;
  }

  p = PMP::centroid(sm);
  p = p - v;
  assert (p.x() > 0.007 && p.x() < 0.008);
  assert (p.y() > -0.14 && p.y() < -0.13);
  assert (p.z() > 0.01 && p.z() < 0.02);


}

template <typename PolygonMesh1, typename PolygonMesh2 >
void test_compare()
{
  typedef typename boost::graph_traits<PolygonMesh1>::face_descriptor        face_descriptor1;
  typedef typename boost::graph_traits<PolygonMesh2>::face_descriptor        face_descriptor2;
  namespace PMP = CGAL::Polygon_mesh_processing;

  PolygonMesh1 mesh1;
  PolygonMesh2 mesh2;
  std::vector<std::pair<face_descriptor1, face_descriptor2> > common;
  common.clear();
  std::vector<face_descriptor1> m1_only;
  std::vector<face_descriptor2> m2_only;
  /*************************
   * triangulated and open *
   * **********************/

  std::ifstream input("data/tri1.off");
  if(! (input >> mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  input.open("data/tri2.off");
  if(! (input >> mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only), CGAL::parameters::default_values(), CGAL::parameters::default_values());
  assert(common.size() == 7);
  assert(m1_only.size() == 11);
  assert(m2_only.size() == 11);
  /*************************
   **** quad and closed ****
   * **********************/
  CGAL::clear(mesh1);
  CGAL::clear(mesh2);
  common.clear();
  m1_only.clear();
  m2_only.clear();
  input.open(CGAL::data_file_path("meshes/cube_quad.off"));
  if(! (input >> mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  input.open("data/cube_quad2.off");
  if(! (input >> mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  std::unordered_map<face_descriptor1, std::size_t> fim1;
  std::unordered_map<face_descriptor2, std::size_t> fim2;
  std::size_t id = 0;
  for(const auto& f : faces(mesh1))
  {
    fim1.insert(std::make_pair(f, id++));
  }
  id = 0;
  for(const auto& f : faces(mesh2))
  {
    fim2.insert(std::make_pair(f, id++));
  }
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only), CGAL::parameters::default_values(), CGAL::parameters::default_values());
  assert(common.size() == 3);
  assert(m1_only.size() == 3);
  assert(fim1[m1_only[0]] == 0);
  assert(fim1[m1_only[1]] == 3);
  assert(fim1[m1_only[2]] == 4);
  assert(m2_only.size() == 3);
  assert(fim2[m2_only[0]] == 0);
  assert(fim2[m2_only[1]] == 3);
  assert(fim2[m2_only[2]] == 4);
  /*************************
   **** tri and hole****
   * **********************/
  CGAL::clear(mesh1);
  CGAL::clear(mesh2);
  common.clear();
  m1_only.clear();
  m2_only.clear();
  input.open("data/tri1.off");
  if(! (input >> mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  input.open("data/tri1-hole.off");
  if(! (input >> mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only), CGAL::parameters::default_values(), CGAL::parameters::default_values());
  assert(common.size() == 17);
  assert(m1_only.size() == 1);
  assert(m2_only.size() == 0);

  /*************************
   **** tri and orient****
   * **********************/
  CGAL::clear(mesh1);
  CGAL::clear(mesh2);
  common.clear();
  m1_only.clear();
  m2_only.clear();
  input.open("data/tri2.off");
  if(! (input >> mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  input.open("data/tri2-out.off");
  if(! (input >> mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    assert (false);
    return;
  }
  input.close();
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only),
                      std::back_inserter(m2_only));
  assert(common.size() == 0);
  assert(m1_only.size() == 18);
  assert(m2_only.size() == 18);

}

int main(int argc, char* argv[])
{
  const std::string filename_polyhedron =
    (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mech-holes-shark.off");
  test_polyhedron<CGAL::Polyhedron_3<Epic>,Epic>(filename_polyhedron);
  test_polyhedron<CGAL::Polyhedron_3<Epec>,Epec>(filename_polyhedron);

  const std::string filename_surface_mesh =
    (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  test_closed_surface_mesh<CGAL::Surface_mesh<Epic::Point_3>,Epic>(filename_surface_mesh);
  test_closed_surface_mesh<CGAL::Surface_mesh<Epec::Point_3>,Epec>(filename_surface_mesh);

  // It won't work with Epec for large meshes as it builds up a deep DAG
  // leading to a stackoverflow when the destructor is called.
  test_centroid<CGAL::Surface_mesh<Epic::Point_3>,Epic>(filename_surface_mesh);
  test_compare<CGAL::Polyhedron_3<Epic>, CGAL::Surface_mesh<Epic::Point_3> >();
  test_compare<CGAL::Polyhedron_3<Epec>, CGAL::Surface_mesh<Epec::Point_3> >();
  test_compare<CGAL::Surface_mesh<Epic::Point_3>, CGAL::Polyhedron_3<Epic> >();
  test_compare<CGAL::Surface_mesh<Epec::Point_3>, CGAL::Polyhedron_3<Epec> >();
  std::cerr << "All done." << std::endl;
  return 0;
}
