#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Bbox_3.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>

#include <boost/foreach.hpp>

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
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
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
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
  {
    if (is_border(h, pmesh) || is_border(opposite(h, pmesh), pmesh))
      continue;
    else
    {
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
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f, pmesh), pmesh))
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
    PMP::parameters::geom_traits(K()));
  std::cout << "mesh area (NP) = " << mesh_area_np << std::endl;
  assert(mesh_area_np > 0);

  CGAL::Bbox_3 bb = PMP::bbox(pmesh);
  std::cout << "bbox x[" << bb.xmin() << "; " << bb.xmax() << "]" << std::endl;
  std::cout << "     y[" << bb.ymin() << "; " << bb.ymax() << "]" << std::endl;
  std::cout << "     z[" << bb.zmin() << "; " << bb.zmax() << "]" << std::endl;

  CGAL::Bbox_3 bb_v;
  BOOST_FOREACH(vertex_descriptor vd, vertices(pmesh))
    bb_v+=PMP::vertex_bbox(vd, pmesh);

  CGAL::Bbox_3 bb_f;
  BOOST_FOREACH(face_descriptor fd, faces(pmesh))
    bb_f+=PMP::face_bbox(fd, pmesh);

  CGAL::Bbox_3 bb_e;
  BOOST_FOREACH(edge_descriptor ed, edges(pmesh))
    bb_e+=PMP::edge_bbox(ed, pmesh);

  assert(bb==bb_v);
  assert(bb==bb_f);
  assert(bb==bb_e);
}

template <typename Polyhedron, typename K>
void test_polyhedron(const char* filename)
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
void test_closed_surface_mesh(const char* filename)
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

int main(int argc, char* argv[])
{
  const char* filename_polyhedron =
    (argc > 1) ? argv[1] : "data/mech-holes-shark.off";
  test_polyhedron<CGAL::Polyhedron_3<Epic>,Epic>(filename_polyhedron);
  test_polyhedron<CGAL::Polyhedron_3<Epec>,Epec>(filename_polyhedron);

  const char* filename_surface_mesh =
    (argc > 1) ? argv[1] : "data/elephant.off";
  test_closed_surface_mesh<CGAL::Surface_mesh<Epic::Point_3>,Epic>(filename_surface_mesh);
  test_closed_surface_mesh<CGAL::Surface_mesh<Epec::Point_3>,Epec>(filename_surface_mesh);

  std::cerr << "All done." << std::endl;
  return 0;
}
