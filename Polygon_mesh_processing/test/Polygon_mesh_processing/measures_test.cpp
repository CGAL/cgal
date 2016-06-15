#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

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

  halfedge_descriptor border_he;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
  {
    if (is_border(h, pmesh))
    {
      border_he = h;
      break;
    }
  }
  FT border_l = PMP::face_border_length(border_he, pmesh);
  std::cout << "length of hole border = " << border_l << std::endl;

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
  
  FT mesh_area = PMP::area(pmesh);
  std::cout << "mesh area = " << mesh_area << std::endl;
  
  FT mesh_area_np = PMP::area(pmesh,
    PMP::parameters::geom_traits(K()));
  std::cout << "mesh area (NP) = " << mesh_area_np << std::endl;


  CGAL::Bbox_3 bb = PMP::bbox_3(pmesh);
  std::cout << "bbox x[" << bb.xmin() << "; " << bb.xmax() << "]" << std::endl;
  std::cout << "     y[" << bb.ymin() << "; " << bb.ymax() << "]" << std::endl;
  std::cout << "     z[" << bb.zmin() << "; " << bb.zmax() << "]" << std::endl;

}

template <typename Polyhedron, typename K>
void test_polyhedron(const char* filename)
{
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
