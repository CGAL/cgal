#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Graphics_scene_options.h>
#include <iostream>

#include "draw_surface_mesh.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3                Point;
typedef CGAL::Surface_mesh<Point>      Mesh;

struct Colored_faces_given_height:
  public CGAL::Graphics_scene_options<Mesh,
                                      typename Mesh::Vertex_index,
                                      typename Mesh::Edge_index,
                                      typename Mesh::Face_index>
{
  Colored_faces_given_height(const Mesh& sm)
  {
    if(sm.is_empty()) return;

    double m_min_z=0., m_max_z=0.;
    bool first=true;
    for(typename Mesh::Vertex_index vi: sm.vertices())
    {
      if(first)
      { m_min_z=sm.point(vi).z(); m_max_z=m_min_z; first=false; }
      else
      {
        m_min_z=(std::min)(m_min_z, sm.point(vi).z());
        m_max_z=(std::max)(m_max_z, sm.point(vi).z());
      }
    }

    this->colored_face=[](const Mesh &, typename Mesh::Face_index)->bool { return true; };

    this->face_color=[m_min_z, m_max_z]
      (const Mesh& sm, typename Mesh::Face_index fi)->CGAL::IO::Color
    {
      double res=0.;
      std::size_t n=0;
      for(typename Mesh::Vertex_index vi: vertices_around_face(sm.halfedge(fi), sm))
      {
        res+=sm.point(vi).z();
        ++n;
      }
      // Random color depending on the "height" of the facet
      CGAL::Random random(static_cast<unsigned int>(30*((res/n)-m_min_z)/(m_max_z-m_min_z)));
      return CGAL::get_random_color(random);
    };
  }
};

int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? CGAL::data_file_path(argv[1]) : CGAL::data_file_path("meshes/elephant.off");

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::draw(sm, Colored_faces_given_height(sm));

  return EXIT_SUCCESS;
}
