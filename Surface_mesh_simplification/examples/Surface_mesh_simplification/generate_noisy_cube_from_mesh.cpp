#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_policies.h>
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <random>


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Bbox_3 Bbox_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef Surface_mesh::Vertex_index vertex_descriptor;
typedef Surface_mesh::Halfedge_index halfedge_desc;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace Subdivision_method_3 = CGAL::Subdivision_method_3;

namespace params = CGAL::parameters;

void varying_noise(Surface_mesh& surface_mesh, K kernel)
{
  // empty bounding box
  Bbox_3 bbox { };

  for (vertex_descriptor v : surface_mesh.vertices())
  {
    bbox += surface_mesh.point(v).bbox(); 
  }

  std::random_device r;
  std::default_random_engine e1(r());

  FT range = abs(bbox.xmax() - bbox.xmin());

  for (vertex_descriptor v : surface_mesh.vertices())
  {
    Point_3& p = surface_mesh.point(v);

    if (p.x() == bbox.xmax() && !(p.y() == bbox.ymax() ||
          p.y() == bbox.ymin() || p.z() == bbox.zmax() || p.z() == bbox.zmin()))
    {
      continue;
    }

    auto t = (p.x() - bbox.xmin()) / (range * 15); 

    std::normal_distribution norm {0.0, t};

    // only change y and z axes
    FT ry = FT(1.2) * norm(e1);
    FT rz = FT(1.2) * norm(e1);

    if (p.x() == bbox.xmax())
    {
      if (p.y() == bbox.ymin())
      {
        ry = -abs(ry);
      }
      else 
      {
        ry = abs(ry);
      }

      if (p.z() == bbox.zmin())
      {
        rz = -abs(rz);
      }
      else 
      {
        rz = abs(rz);
      }
    }

    p = {p.x(), p.y() + ry, p.z() + rz};
  }
}

// input mesh should be some sort of a meshed cube, unexpected things may occur when not
// the code then applies noise to the vertices depending on x coordinate. 
int main(int argc, char** argv)
{

  Surface_mesh surface_mesh { };
  K kernel { };

  if (argc != 3) 
  {
    std::cerr << "Please provide exactly one input and one output filename as arguments" << std::endl;
    return EXIT_FAILURE;
  }
  std::string filename = argv[1];

  std::ifstream is(filename);

  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  varying_noise(surface_mesh, kernel);

  CGAL::IO::write_polygon_mesh(argv[2], surface_mesh, 
      CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
