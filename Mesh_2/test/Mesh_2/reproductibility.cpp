//#define CGAL_MESH_2_DEBUG_BAD_FACES
//#define CGAL_MESH_2_DEBUG_CLUSTERS

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;
using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
using Mesher =  CGAL::Delaunay_mesher_2<CDT, Criteria>;

using Vertex_handle = CDT::Vertex_handle;
using Point = CDT::Point;

int main(int, char**)
{
  auto triangulate = [](int index)
  {
    CDT cdt;

    Vertex_handle va = cdt.insert(Point(-0.74397572, -0.54545455));
    Vertex_handle vb = cdt.insert(Point(-0.13526831, -1));
    Vertex_handle vc = cdt.insert(Point(0.067634156, -1));
    Vertex_handle vd = cdt.insert(Point(0.33817078, -0.54545455));
    Vertex_handle ve = cdt.insert(Point(0.74397572, 0.27272727));
    Vertex_handle vf = cdt.insert(Point(0.74397572, 0.54545455));
    Vertex_handle vg = cdt.insert(Point(0.067634156, 1));
    Vertex_handle vh = cdt.insert(Point(-0.13526831, 1));
    Vertex_handle vi = cdt.insert(Point(-0.74397572, -0.18181818));

    cdt.insert_constraint(va, vb);
    cdt.insert_constraint(vb, vc);
    cdt.insert_constraint(vc, vd);
    cdt.insert_constraint(vd, ve);
    cdt.insert_constraint(ve, vf);
    cdt.insert_constraint(vf, vg);
    cdt.insert_constraint(vg, vh);
    cdt.insert_constraint(vh, vi);
    cdt.insert_constraint(vi, va);

    const std::vector<Point> points{
      Point(0.65605132, 0.43821259),
      Point(0.23073753, -0.4476739),
      Point(-0.037496007, -0.93636364),
      Point(-0.00095596601, 0.88181818),
      Point(-0.62452925, -0.30720903),
      Point(-0.69663181, -0.45045525),
    };

    cdt.insert(points.cbegin(), points.cend());

    std::cout << "Meshing: " << index << std::endl;

    std::cout << "Number of vertices before: " << cdt.number_of_vertices() << std::endl;

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125, 0.05*std::sqrt(2)));

    mesher.refine_mesh();

    std::cout << "Number of vertices after: " << cdt.number_of_vertices() << std::endl;

    std::stringstream ss;
    ss << cdt;

    return ss.str();
  };

  const std::string ref_cdts = triangulate(0);

  for (int i = 1; i < 20; ++i)
  {
    const std::string cdts = triangulate(i);
    if (ref_cdts != cdts)
      return 1;
  }

  return 0;
}
