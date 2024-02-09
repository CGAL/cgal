#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#if 0
#include <CGAL/Polygon_mesh_processing/walk_to_select.h>
#else

#endif

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


int main(int argc, char** argv)
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  std::string filename_poly = (argc > 2) ? std::string(argv[2])
    : CGAL::data_file_path("XXXXX");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<std::vector<K::Point_2>> polygons;
  std::ifstream in(filename_poly);
  if (!in)
  {
    std::cerr << "Error cannot open " << filename_poly << "\n";
    return 1;
  }

  int nb_pt;
  K::Point_3 pt;
  CGAL::Bbox_2 bb2;
  while (in >> nb_pt)
  {
    polygons.emplace_back();
    polygons.back().reserve(nb_pt-1);
    for (int i=0; i<nb_pt-1; ++i)
    {
      if (!in)
      {
        std::cerr << "Error reading input polygons\n";
        return 1;
      }
      in >> pt;
      polygons.back().emplace_back(pt.x(), pt.y());
      bb2+=polygons.back().back().bbox();
    }
    in >> pt;
    if (!in)
    {
      std::cerr << "Error reading input polygons\n";
      return 1;
    }
    // check if last point is duplicated
    if (polygons.back().back().x()!=pt.x() || polygons.back().back().y()!=pt.y())
    {
      polygons.back().emplace_back(pt.x(), pt.y());
      bb2+=polygons.back().back().bbox();
    }
    if (!in) break;
  }

  std::cout << polygons.size() << " polygons read\n";

  // tracing center
  std::size_t nb_faces = faces(mesh).size();
  Mesh::Face_index f = *std::next(faces(mesh).begin(), (2154)%nb_faces);
  Face_location center(f, CGAL::make_array(0.3,0.3,0.4));

  // convert polygons to polar coordinates
  typename K::Point_2 center_2((bb2.xmax()+bb2.xmin())/2., (bb2.ymax()+bb2.ymin())/2.);
  double diag = std::sqrt( CGAL::square(bb2.xmin()-bb2.xmax()) + CGAL::square(bb2.xmin()-bb2.xmax()) );
  const double expected_diag = 0.45; // user parameter for scaling
  const double scaling = expected_diag/diag;

  std::ofstream out("geodesic_polygon.polylines.txt");
  out << std::setprecision(17);

  K::Point_3 center_pt = PMP::construct_point(center, mesh);
  std::cout << "center = " << center_pt << "\n";
  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);

  for (const std::vector<K::Point_2>& polygon : polygons)
  {
    std::vector<std::pair<double, double>> polar_coords =
      PMP::convert_polygon_to_polar_coordinates<K>(polygon, center_2);

    std::vector<K::Vector_2> directions;
    std::vector<K::FT> lens;
    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second), std::sin(coord.second));
    }

    // last point is duplicated
    std::vector<K::Point_3> out_polygon = PMP::trace_geodesic_polygon<K>(center,directions,lens,mesh, solver);
    out << out_polygon.size();
    for (auto p : out_polygon)
      out << " " << p;
    out << std::endl;
  }

  // second method
  out.close();
  out.open("geodesic_polygons.polylines.txt");
  out << std::setprecision(17);

  std::vector<std::vector<K::Point_3>> polygons_3
    = PMP::trace_geodesic_polygons<K>(center, polygons, scaling, mesh, solver);

  for (const auto& polygon : polygons_3)
  {
    out << polygon.size();
    for (auto p : polygon)
      out << " " << p;
    out << std::endl;
  }

  // third method
  out.close();
  out.open("geodesic_label.polylines.txt");
  out << std::setprecision(17);

  polygons_3.clear();
  polygons_3 = PMP::trace_geodesic_label<K>(center, polygons, scaling, mesh, solver);

  for (const auto& polygon : polygons_3)
  {
    out << polygon.size();
    for (auto p : polygon)
      out << " " << p;
    out << std::endl;
  }

  return 0;
}
