#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <nanosvg.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


int main(int argc, char** argv)
{
  std::string mesh_filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  std::string svg_filename = (argc > 2) ? std::string(argv[2])
    : CGAL::data_file_path("polylines_2/nano.svg");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input mesh." << std::endl;
    return 1;
  }

  NSVGimage* g_image = nsvgParseFromFile(svg_filename.c_str(), "px", 96.0f);
  if (g_image == NULL) {
    printf("Could not open SVG image.\n");
    return 1;
  }

  // extract control points
  std::vector< std::array<K::Point_2, 4> > bezier_curves;
  CGAL::Bbox_2 bb2;

  // in SVG's the y axis points downward, so we must take the opposite y coordinates
  for (NSVGshape* shape = g_image->shapes; shape != NULL; shape = shape->next)
  {
    for (NSVGpath* path = shape->paths; path != NULL; path = path->next)
    {
      CGAL::Bbox_2 path_bbox(path->bounds[0], -path->bounds[1],
                             path->bounds[2], -path->bounds[3]);
      bb2+=path_bbox;

      float* pts=path->pts;
      int npts=path->npts;

      for (int i=0; i<npts-1; i += 3)
      {
        bezier_curves.emplace_back();
        float* p = &pts[i*2];
        bezier_curves.back()[0]=K::Point_2(p[0],-p[1]);
        bezier_curves.back()[1]=K::Point_2(p[2],-p[3]);
        bezier_curves.back()[2]=K::Point_2(p[4],-p[5]);
        bezier_curves.back()[3]=K::Point_2(p[6],-p[7]);
      }
    }
  }

  nsvgDelete(g_image);

  std::cout << "#Bezier curves read: " << bezier_curves.size() << "\n";

  // convert control points to polar coordinates
  typename K::Point_2 center_2((bb2.xmax()+bb2.xmin())/2., (bb2.ymax()+bb2.ymin())/2.);
  double diag = std::sqrt( CGAL::square(bb2.xmin()-bb2.xmax()) + CGAL::square(bb2.xmin()-bb2.xmax()) );
  const double expected_diag = 0.45; // user parameter for scaling
  const double scaling = expected_diag/diag;

  //TODO: do the scaling at read time!

  std::vector<std::array<K::Vector_2, 4>> directions;
  std::vector<std::array<K::FT, 4>> lengths;
  directions.reserve(bezier_curves.size());
  lengths.reserve(bezier_curves.size());

  for (const std::array<K::Point_2, 4>& bezier  : bezier_curves)
  {
    std::vector<std::pair<double, double>> polar_coords =
      PMP::convert_polygon_to_polar_coordinates<K>(bezier, center_2);

    directions.emplace_back();
    lengths.emplace_back();

    assert(polar_coords.size()==4);

    for (int i=0;i<4; ++i)
    {
      lengths.back()[i] = scaling * polar_coords[i].first;
      directions.back()[i]=K::Vector_2(std::cos(polar_coords[i].second), std::sin(polar_coords[i].second));
    }
  }

  // trace bezier curves
  std::size_t nb_faces = faces(mesh).size();
  Mesh::Face_index f = *std::next(faces(mesh).begin(), (2154)%nb_faces);
  Face_location center(f, CGAL::make_array(0.3,0.3,0.4));

  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);

  std::vector< std::vector<Face_location> > res =
    PMP::trace_bezier_curves<K>(center, directions, lengths, 4, mesh, solver);

  // write result
  std::ofstream out("svg.polylines.txt");
  out << std::setprecision(17);
  for (const auto& b : res)
  {
    std::vector<K::Point_3> poly;
    poly.reserve(b.size());
    PMP::convert_path_to_polyline(b, mesh, std::back_inserter(poly));


    out << poly.size();
    for (const K::Point_3& p : poly)
      out << " " << p;
    out << "\n";
  }

  return 0;
}
