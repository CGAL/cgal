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

std::vector<Face_location>
get_supporting_curve(std::string svg_filename,
                     const Mesh& mesh,
                     Face_location center, const PMP::Dual_geodesic_solver<double>& solver)
{
  std::vector<Face_location> res;

  NSVGimage* g_image = nsvgParseFromFile(svg_filename.c_str(), "px", 96.0f);
  if (g_image == NULL) {
    printf("Could not open SVG image.\n");
    return res;
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

  CGAL_assertion_code(K::Point_2 prev = bezier_curves[0][0];)
  for (const std::array<K::Point_2, 4>& bezier  : bezier_curves)
  {
    CGAL_assertion(bezier[0]==prev);
    CGAL_assertion_code(prev=bezier[3];)
    std::vector<std::pair<double, double>> polar_coords =
      PMP::convert_polygon_to_polar_coordinates<K>(bezier, center_2);

    int start_id=directions.empty()?0:1;

    directions.emplace_back();
    lengths.emplace_back();

    assert(polar_coords.size()==4);

    for (int i=start_id;i<4; ++i)
    {
      lengths.back()[i] = scaling * polar_coords[i].first;
      directions.back()[i]=K::Vector_2(std::cos(polar_coords[i].second), std::sin(polar_coords[i].second));
    }
    if (start_id==1)
    {
      lengths.back()[0] = lengths[lengths.size()-2][3];
      directions.back()[0] = directions[directions.size()-2][3];
    }
  }

  // trace bezier curves
  std::vector< std::vector<Face_location> > resi =
    PMP::trace_bezier_curves<K>(center, directions, lengths, 6, mesh, solver);

  //TODO: here we assume that curves are parameterized in the same order and are consecutives
  for (const std::vector<Face_location>& r : resi)
    res.insert(res.end(), r.begin(), std::prev(r.end()));
  res.push_back(resi.back().back());

  std::reverse(res.begin(), res.end()); // TODO: should be an option!

  return res;
}

int main(int argc, char** argv)
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  std::string svg_filename = (argc > 2) ? std::string(argv[2])
    : CGAL::data_file_path("polylines_2/Archimedean_spiral.svg");

  std::string filename_poly = (argc > 3) ? std::string(argv[3])
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

  K::Point_3 center_pt = PMP::construct_point(center, mesh);
  std::cout << "center = " << center_pt << "\n";
  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);


  // get supporting curve
  std::vector<Face_location> supporting_curve = get_supporting_curve(svg_filename, mesh, center, solver);
  if (supporting_curve.empty()) return 1;

  std::cout <<"supporting_curve generated!\n";
  std::ofstream support_out("svg_support.polylines.txt");

  std::vector<K::Point_3> support_poly;
  support_poly.reserve(supporting_curve.size());
  PMP::convert_path_to_polyline(supporting_curve, mesh, std::back_inserter(support_poly));

  support_out << std::setprecision(17) << support_poly.size();
  for (auto p : support_poly)
    support_out << " " << p;
  support_out << "\n";
  support_out.close();

  // convert polygons to polar coordinates
  typename K::Point_2 center_2((bb2.xmax()+bb2.xmin())/2., (bb2.ymax()+bb2.ymin())/2.);
  double diag = std::sqrt( CGAL::square(bb2.xmin()-bb2.xmax()) + CGAL::square(bb2.xmin()-bb2.xmax()) );
  const double expected_diag = 2.1; // user parameter for scaling
  const double scaling = expected_diag/diag;


  std::ofstream out("label_on_curve.polylines.txt");
  out << std::setprecision(17);

  std::vector<std::vector<Face_location>> polygons_3;
  polygons_3 = PMP::trace_geodesic_label_along_curve<K>(supporting_curve, polygons, scaling, 0., true, mesh, solver);

  for (const auto& polygon_path : polygons_3)
  {
    std::vector<K::Point_3> poly;
    poly.reserve(polygon_path.size());
    PMP::convert_path_to_polyline(polygon_path, mesh, std::back_inserter(poly));

    out << poly.size();
    for (auto p : poly)
      out << " " << p;
    out << std::endl;
  }

  return 0;
}
