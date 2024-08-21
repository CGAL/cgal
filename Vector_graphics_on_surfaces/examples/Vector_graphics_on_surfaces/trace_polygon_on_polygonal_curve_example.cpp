#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;

std::vector<Face_location>
get_supporting_curve(std::string filename,
                     const Mesh& mesh)
{
  std::ifstream in(filename);
  if (!in)
  {
    std::cerr << "ERROR: cannot open " << filename << "\n";
    exit(1);
  }

  int nbp;
  in >> nbp;
  std::cout << "Support polyline size: " << nbp << "\n";
  std::vector<K::Point_3> polyline(nbp);

  for (int i=0;i<nbp; ++i)
    in >> polyline[i];

  std::vector<Face_location> face_locations(nbp);

  using AABB_face_graph_primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using AABB_face_graph_traits = CGAL::AABB_traits_3<K, AABB_face_graph_primitive>;

  CGAL::AABB_tree<AABB_face_graph_traits> tree;
  PMP::build_AABB_tree(mesh, tree);

  //Remarks on snapping_tolerance of locate: it is actually a value used as is in barycentric coordinates with no relation to the size of the triangle
  //      and snapping is only done when at least one coordinate is negative (clamping would be a better name, snapping is what I'm doing below).
  //TODO: the following works only if we don't have a lfs smaller than snap_tol
  double snap_tol = 0.00001;
  for(int i=0; i<nbp; ++i)
  {
    K::Sphere_3 query(polyline[i], snap_tol*snap_tol);
    std::vector<Mesh::Face_index> primitives;
    tree.all_intersected_primitives(query, std::back_inserter(primitives));

    switch (primitives.size())
    {
      case 0:
        std::cerr <<"ERROR points too far!\n";
        exit(1);
      break;
      case 1:
        face_locations[i] = PMP::locate_in_face(polyline[i],primitives[0], mesh);
      break;
      case 2:
      {
        // TODO: we assume we don't have any boundary to it's clear it's on an edge
        Mesh::Halfedge_index h1 = halfedge(primitives[0], mesh),
                             h2 = halfedge(primitives[1], mesh);
        bool found=false;
        for (int j=0;j<3;++j)
        {
          for (int k=0;k<3; ++k)
          {
            if (h1==opposite(h2, mesh))
            {
              found = true;
              break;
            }
            h2=next(h2, mesh);
          }
          if (found) break;
          h1=next(h1, mesh);
        }
        CGAL_assertion(h1 == opposite(h2, mesh));

        const K::Point_3& src = mesh.point(source(h1, mesh));
        const K::Point_3& tgt= mesh.point(target(h1, mesh));
        double t = std::sqrt(CGAL::squared_distance(polyline[i], src)/CGAL::squared_distance(src,tgt));
        face_locations[i] = PMP::locate_on_halfedge(h1, t, mesh);
      }
      break;
      default:
      {
        std::vector<Mesh::Vertex_index> vrts;
        Mesh::Halfedge_index h=halfedge(primitives[0], mesh);
        vrts.push_back(source(h, mesh));
        vrts.push_back(target(h, mesh));
        vrts.push_back(target(next(h, mesh), mesh));
        std::sort(vrts.begin(), vrts.end());
        for(std::size_t k=1;k<primitives.size(); ++k)
        {
          std::vector<Mesh::Vertex_index> tmp;
          h=halfedge(primitives[k], mesh);
          for (int j=0;j<3;++j)
          {
            Mesh::Vertex_index v=target(h, mesh);
            if (std::binary_search(vrts.begin(), vrts.end(), v))
              tmp.push_back(v);
            h=next(h, mesh);
          }
          tmp.swap(vrts);
          std::sort(vrts.begin(), vrts.end());
        }
        CGAL_assertion(vrts.size()==1);

        int offset=0;
        h=halfedge(primitives[0], mesh);
        if (target(h, mesh)==vrts[0])
          offset=1;
        else
          if (target(next(h, mesh), mesh)==vrts[0])
            offset=2;
          else
            CGAL_assertion(source(h, mesh)==vrts[0]);
        std::array<double, 3> bary = CGAL::make_array(0.,0.,0.);
        bary[offset]=1.;
        face_locations[i] = std::make_pair(primitives[0], bary);
      }
    }
    // face_locations[i] = PMP::locate_with_AABB_tree(polyline[i], tree, mesh, CGAL::parameters::snapping_tolerance(snap_tol));
    for (int k=0;k<3;++k)
      if (face_locations[i].second[k]<0) face_locations[i].second[k]=0;
  }

  // DEBUG code to check that the polyline has correctly been converted to face locations (it's not a path as if vertex is hit we don't have the face continuity property)
  for (int i=0;i<nbp-1; ++i)
  {
    bool OK = PMP::locate_in_common_face(face_locations[i], face_locations[i+1], mesh);
    if (!OK)
    {
      std::cout << std::setprecision(17);
      std::cout << "issue with " << PMP::construct_point(face_locations[i], mesh) << " " <<  PMP::construct_point(face_locations[i+1], mesh) << "\n";
      std::cout << "   input = " << polyline[i] << " " <<  polyline[i+1] << "\n";
      std::cout << "           (" << face_locations[i].first << " , ("
                                  << face_locations[i].second[0] << ","
                                  << face_locations[i].second[1] << ","
                                  << face_locations[i].second[2]
                                  << "))  --  ("
                                  <<  face_locations[i+1].first << " , ("
                                  << face_locations[i+1].second[0] << ","
                                  << face_locations[i+1].second[1] << ","
                                  << face_locations[i+1].second[2] << "))\n";

    }
    CGAL_assertion(OK);
  }


  return face_locations;
}

int main(int argc, char** argv)
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  std::string support_filename = (argc > 2) ? std::string(argv[2])
    : CGAL::data_file_path("XXXXX");

  std::string filename_poly = (argc > 3) ? std::string(argv[3])
    : CGAL::data_file_path("XXXXX");

  std::size_t nb_copies =  atoi(argv[4]);

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


  // duplicate them a bit to get longer text while getting another one
  CGAL::Bbox_2 gbox;
  for (const auto& polygon : polygons)
    gbox+=CGAL::bbox_2(polygon.begin(), polygon.end());
  K::Vector_2 delta(gbox.xmax()-gbox.xmin(), 0);
  std::size_t nbpoly=polygons.size();
  polygons.reserve(nbpoly*(nb_copies+1));
  for (std::size_t c=0; c<nb_copies; ++c)
  {
    polygons.insert(polygons.end(), polygons.begin(), polygons.begin()+nbpoly);
    for (std::size_t ip=nbpoly*(c+1); ip<nbpoly*(c+2); ++ip)
    {
      for (K::Point_2& p : polygons[ip])
        p = p + delta * (c+1);
    }
  }


  std::ofstream debug("/tmp/input_polygons.polylines.txt");
  for (const auto& poly : polygons)
  {
    debug << poly.size();
    for (const auto &p : poly)
      debug << " " << p << " 0";
    debug << std::endl;
  }
  debug.close();

  std::cout << "Total number of polygons (after copy): " << polygons.size() << "\n";

  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);

  // get supporting curve
  std::vector<Face_location> supporting_curve = get_supporting_curve(support_filename, mesh);
  if (supporting_curve.empty()) return 1;


  std::ofstream support_out("support.polylines.txt");

  std::vector<K::Point_3> support_poly;
  support_poly.reserve(supporting_curve.size());
  PMP::convert_path_to_polyline(supporting_curve, mesh, std::back_inserter(support_poly));

  support_out << std::setprecision(17) << support_poly.size();
  for (auto p : support_poly)
    support_out << " " << p;
  support_out << "\n";
  support_out.close();
  std::cout <<"supporting_curve generated!\n";

  // convert polygons to polar coordinates
  const double expected_height = 0.025; // user parameter for scaling
  const double scaling = expected_height/(gbox.ymax()-gbox.ymin());


  std::ofstream out("label_on_curve.polylines.txt");
  out << std::setprecision(17);

  std::vector<std::vector<Face_location>> polygons_3;
  polygons_3 = PMP::trace_geodesic_label_along_curve<K>(supporting_curve, polygons, scaling, 0., false, mesh, solver);

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
