#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

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
    if (polygon.front()==polygon.back()) polar_coords.pop_back();

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
    std::vector<Face_location> out_polygon_path = PMP::trace_geodesic_polygon<K>(center,directions,lens,mesh, solver);
    std::vector<K::Point_3> poly;
    poly.reserve(out_polygon_path.size());
    PMP::convert_path_to_polyline(out_polygon_path, mesh, std::back_inserter(poly));

    out << poly.size();
    for (auto p : poly)
      out << " " << p;
    out << std::endl;
  }

  // second method
  out.close();
  out.open("geodesic_polygons.polylines.txt");
  out << std::setprecision(17);

  std::vector<std::vector<Face_location>> polygons_3
    = PMP::trace_geodesic_polygons<K>(center, polygons, scaling, mesh, solver);

  for (const auto& polygon : polygons_3)
  {
    out << polygon.size();
    for (auto p : polygon)
      out << " " << PMP::construct_point(p, mesh);
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
      out << " " << PMP::construct_point(p, mesh);
    out << std::endl;
  }

  // now refine the input mesh
  std::vector<Mesh::Halfedge_index> cst_hedges;
  auto vnm = mesh.add_property_map<Mesh::Vertex_index, K::Vector_3>("vnm", K::Vector_3(0,0,0)).first;
  auto fnm = mesh.add_property_map<Mesh::Face_index, K::Vector_3>("fnm", K::Vector_3(0,0,0)).first;
  using VNM = decltype(vnm);

  PMP::compute_normals(mesh, vnm, fnm);
  PMP::refine_mesh_along_paths<K>(polygons_3, mesh, vnm, fnm, std::back_inserter(cst_hedges));

  std::ofstream("mesh_refined.off") << std::setprecision(17) << mesh;

  std::ofstream cst_edges("refinement_edges.polylines.txt");
  cst_edges.precision(17);
  for (Mesh::Halfedge_index h : cst_hedges)
    cst_edges << "2 " << mesh.point(source(h,mesh)) << " " << mesh.point(target(h,mesh)) << "\n";


  auto ecm = mesh.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  for (Mesh::Halfedge_index h : cst_hedges)
    ecm[edge(h, mesh)]=true;


  // face index for doing flood fill and mark inside-out
  Mesh::Face_index out_face(2612);
  std::vector<int> in_out(num_faces(mesh), -1);

  bool inorout=false;
  std::vector<Mesh::Face_index> queue, next_queue;
  queue.push_back(out_face);

  while(!queue.empty())
  {
    Mesh::Face_index f = queue.back();
    queue.pop_back();
    if (in_out[f]==-1)
    {
      in_out[f]=inorout?1:0;
      Mesh::Halfedge_index h=halfedge(f, mesh);
      for (int i=0; i<3; ++i)
      {
        Mesh::Face_index nf = face(opposite(h, mesh), mesh);
        if (nf!=boost::graph_traits<Mesh>::null_face() && in_out[nf]==-1)
        {
          if (ecm[edge(h,mesh)])
            next_queue.push_back(nf);
          else
            queue.push_back(nf);
        }
        h=next(h, mesh);
      }
    }
    if (queue.empty())
    {
      queue.swap(next_queue);
      inorout=!inorout;
    }
  }



  struct Visitor
    : public PMP::Corefinement::Default_visitor<Mesh>
  {
    VNM vnm;
    Visitor(VNM vnm) : vnm(vnm) {}

    std::vector<std::pair<Mesh::Halfedge_index, Mesh::Halfedge_index> > hedge_map;
    void after_edge_duplicated(Mesh::Halfedge_index h, Mesh::Halfedge_index new_hedge, const Mesh&)
    {
      hedge_map.emplace_back(h, new_hedge);
    }

    void after_vertex_copy(Mesh::Vertex_index v, const Mesh&, Mesh::Vertex_index nv, const Mesh&)
    {
      put(vnm, nv, get(vnm, v));
    }
  };
  Visitor visitor(vnm);

  PMP::internal::split_along_edges(mesh, ecm, mesh.points(), visitor);

  // TODO should actually only handle interior vertices...
  double delta = -0.005;
  for (const auto& ph : visitor.hedge_map)
  {
    Mesh::Halfedge_index h1 = ph.first;
    Mesh::Halfedge_index h2 = ph.second;
    if (is_border(h1, mesh)) h1=opposite(h1, mesh);
    if (is_border(h2, mesh)) h2=opposite(h2, mesh);
    Mesh::Halfedge_index h = in_out[face(h1, mesh)]==1 ? h1 : h2;

    Mesh::Vertex_index v = target(h, mesh);
    K::Vector_3 n = get(vnm, v);
    mesh.point(v) = mesh.point(v)+delta*n;
  }

  std::vector<Mesh::Halfedge_index> b1(visitor.hedge_map.size());
  std::vector<Mesh::Halfedge_index> b2(visitor.hedge_map.size());
  for (std::size_t i=0; i<visitor.hedge_map.size(); ++i)
  {
    Mesh::Halfedge_index h1 = visitor.hedge_map[i].first;
    Mesh::Halfedge_index h2 = visitor.hedge_map[i].second;
    if (is_border(h1, mesh)) h1=opposite(h1, mesh);
    if (is_border(h2, mesh)) h2=opposite(h2, mesh);
    if (in_out[face(h1, mesh)]==1) std::swap(h1,h2);

    b1[i]=opposite(h1, mesh);
    b2[i]=opposite(h2, mesh);

  }

  PMP::extrude_impl::create_strip(b1, b2, mesh);

  std::ofstream("mesh_refined_split.off") << std::setprecision(17) << mesh;

  return 0;
}
