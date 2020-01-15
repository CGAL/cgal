// #define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
// #define CGAL_PMP_SNAP_DEBUG_PP
// #define CGAL_PMP_SNAP_DEBUG_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/PLY_reader.h>

#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <boost/property_map/function_property_map.hpp>

#include <iostream>
#include <fstream>

//typedef CGAL::Simple_cartesian<double>                                Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel           Kernel;

typedef Kernel::FT                                                    FT;
typedef Kernel::Point_3                                               Point_3;
typedef Kernel::Vector_3                                              Vector_3;
typedef CGAL::Surface_mesh<Point_3>                                   Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor            face_descriptor;

typedef CGAL::dynamic_face_property_t<std::size_t>                    Face_patch_id_tag;
typedef boost::property_map<Surface_mesh, Face_patch_id_tag>::type    Face_patch_map;

namespace PMP = CGAL::Polygon_mesh_processing;

#define MARK_COMPATIBLE_FACES
#define CLEAN_INPUT_MESH
#define REMOVE_DEGENERATE_FACES

const double snapping_tolerance = 0.001;
const double squared_snapping_tolerance = CGAL::square(snapping_tolerance);

template <typename K, typename Mesh>
bool read_mesh(const char* filename, Mesh& sm)
{
  typedef typename K::Point_3                                   Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::vector<int> > faces;
    if(!CGAL::read_STL(in, points, faces))
    {
      std::cerr << "Error: cannot read STL mesh\n";
      return false;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "off")
  {
    if(!(in >> sm))
    {
      std::cerr << "Cannot read .OFF as polygon mesh, reading .OFF as a polygon soup...\n";

      in.clear(); // clear fail and eof bits
      in.seekg(0, std::ios::beg); // back at the start

      std::vector<Point> points;
      std::vector<std::vector<int> > faces;
      if(!CGAL::read_OFF(in, points, faces))
      {
        std::cerr << "Error: cannot read OFF mesh as a polygon soup\n";
        return false;
      }

      std::cout << "Cleaning polygon soup..." << std::endl;
      PMP::repair_polygon_soup(points, faces);

      if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
        std::cerr << "W: File does not describe a polygon mesh" << std::endl;

      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
    }
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "ply")
  {
    std::vector<Point> points;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;

    if(!(CGAL::read_PLY(in, points, polygons, fcolors, vcolors)))
    {
      std::cerr << "Error: cannot read PLY mesh\n";
      return false;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, polygons);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return false;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " nv "
                              << num_edges(sm) << " ne "
                              << num_faces(sm) << " nf" << std::endl;

  if(!CGAL::is_triangle_mesh(sm))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return false;
  }

  return true;
}

// The goal is to remove needle-like very thin triangles that do not really contribute to the result,
// but make life a lot harder
void clean_mesh(Surface_mesh& sm)
{
  std::size_t nb_removed_CC = PMP::remove_connected_components_of_negligible_size(
                                sm, CGAL::parameters::area_threshold(0.005));
  std::cout << "Removed " << nb_removed_CC << " tiny CC(s)" << std::endl;

  int re = 0;
  std::set<face_descriptor> faces_to_remove;

  for(edge_descriptor e : edges(sm))
  {
    if(!is_border(e, sm))
      continue;

    if(CGAL::squared_distance(sm.point(source(e, sm)), sm.point(target(e, sm))) < squared_snapping_tolerance)
    {
      ++re;
      if(CGAL::Euler::does_satisfy_link_condition(e, sm))
      {
        CGAL::Euler::collapse_edge(e, sm);
      }
      else
      {
        halfedge_descriptor h = halfedge(e, sm);
        if(!is_border(h, sm))
          h = opposite(h, sm);

        if(PMP::internal::border_size(h, sm) == 3) // lone face
          faces_to_remove.insert(face(opposite(h, sm), sm)); // can't remove yet since we're looping edges
        else
          std::cerr << "Can't collapse " << e << " ?" << std::endl;
      }
    }
  }

  for(face_descriptor f : faces_to_remove)
    CGAL::Euler::remove_face(halfedge(f, sm), sm);

  std::cout << "Removed " << faces_to_remove.size() << " tiny faces" << std::endl;
  std::cout << "Removed " << re << " tiny edges" << std::endl;

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  std::ofstream("results/cleaned.off") << std::setprecision(17) << sm;
#endif
}

template <typename FaceContainer, typename TriangleMesh>
void dump_cc(const FaceContainer& cc_faces,
             const TriangleMesh& mesh,
             const std::string filename)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  std::ofstream out(filename);
  out.precision(17);

  out << "OFF\n";
  out << 3*cc_faces.size() << " " << cc_faces.size() << " 0\n";

  for(const face_descriptor f : cc_faces)
  {
    out << mesh.point(source(halfedge(f, mesh), mesh)) << "\n";
    out << mesh.point(target(halfedge(f, mesh), mesh)) << "\n";
    out << mesh.point(target(next(halfedge(f, mesh), mesh), mesh)) << "\n";
  }

  int id = 0;
  for(const face_descriptor f : cc_faces)
  {
    CGAL_USE(f);
    out << "3 " << id << " " << id+1 << " " << id+2 << "\n";
    id += 3;
  }

  out.close();
}

// Quite naive currently as it simply groups together faces that have the same the normal
// and assumes that patches that will be snapped together are all planar
std::size_t mark_compatible_patches(Face_patch_map fpmap,
                                    const Surface_mesh& sm)
{
  std::vector<Vector_3> distinct_patch_normals;

  for(face_descriptor f : faces(sm))
    put(fpmap, f, std::size_t(-1));

  std::size_t new_id = 0;
  for(face_descriptor f : faces(sm))
  {
    if(get(fpmap, f) != std::size_t(-1)) // already got a patch id, nothing to do
      continue;

    Vector_3 fn = PMP::compute_face_normal(f, sm);
    const FT sq_fn_l = fn.squared_length();

    std::size_t id = new_id;
    for(std::size_t i=0, s=distinct_patch_normals.size(); i<s; ++i)
    {
      const Vector_3& n = distinct_patch_normals[i];
      const FT sq_n_l = n.squared_length();
      const FT sp = CGAL::scalar_product(fn, n);

      // less than 10 degrees difference ==> not distinct ==> same patch
      const FT sq_cos = 0.99922859349937227; // CGAL::square(std::cos(10/360));

      if(CGAL::square(sp) >= sq_fn_l * sq_n_l * sq_cos)
      {
        id = i;
        break;
      }
    }

    if(id == new_id)
    {
      distinct_patch_normals.push_back(fn);
      ++new_id;
    }

    put(fpmap, f, id);
  }

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
    for(std::size_t i=0; i<new_id; ++i)
    {
      std::vector<face_descriptor> cc_fs;
      for(face_descriptor f : faces(sm))
        if(get(fpmap, f) == i)
          cc_fs.push_back(f);

      std::stringstream oss;
      oss << "results/cc_" << i << ".off" << std::ends;
      dump_cc(cc_fs, sm, oss.str().c_str());
    }
#endif

  return new_id;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);

  Surface_mesh sm;
  if(argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " input_mesh" << std::endl;
    return EXIT_FAILURE;
  }

  const char* filename = argv[1];
  if(!read_mesh<Kernel>(filename, sm))
    return EXIT_FAILURE;

  Surface_mesh::Property_map<vertex_descriptor, double> tolerance_map;
  tolerance_map = sm.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(sm))
    put(tolerance_map, v, snapping_tolerance);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // Remove small edges, small faces, and small CCs
#ifdef CLEAN_INPUT_MESH
  clean_mesh(sm);
#endif

  std::chrono::steady_clock::time_point clean_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (clean): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(clean_time - start_time).count()
            << "ms" << std::endl;

  // Color the faces of the mesh into compatible components
#ifdef MARK_COMPATIBLE_FACES
  Face_patch_map fpmap = get(Face_patch_id_tag(), sm);
  std::size_t patch_n = mark_compatible_patches(fpmap, sm);
  std::cout << patch_n << " different patch(es)" << std::endl;
#endif

  std::chrono::steady_clock::time_point segmentation_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (segmentation): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(segmentation_time - clean_time).count()
            << "ms" << std::endl;

  // Snap
  std::size_t nb_snapped = PMP::experimental::snap_borders<CGAL::Parallel_tag>(
                             sm, tolerance_map, CGAL::parameters::face_patch_map(fpmap));
  std::cout << "#snapped: " << nb_snapped << std::endl;

  std::chrono::steady_clock::time_point snap_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (snap): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(snap_time - segmentation_time).count()
            << "ms" << std::endl;

#ifdef REMOVE_DEGENERATE_FACES
  PMP::remove_degenerate_faces(sm);
#endif

  std::chrono::steady_clock::time_point degen_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (remove degen faces): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(snap_time - degen_time).count()
            << "ms" << std::endl;

  // Stitch
  std::cout << "Stitch, #ne: " << edges(sm).size() << std::endl;
  PMP::stitch_borders(sm);

  std::chrono::steady_clock::time_point stitch_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (stitch): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stitch_time - snap_time).count()
            << "ms" << std::endl;

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  std::cout << "Total time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "#border: " << PMP::number_of_borders(sm) << std::endl;
  std::cout << "Done!" << std::endl;

  std::ofstream("results/snapped.off") << std::setprecision(17) << sm;

  return EXIT_SUCCESS;
}
