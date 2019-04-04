#include <fstream>

//#define CGAL_KSR_VERBOSE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Kinetic_shape_reconstruction_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Segment_2 Segment_2;

typedef CGAL::Surface_mesh<Point_2> Mesh;

typedef CGAL::Kinetic_shape_reconstruction_2<Kernel> Reconstruction;


int main (int argc, char** argv)
{
  CGAL::Random rand(0);
  std::vector<Segment_2> segments;

  unsigned int nb_lines = 30;
  if (argc > 1)
    nb_lines = std::atoi(argv[1]);
  unsigned int k = 2;
  if (argc > 2)
    k = std::atoi(argv[2]);
#define REGULAR_CASE
#ifdef REGULAR_CASE
  segments.push_back (Segment_2(Point_2 (0, 1), Point_2 (0, 3)));
  segments.push_back (Segment_2(Point_2 (0, 5), Point_2 (0, 7)));
  segments.push_back (Segment_2(Point_2 (4, 1), Point_2 (4, 3)));
  segments.push_back (Segment_2(Point_2 (4, 6), Point_2 (4, 7)));
  segments.push_back (Segment_2(Point_2 (1, 0), Point_2 (3, 0)));
  segments.push_back (Segment_2(Point_2 (2, 4), Point_2 (3, 4)));
  segments.push_back (Segment_2(Point_2 (1.2, 8), Point_2 (2.5, 8)));
#else
  
  for (unsigned int i = 0; i < nb_lines; ++ i)
  {
    Point_2 source (rand.get_double(0, 5), rand.get_double(0, 5));
    Vector_2 vec (rand.get_double(-0.5, 0.5), rand.get_double(-0.5, 0.5));
    Point_2 target = source + vec;
    segments.push_back (Segment_2(source, target));
  }
#endif
  
  std::ofstream input_file ("input.polylines.txt");
  for (const Segment_2& s : segments)
    input_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;

  Reconstruction reconstruction;


  reconstruction.partition (segments, CGAL::Identity_property_map<Segment_2>(), k);

  segments.clear();
  
  reconstruction.output_partition_edges_to_segment_soup (std::back_inserter (segments));

  std::ofstream output_file ("output.polylines.txt");
  for (const Segment_2& s : segments)
    output_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;

  if (!reconstruction.check_integrity(true))
  {
    std::cerr << "Integrity of reconstruction failed" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh mesh;

  if (reconstruction.output_partition_cells_to_face_graph(mesh))
  {
    std::cerr << mesh.number_of_vertices() << " vertices and " << mesh.number_of_faces() << " faces" << std::endl;
    
    std::ofstream output_shapes_file ("out.ply");
    output_shapes_file << "ply" << std::endl
                       << "format ascii 1.0" << std::endl
                       << "element vertex " << mesh.number_of_vertices() << std::endl
                       << "property double x" << std::endl
                       << "property double y" << std::endl
                       << "property double z" << std::endl
                       << "element face " << mesh.number_of_faces() << std::endl
                       << "property list uchar int vertex_index" << std::endl
                       << "property uchar red" << std::endl
                       << "property uchar green" << std::endl
                       << "property uchar blue" << std::endl
                       << "end_header" << std::endl;
    for (const auto& vindex : vertices(mesh))
      output_shapes_file << mesh.point(vindex) << " 0" << std::endl;
    for (const auto& findex : faces(mesh))
    {
      output_shapes_file << degree(findex, mesh);
      for (const auto& hindex : CGAL::halfedges_around_face(halfedge(findex,mesh),mesh))
        output_shapes_file << " " << int(target(hindex,mesh));
      output_shapes_file << " " << rand.get_int(64,192)
                         << " " << rand.get_int(64,192)
                         << " " << rand.get_int(64,192) << std::endl;
    }
  }
  else
    std::cerr << "Invalid face graph" << std::endl;
  
  return EXIT_SUCCESS;
}
