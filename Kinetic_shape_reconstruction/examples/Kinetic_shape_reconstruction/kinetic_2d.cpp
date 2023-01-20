#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_2.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/Random.h>

#include <vector>
#include <fstream>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point_2   = typename Kernel::Point_2;
using Point_3   = typename Kernel::Point_3;
using Vector_2  = typename Kernel::Vector_2;
using Segment_2 = typename Kernel::Segment_2;

using Transform    = CGAL::Aff_transformation_2<Kernel>;
using Surface_mesh = CGAL::Surface_mesh<Point_2>;
using KSR          = CGAL::Kinetic_shape_reconstruction_2<Kernel>;

void add_regular_case(std::vector<Segment_2>& segments, CGAL::Random& rand) {

  const std::size_t size_before = segments.size();
  segments.push_back(Segment_2(Point_2(0.0, 1.0), Point_2(0.0, 3.0)));
  segments.push_back(Segment_2(Point_2(0.0, 5.0), Point_2(0.0, 7.0)));
  segments.push_back(Segment_2(Point_2(4.0, 1.0), Point_2(4.0, 3.0)));
  segments.push_back(Segment_2(Point_2(4.0, 6.0), Point_2(4.0, 7.0)));
  segments.push_back(Segment_2(Point_2(1.0, 0.0), Point_2(3.0, 0.0)));
  segments.push_back(Segment_2(Point_2(2.0, 4.0), Point_2(3.0, 4.0)));
  segments.push_back(Segment_2(Point_2(1.2, 8.0), Point_2(2.5, 8.0)));

  // Random rotation.
  const double sine = rand.get_double(-1.1);
  const double cosine = CGAL::sqrt(1.0 - sine * sine);
  const Transform rotate(CGAL::Rotation(), sine, cosine);
  const Transform scale(CGAL::Scaling(), rand.get_double(0.1, 10));
  const Transform translate(CGAL::Translation(),
    Vector_2(rand.get_double(-5, 5), rand.get_double(-5, 5)));
  const Transform transform = scale * rotate * translate;

  for (std::size_t i = size_before; i < segments.size(); ++i) {
    const Point_2 source = transform.transform(segments[i].source());
    const Point_2 target = transform.transform(segments[i].target());
    segments[i] = Segment_2(source, target);
  }
}

int main(int argc, char** argv) {

  CGAL::Random rand(0);
  std::vector<Segment_2> segments;
  #define REGULAR_CASE

  unsigned int nb_lines = 30;
  if (argc > 1) {
    nb_lines = std::atoi(argv[1]);
  }
  unsigned int k = 2;
  if (argc > 2) {
    k = std::atoi(argv[2]);
  }

  #ifdef REGULAR_CASE
    add_regular_case(segments, rand);
  #else
    for (unsigned int i = 0; i < nb_lines; ++i) {
      const Point_2 source(rand.get_double(0, 5), rand.get_double(0, 5));
      const Vector_2 vec(rand.get_double(-0.5, 0.5), rand.get_double(-0.5, 0.5));
      const Point_2 target = source + vec;
      segments.push_back(Segment_2(source, target));
    }
  #endif

  std::ofstream input_file("input.polylines.txt");
  for (const Segment_2& segment : segments) {
    input_file << "2 " << segment.source() << " 0 " << segment.target() << " 0" << std::endl;
  }

  KSR ksr;
  ksr.partition(segments, CGAL::Identity_property_map<Segment_2>(), k, 2);

  segments.clear();
  ksr.output_raw_partition_edges_to_segment_soup(std::back_inserter(segments));
  std::ofstream raw_output_file("output_raw.polylines.txt");
  for (const Segment_2& segment : segments) {
    raw_output_file << "2 " << segment.source() << " 0 " << segment.target() << " 0" << std::endl;
  }

  segments.clear();
  ksr.output_partition_edges_to_segment_soup(std::back_inserter(segments));
  std::ofstream output_file("output.polylines.txt");
  for (const Segment_2& segment : segments) {
    output_file << "2 " << segment.source() << " 0 " << segment.target() << " 0" << std::endl;
  }

  if (!ksr.check_integrity(true)) {
    std::cerr << "ERROR: KSR INTEGRITY FAILED!" << std::endl;
    return EXIT_FAILURE;
  }

  Surface_mesh mesh;
  if (ksr.output_partition_cells_to_face_graph(mesh)) {
    std::cout << mesh.number_of_vertices() <<
    " vertices and " << mesh.number_of_faces() << " faces" << std::endl;

    std::ofstream output_shapes_file("ksr.ply");
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

    for (const auto& vindex : vertices(mesh)) {
      output_shapes_file << mesh.point(vindex) << " 0" << std::endl;
    }

    for (const auto& findex : faces(mesh)) {
      output_shapes_file << degree(findex, mesh);
      for (const auto& hindex : CGAL::halfedges_around_face(halfedge(findex,mesh), mesh)) {
        output_shapes_file << " " << int(target(hindex,mesh));
      }
      output_shapes_file
        << " " << rand.get_int(64,192)
        << " " << rand.get_int(64,192)
        << " " << rand.get_int(64,192) << std::endl;
    }
  } else {
    std::cerr << "ERROR: INVALID FACE GRAPH!" << std::endl;
  }
  return EXIT_SUCCESS;
}
