#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Kinetic_shape_reconstruction_2.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Segment_2 Segment_2;

typedef CGAL::Kinetic_shape_reconstruction_2<Kernel> Reconstruction;


int main (int argc, char** argv)
{
  CGAL::Random rand(0);
  std::vector<Segment_2> segments;
  
  for (std::size_t i = 0; i < 10; ++ i)
  {
    Point_2 source (rand.get_double(0, 5), rand.get_double(0, 5));
    Vector_2 vec (rand.get_double(-0.5, 0.5), rand.get_double(-0.5, 0.5));
    Point_2 target = source + vec;
    segments.push_back (Segment_2(source, target));
  }

  std::ofstream input_file ("input.polylines.txt");
  for (const Segment_2& s : segments)
    input_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;

  Reconstruction reconstruction;

  unsigned int k = 2;
  if (argc > 1)
    k = std::atoi(argv[1]);

  reconstruction.partition (segments, CGAL::Identity_property_map<Segment_2>(), k);

  segments.clear();
  
  reconstruction.output_partition_edges_to_segment_soup (std::back_inserter (segments));

  std::ofstream output_file ("output.polylines.txt");
  for (const Segment_2& s : segments)
    output_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;
  
  return EXIT_SUCCESS;
}
