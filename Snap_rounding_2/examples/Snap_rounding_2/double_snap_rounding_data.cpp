
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/snap_rounding_2.h>
#include <CGAL/IO/OBJ.h>

#include <CGAL/Real_timer.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;

using FT = K::FT;
using Point_2 = K::Point_2;
using Segment_2 = K::Segment_2;

struct Obj_data {
  std::vector<Point_2> points;
  std::vector<std::vector<std::size_t> > polylines;
};

static int run_mesh(const Obj_data& data)
{
  std::vector< Segment_2 > segs;
  std::vector< Segment_2 > out;

  for(const std::vector<std::size_t>& id_pl : data.polylines)
    for(std::size_t pid = 1; pid < id_pl.size(); ++pid)
      segs.emplace_back(data.points[id_pl[pid - 1]], data.points[id_pl[pid]]);

  CGAL::Real_timer t;
  t.start();
  CGAL::snap_rounding_2(segs, std::back_inserter(out));
  t.stop();
  std::cout << "Running time: " << t.time() << std::endl;
  std::cout << "Input_size: "<< segs.size() << " , Output size: " << out.size() << "s" << std::endl;
  std::cout << "Do output intersect: " << CGAL::do_curves_intersect(out.begin(), out.end()) << "\n" << std::endl;
  return EXIT_SUCCESS;
}

int main(int argc, char*argv[] )
{
  std::cerr.precision(17);
  std::cout.precision(17);
  std::clog.precision(17);
  std::string path = argc > 1 ? argv[1] : CGAL::data_file_path("2d_segments/mini.obj");
  std::ifstream in(path);
  if(!in){
    std::cout << "File not found: " << path << std::endl;
    return EXIT_FAILURE;
  }

  Obj_data data;
  std::vector<std::vector<std::size_t> > unused_id_polygons;
  bool success = CGAL::IO::internal::read_OBJ(in, data.points, data.polylines, unused_id_polygons);
  if(!success)
    return EXIT_FAILURE;

  std::cout << path << std::endl;
  return run_mesh(data);
}
