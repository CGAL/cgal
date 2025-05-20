#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

bool treat_dataset(const std::string& dataset)
{
  std::string base = "";

  Mesh sm;
  auto fcm = sm.add_property_map<face_descriptor, CGAL::IO::Color>("f:color").first;
  auto fwm = sm.add_property_map<face_descriptor, double>("f:weight").first;

  std::size_t current_color_id = 0;
  std::vector<CGAL::IO::Color> palette {{ CGAL::IO::red(),
                                          CGAL::IO::green(),
                                          CGAL::IO::blue(),
                                          CGAL::IO::purple(),
                                          CGAL::IO::orange(),
                                          CGAL::IO::deep_blue(),
                                          CGAL::IO::yellow(),
                                          CGAL::IO::black(),
                                          CGAL::IO::violet(),
                                          CGAL::IO::gray(),
                                          CGAL::IO::white() }};
  std::map<double, CGAL::IO::Color> weight_colors;

  std::size_t weight_id = 0;
  for(;;)
  {
    std::string weight_filename = base + "offset_" + std::to_string(weight_id) + ".txt";
    std::string ply_filename = base + "face_" + std::to_string(weight_id) + ".ply";

    std::ifstream in(weight_filename);
    if(!in)
    {
      std::cerr << "No offset file at " << weight_filename << " -- stopping" << std::endl;
      break;
    }

    double weight_value;
    if(!(in >> weight_value))
    {
      std::cerr << "Error: failed to read offset for face #" << weight_id << std::endl;
      return false;
    }

    CGAL::IO::Color color = palette.at(current_color_id);
    auto res = weight_colors.emplace(weight_value, color);
    if(res.second) // never seen that weight before
      current_color_id = (current_color_id + 1) % palette.size();
    color = res.first->second;

    Mesh local_sm;
    if(!CGAL::IO::read_PLY(ply_filename, local_sm))
    {
      std::cerr << "Failed to read " << ply_filename << std::endl;
      return false;
    }

    std::cout << "Local SM with " << num_faces(local_sm) << " faces" << std::endl;

    std::unordered_map<face_descriptor, face_descriptor> f2f;
    CGAL::copy_face_graph(local_sm, sm, CGAL::parameters::face_to_face_output_iterator(std::inserter(f2f, f2f.end())));
    for(const auto& e : f2f)
    {
      put(fwm, e.second, weight_value);
      put(fcm, e.second, color);
    }

    ++weight_id;
  }

  std::ofstream out("colored_input_0.ply");
  if(!out)
  {
    std::cerr << "Error: failed to open out file" << std::endl;
    return false;
  }

  CGAL::IO::write_PLY(out, sm, CGAL::parameters::use_binary_mode(false)
                                                .stream_precision(17)
                                                .face_color_map(fcm));

  return true;
}

int main(int, char**)
{
  treat_dataset("tmp");
  treat_dataset("acute");
  treat_dataset("merging");
  treat_dataset("performance");

  return EXIT_SUCCESS;
}
