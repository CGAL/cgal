#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <array>
#include <iostream>
#include <optional>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

using Output_array = std::array<std::optional<Mesh*>, 4>;

enum Output_type { OUTPLACE, INPLACE_TM1, INPLACE_TM2 };

bool run_case(const Mesh& input, const Vector_3& translation,
              const std::size_t operation, const Output_type output_type)
{
  Mesh tm1 = input;
  Mesh tm2 = input;
  for(auto v : tm2.vertices())
    tm2.point(v) = tm2.point(v) + translation;

  Mesh out;
  Output_array output = {std::nullopt, std::nullopt, std::nullopt, std::nullopt};

  switch(output_type)
  {
    case OUTPLACE:       output[operation] = &out; break;
    case INPLACE_TM1:    output[operation] = &tm1; break;
    case INPLACE_TM2:    output[operation] = &tm2; break;
  }

  const auto success = PMP::corefine_and_compute_boolean_operations(
    tm1, tm2, output,
    CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

  Mesh* result = output[operation].value();
  return success[operation] && result != nullptr && !result->is_empty();
}

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_file>\n";
    return 1;
  }

  Mesh input;
  if(!CGAL::IO::read_polygon_mesh(argv[1], input) || input.is_empty())
    return 1;

  const CGAL::Bbox_3 bb = PMP::bbox(input);
  const Vector_3 translation((bb.xmax()-bb.xmin()) * 0.2,
                             (bb.ymax()-bb.ymin()) * 0.2,
                             (bb.zmax()-bb.zmin()) * 0.2);

  const char* output_names[] = {"outplace", "inplace_tm1", "inplace_tm2"};
  for(std::size_t op=0; op<4; ++op)
  {
    for(int type=OUTPLACE; type<=INPLACE_TM2; ++type)
    {
      const bool ok = run_case(input, translation, op, static_cast<Output_type>(type));
      std::cout << "output_" << op << "." << output_names[type]
                << ": " << (ok ? "OK" : "FAILED") << '\n';
      if(!ok)
        return 1;
    }
  }
  return 0;
}
