#include "split_polylines.h"

#include <CGAL/Mesh_3/polylines_to_protect.h>

auto split_polylines(const Polylines_container& input)
  -> std::vector<Polylines_container::value_type>
{
  std::vector<Polylines_container::value_type> result;
  CGAL::polylines_to_protect(result, input.begin(), input.end());
  return result;
}
