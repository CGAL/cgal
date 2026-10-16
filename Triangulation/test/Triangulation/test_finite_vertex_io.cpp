#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/IO/io.h>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <vector>

template <typename Dimension>
void check_io(int dimension, int count, CGAL::IO::Mode mode)
{
  using Traits = CGAL::Epick_d<Dimension>;
  using Triangulation = CGAL::Triangulation<Traits>;
  using Point = typename Triangulation::Point;

  Triangulation original(dimension);
  std::vector<Point> points;
  std::vector<double> coordinates(dimension, 0);
  for(int i = 0; i < count; ++i)
  {
    std::fill(coordinates.begin(), coordinates.end(), 0);
    if(i > 0)
      coordinates[i - 1] = 1;
    points.emplace_back(dimension, coordinates.begin(), coordinates.end());
    original.insert(points.back());
  }
  assert(original.is_valid());
  std::fill(coordinates.begin(), coordinates.end(), 100);
  original.infinite_vertex()->set_point(Point(dimension, coordinates.begin(), coordinates.end()));

  std::stringstream output;
  CGAL::IO::set_mode(output, mode);
  output << original;

  if(count > 0)
  {
    std::istringstream vertices(output.str());
    CGAL::IO::set_mode(vertices, mode);
    int stored_dimension;
    std::size_t stored_count;
    if(mode == CGAL::IO::ASCII)
      vertices >> stored_dimension >> stored_count;
    else
    {
      CGAL::read(vertices, stored_dimension);
      CGAL::read(vertices, stored_count);
    }
    assert(stored_dimension == original.current_dimension());
    assert(stored_count == points.size());
    for(std::size_t i = 0; i < stored_count; ++i)
    {
      Point point;
      vertices >> point;
      assert(vertices);
      assert(std::find(points.begin(), points.end(), point) != points.end());
    }
  }

  Triangulation restored(dimension);
  output >> restored;
  assert(output);
  assert(restored.is_valid());
  assert(restored.current_dimension() == original.current_dimension());
  assert(restored.number_of_vertices() == original.number_of_vertices());
  assert(restored.number_of_full_cells() == original.number_of_full_cells());
  for(auto vertex = restored.finite_vertices_begin(); vertex != restored.finite_vertices_end(); ++vertex)
    assert(std::find(points.begin(), points.end(), vertex->point()) != points.end());
}

int main()
{
  for(const auto mode : {CGAL::IO::ASCII, CGAL::IO::BINARY})
    for(int count = 0; count <= 4; ++count)
    {
      check_io<CGAL::Dimension_tag<3>>(3, count, mode);
      check_io<CGAL::Dynamic_dimension_tag>(3, count, mode);
    }
}
