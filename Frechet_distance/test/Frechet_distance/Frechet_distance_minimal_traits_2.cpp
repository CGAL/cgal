#include <CGAL/Dimension.h>

struct MinimalFrechetTraits {

  using Dimension = CGAL::Dimension_tag<2>;
  using FT = double;

  struct Point {
    Point(double, double) {}

    double operator[](int) const
    {
        return 0;
    }
  };

  struct Squared_distance
  {
    double operator()(Point, Point) const
    {
      return 0;
    }
  };

};


#include <CGAL/Frechet_distance.h>

#include <vector>

int main()
{
  std::vector<MinimalFrechetTraits::Point> curve;
  /* bool decision = */ CGAL::is_Frechet_distance_larger<MinimalFrechetTraits>(curve, curve, 0.1);
  return 0;
}
