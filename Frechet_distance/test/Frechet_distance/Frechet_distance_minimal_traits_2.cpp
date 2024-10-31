#include <CGAL/Dimension.h>

struct MinimalFrechetTraits {

  using Dimension = CGAL::Dimension_tag<2>;
  using FT = double;

  struct Point_d {
    Point_d(double, double) {}

    double operator[](int) const
    {
        return 0;
    }
  };

  struct Compute_squared_distance_d
  {
    double operator()(Point, Point) const
    {
      return 0;
    }
  };

  struct Construct_bbox_d
  {
    CGAL::Bbox_2 operator()(Point_d) const
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
