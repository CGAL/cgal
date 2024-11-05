#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>

struct MinimalFrechetTraits {

  using Dimension = CGAL::Dimension_tag<2>;
  using FT = double;

  struct Point_d
  {
    Point_d(double, double) {}

    //~ double operator[](int) const
    //~ {
        //~ return 0;
    //~ }
  };

  struct Compute_squared_distance_d
  {
    double operator()(Point_d, Point_d) const
    {
      return 0;
    }
  };

  struct Construct_bbox_d
  {
    CGAL::Bbox_2 operator()(Point_d) const
    {
      return CGAL::Bbox_2();
    }
  };

  struct Cartesian_const_iterator_d
  {
    FT operator*() { return 0; }
    Cartesian_const_iterator_d& operator++() { return *this;}
    Cartesian_const_iterator_d operator++(int) { return Cartesian_const_iterator_d(); }
  };

  struct Construct_cartesian_const_iterator_d
  {
    Cartesian_const_iterator_d operator()(Point_d){ return Cartesian_const_iterator_d(); }
    Cartesian_const_iterator_d operator()(Point_d, int){ return Cartesian_const_iterator_d(); }
  };

};


#include <CGAL/Frechet_distance.h>

#include <vector>

int main()
{
  std::vector<MinimalFrechetTraits::Point_d> curve;
  /* bool decision = */ CGAL::is_Frechet_distance_larger<MinimalFrechetTraits>(curve, curve, 0.1);
  return 0;
}
