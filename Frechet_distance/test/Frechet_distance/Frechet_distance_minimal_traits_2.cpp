#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>

template <class FT_>
struct MinimalFrechetTraits {
  using Dimension = CGAL::Dimension_tag<2>;
  using FT = FT_;

  struct Point_d
  {};

  struct Construct_bbox_d
  {
    Construct_bbox_d() = delete;
    Construct_bbox_d(int){}
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
    Construct_cartesian_const_iterator_d() = delete;
    Construct_cartesian_const_iterator_d(int){}

    Cartesian_const_iterator_d operator()(Point_d){ return Cartesian_const_iterator_d(); }
    Cartesian_const_iterator_d operator()(Point_d, int){ return Cartesian_const_iterator_d(); }
  };

  Construct_bbox_d construct_bbox_d_object() const
  {
    return Construct_bbox_d(0);
  }

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const
  {
    return Construct_cartesian_const_iterator_d(0);
  }
};


#include <CGAL/Frechet_distance.h>

#include <vector>

int main()
{
  namespace params = CGAL::parameters;

  {
    using Traits = MinimalFrechetTraits<double>;
    std::vector<Traits::Point_d> curve;
    /* bool decision = */ CGAL::is_Frechet_distance_larger(
        curve, curve, 0.1, params::force_filtering(std::true_type())
                                  .geom_traits(Traits()));
  }

  {
    using Traits = MinimalFrechetTraits<CGAL::Exact_rational>;
    std::vector<Traits::Point_d> curve;
    /* bool decision = */ CGAL::is_Frechet_distance_larger(
        curve, curve, 0.1, params::force_filtering(std::true_type())
                                  .geom_traits(Traits()));
  }

  return 0;
}
