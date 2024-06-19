#include <CGAL/tags.h>

struct MinimalFrechetTraits {

  static const int dimension = 2;

  struct Kernel {
    typedef double FT;
    enum { Has_filtered_predicates = false };
  };

  struct Point {
    double operator[](int) const
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
  bool decision = CGAL::Frechet_distance_at_most<MinimalFrechetTraits>(curve, curve, 0.1);
  return 0;
}
