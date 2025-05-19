#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Timer.h>

#include "read_polygon.h"

#include <string.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

bool are_equal(const Polygon_with_holes_2& ph1,
               const Polygon_with_holes_2& ph2)
{
  std::list<Polygon_with_holes_2> sym_diff;
  CGAL::symmetric_difference(ph1, ph2, std::back_inserter(sym_diff));
  return sym_diff.empty();
}

typedef enum {
  REDUCED_CONVOLUTION,
  FULL_CONVOLUTION,
  SSAB_DECOMP,
  OPT_DECOMP,
  HM_DECOMP,
  GREENE_DECOMP,
  TRIANGULATION_DECOMP,
  VERTICAL_DECOMP
} Strategy;

static const char* strategy_names[] = {
  "reduced convolution",
  "full convolution",
  "small-side angle-bisector decomposition",
  "optimal convex decomposition",
  "Hertel-Mehlhorn decomposition",
  "Greene decomposition",
  "Triangulation",
  "Vertical decomposition"
};

Polygon_with_holes_2 compute_minkowski_sum_2(Polygon_2& p, Polygon_2& q,
                                             Strategy strategy)
{
  switch (strategy) {
   case REDUCED_CONVOLUTION:
    return minkowski_sum_by_reduced_convolution_2(p, q);

   case FULL_CONVOLUTION:
    return minkowski_sum_by_full_convolution_2(p, q);

   case SSAB_DECOMP:
    {
     CGAL::Small_side_angle_bisector_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }

   case OPT_DECOMP:
    {
     CGAL::Optimal_convex_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }

   case HM_DECOMP:
    {
     CGAL::Hertel_Mehlhorn_convex_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }

   case GREENE_DECOMP:
    {
     CGAL::Greene_convex_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }

   case TRIANGULATION_DECOMP:
    {
     CGAL::Polygon_triangulation_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }

   default: // VERTICAL_DECOMP
    {
     CGAL::Polygon_vertical_decomposition_2<Kernel> decomp;
     return minkowski_sum_2(p, q, decomp);
    }
  }
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [method flag] [polygon files]..."
              << std::endl;
    std::cerr << "For the method flag, use a subset of the letters 'rfsohg'."
              << std::endl;
    std::cerr << "The program will compute the Minkowski sum of the first "
              << "and second polygon, of the third and fourth, and so on."
              << std::endl;
    return 1;
  }

  Polygon_2 p, q;
  CGAL::Timer timer;

  std::list<Strategy> strategies;
  for (std::size_t i = 0; i < strlen(argv[1]); ++i) {
    switch (argv[1][i]) {
      case 'r':
        strategies.push_back(REDUCED_CONVOLUTION);
        break;
      case 'f':
        strategies.push_back(FULL_CONVOLUTION);
        break;
      case 's':
        strategies.push_back(SSAB_DECOMP);
        break;
      case 'o':
        strategies.push_back(OPT_DECOMP);
        break;
      case 'h':
        strategies.push_back(HM_DECOMP);
        break;
      case 'g':
        strategies.push_back(GREENE_DECOMP);
        break;
      case 't':
        strategies.push_back(TRIANGULATION_DECOMP);
        break;
      case 'v':
        strategies.push_back(VERTICAL_DECOMP);
        break;
      default:
        std::cerr << "Unknown flag '" << argv[1][i] << "'" << std::endl;
        return 1;
    }
  }

  int i = 2;
  while (i+1 < argc) {
    std::cout << "Testing " << argv[i] << " + " << argv[i+1] << std::endl;
    if (! read_polygon(argv[i], p)) return -1;
    if (! read_polygon(argv[i+1], q)) return -1;

    bool compare = false;
    Polygon_with_holes_2 reference;
    std::list<Strategy>::iterator it;
    for (it = strategies.begin(); it != strategies.end(); ++it) {
      std::cout << "Using " << strategy_names[*it] << ": ";
      timer.reset();
      timer.start();
      Polygon_with_holes_2 result = compute_minkowski_sum_2(p, q, *it);
      timer.stop();
      std::cout << timer.time() << " s " << std::flush;

      if (compare) {
        if (are_equal(reference, result)) std::cout << "(OK)";
        else {
          std::cout << "(ERROR: different result)";
          return 1;
        }
      }
      else {
        compare = true;
        reference = result;

        std::size_t n = result.outer_boundary().size();
        Polygon_with_holes_2::Hole_const_iterator it = result.holes_begin();
        while(it != result.holes_end()) n += (*it++).size();

        std::cout << std::endl << "Result has " << n << " vertices and "
                  << result.number_of_holes() << " holes.";
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;
    i += 2;
  }

  return 0;
}
