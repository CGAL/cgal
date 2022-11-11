#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>

#include "read_polygon.h"

#include <string.h>
#include <list>
#include <CGAL/Timer.h>

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
  VERTICAL_DECOMPOSITION,
  TRIANGULATION_DECOMPOSITION,
  VERTICAL_AND_ANGLE_BISECTOR_DECOMPOSITION,
  TRIANGULATION_AND_ANGLE_BISECTOR_DECOMPOSITION,
  OPTIMAL_DECOMPOSITION,
} Strategy;

static const char* strategy_names[] = {
  "reduced convolution",
  "vertical decomposition",
  "constrained triangulation decomposition",
  "vertical and angle bisector decomposition",
  "constrained triangulation and angle bisector decomposition",
  "optimal decomposition"
};

Polygon_with_holes_2 compute_minkowski_sum_2(Polygon_with_holes_2& p,
                                             Polygon_with_holes_2& q,
                                             Strategy strategy)
{
  switch (strategy) {
   case REDUCED_CONVOLUTION:
     return CGAL::minkowski_sum_by_reduced_convolution_2(p, q);

   case VERTICAL_DECOMPOSITION:
    {
     CGAL::Polygon_vertical_decomposition_2<Kernel> decomp;
     return CGAL::minkowski_sum_2(p, q, decomp);
    }

   case TRIANGULATION_DECOMPOSITION:
    {
     CGAL::Polygon_triangulation_decomposition_2<Kernel> decomp;
     return CGAL::minkowski_sum_2(p, q, decomp);
    }

   case VERTICAL_AND_ANGLE_BISECTOR_DECOMPOSITION:
    {
     typedef CGAL::Small_side_angle_bisector_decomposition_2<Kernel>
                                                No_holes_decomposition;
     typedef CGAL::Polygon_vertical_decomposition_2<Kernel>
                                                With_holes_decomposition;

     if (0 == p.number_of_holes()) {
       const Polygon_2& pnh = p.outer_boundary();
       No_holes_decomposition decomp_no_holes;
       if  (0 == q.number_of_holes()) {
         const Polygon_2& qnh = q.outer_boundary();
         return CGAL::minkowski_sum_2(pnh, qnh, decomp_no_holes, decomp_no_holes);
       }

       With_holes_decomposition decomp_with_holes;
       return CGAL::minkowski_sum_2(pnh, q, decomp_no_holes, decomp_with_holes);
     }

     With_holes_decomposition decomp_with_holes;
     if (0 == q.number_of_holes()) {
       const Polygon_2& qnh = q.outer_boundary();
       No_holes_decomposition decomp_no_holes;
       return CGAL::minkowski_sum_2(p, qnh, decomp_with_holes, decomp_no_holes);
     }

     return CGAL::minkowski_sum_2(p, q, decomp_with_holes, decomp_with_holes);
    }

   case TRIANGULATION_AND_ANGLE_BISECTOR_DECOMPOSITION:
    {
     typedef CGAL::Small_side_angle_bisector_decomposition_2<Kernel>
                                                No_holes_decomposition;
     typedef CGAL::Polygon_triangulation_decomposition_2<Kernel>
                                                With_holes_decomposition;
     if (0 == p.number_of_holes()) {
       const Polygon_2& pnh = p.outer_boundary();
       No_holes_decomposition decomp_no_holes;
       if (0 == q.number_of_holes()) {
         const Polygon_2& qnh = q.outer_boundary();
         return CGAL::minkowski_sum_2(pnh, qnh, decomp_no_holes, decomp_no_holes);
       }

       With_holes_decomposition decomp_with_holes;
       return CGAL::minkowski_sum_2(pnh, q, decomp_no_holes, decomp_with_holes);
     }

     With_holes_decomposition decomp_with_holes;
     if (0 == q.number_of_holes()) {
       const Polygon_2& qnh = q.outer_boundary();
       No_holes_decomposition decomp_no_holes;
       return CGAL::minkowski_sum_2(p, qnh, decomp_with_holes, decomp_no_holes);
     }

     return CGAL::minkowski_sum_2(p, q, decomp_with_holes, decomp_with_holes);
    }

   case OPTIMAL_DECOMPOSITION:
    {
     CGAL::Small_side_angle_bisector_decomposition_2<Kernel> decomp_no_holes;
     CGAL::Polygon_triangulation_decomposition_2<Kernel> decomp_with_holes;
     return CGAL::minkowski_sum_by_decomposition_2(p, q,
                                                   decomp_no_holes,
                                                   decomp_with_holes);
    }

   default:
    std::cerr << "Invalid strategy" << std::endl;
    return Polygon_with_holes_2();
  }
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [-method flag] [polygon files]..."
              << std::endl;
    std::cerr << "For the method flag, use a subset of the letters 'rfsohg'."
              << std::endl;
    std::cerr << "The program will compute the Minkowski sum of the first "
              << "and second polygon, of the third and fourth, and so on."
              << std::endl;
    return 1;
  }

  Polygon_with_holes_2 p, q;
  CGAL::Timer timer;
  timer.start();

  std::list<Strategy> strategies;

  int i(1);
  while (i < argc) {
    if (argv[i][0] == '-') {
      strategies.clear();
      for (std::size_t j = 1; j < strlen(argv[i]); ++j) {
        switch (argv[i][j]) {
         case 'r': strategies.push_back(REDUCED_CONVOLUTION); break;
         case 'v': strategies.push_back(VERTICAL_DECOMPOSITION); break;
         case 't': strategies.push_back(TRIANGULATION_DECOMPOSITION); break;
         case 'w':
          strategies.push_back(VERTICAL_AND_ANGLE_BISECTOR_DECOMPOSITION);
          break;

         case 'u':
          strategies.push_back(TRIANGULATION_AND_ANGLE_BISECTOR_DECOMPOSITION);
          break;

         case 'd': strategies.push_back(OPTIMAL_DECOMPOSITION); break;
         default:
          std::cerr << "Unknown flag '" << argv[i][j] << "'" << std::endl;
          return -1;
        }
      }
      ++i;
      continue;
    }

    std::cout << "Testing " << argv[i] << " + " << argv[i+1] << std::endl;
    if (!read_polygon(argv[i], p)) return -1;
    if (!read_polygon(argv[i+1], q)) return -1;

    bool compare = false;
    Polygon_with_holes_2 reference;
    std::list<Strategy>::iterator it;
    for (it = strategies.begin(); it != strategies.end(); ++it) {
      std::cout << "Using " << strategy_names[*it] << ": ";
      timer.reset();
      Polygon_with_holes_2 result = compute_minkowski_sum_2(p, q, *it);
      double secs = timer.time();
      std::cout << secs << " s " << std::flush;

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
