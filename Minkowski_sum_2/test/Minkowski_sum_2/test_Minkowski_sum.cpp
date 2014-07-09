#include <CGAL/Arithmetic_kernel.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Arithmetic_kernel::Rational         Rational;

#include <CGAL/Cartesian.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/Polygon_set_2.h>

#include "read_polygon.h"
#include <cstring>

#include <list>

typedef CGAL::Cartesian<Rational>                   Kernel;
typedef Kernel::Point_2                             Point_2;
typedef Kernel::Segment_2                           Segment_2;
typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

typedef CGAL::Polygon_set_2<Kernel>                 Polygon_set_2;
typedef Polygon_set_2::Arrangement_2                Arrangement_2;

// Merge mergable edges
void simplify(Arrangement_2& arr)
{
  Arrangement_2::Vertex_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() != 2) continue;
    Arrangement_2::Halfedge_around_vertex_circulator eit =
      vit->incident_halfedges();
    const Arrangement_2::Geometry_traits_2* traits = arr.geometry_traits();
    if (traits->are_mergeable_2_object()(eit->curve(), eit->next()->curve())) {
      Arrangement_2::Geometry_traits_2::X_monotone_curve_2 c;
      traits->merge_2_object()(eit->curve(), eit->next()->curve(), c);
      arr.merge_edge(eit, eit->next(), c);
    }
  }
}

/*! Check if two polygons with holes are the same. */
bool are_equal(Polygon_set_2& ps1, const Polygon_set_2& ps2)
{
  ps1.symmetric_difference(ps2);
  return (ps1.is_empty());
}

/*! The main program. */
int main(int argc, char* argv[])
{
  // Read the input file. Because of the structure of the *.cmd file
  // (which is concatenated to the command line) we need to get all the
  // inputs in one command line. This is the reason we read triplets/pairs of
  // arguments. Each triplet/double is one input for the program.
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << ". The input are triplets of:"
	      << " <polygon#1> <polygon#2> [decomposition flags]"
	      << std::endl;
    return (1);
  }

  int i = 1;
  while (i < argc) {
    // Read the polygons from the input files.
    Polygon_2 pgn1, pgn2;

    if (! read_polygon(argv[i], pgn1)) {
      std::cerr << "Failed to read: <" << argv[i] << ">." << std::endl;
      return (1);
    }

    if (! read_polygon(argv[i+1], pgn2)) {
      std::cerr << "Failed to read: <" << argv[i+1] << ">." << std::endl;
      return (1);
    }

    std::cout << "Testing " << argv[i] << " and " << argv[i+1] << std::endl;

    // Read the decomposition flags.
    bool use_ssab = true;
    bool use_opt = true;
    bool use_hm = true;
    bool use_greene = true;
    bool use_vertical = true;

    if (((i+2) < argc) && (argv[i+2][0] == '-')) {
      use_ssab = (std::strchr(argv[i+2], 's') != NULL);
      use_opt = (std::strchr(argv[i+2], 'o') != NULL);
      use_hm = (std::strchr(argv[i+2], 'h') != NULL);
      use_greene = (std::strchr(argv[i+2], 'g') != NULL);
      use_vertical = (std::strchr(argv[i+2], 'v') != NULL);
    }

    // Compute the Minkowski sum using the convolution method.
    std::cout << "Using the convolution method ... " << std::flush;
    Polygon_with_holes_2 sum_conv = minkowski_sum_2(pgn1, pgn2);
    std::cout << "simplifying ... " << std::flush;
    Polygon_set_2 ps_conv;
    ps_conv.insert(sum_conv);
    Arrangement_2& arr = ps_conv.arrangement();
    simplify(arr);
    std::cout << "Done." << std::endl;

    // Define auxiliary polygon-decomposition objects.
    CGAL::Small_side_angle_bisector_decomposition_2<Kernel>  ssab_decomp;
    CGAL::Optimal_convex_decomposition_2<Kernel>             opt_decomp;
    CGAL::Hertel_Mehlhorn_convex_decomposition_2<Kernel>     hm_approx_decomp;
    CGAL::Greene_convex_decomposition_2<Kernel>              greene_decomp;
    CGAL::Polygon_vertical_decomposition_2<Kernel>           vertical_decomp;

    if (use_ssab) {
      std::cout << "Using the small-side angle-bisector decomposition ... "
                << std::flush;
      Polygon_with_holes_2 sum = minkowski_sum_2(pgn1, pgn2, ssab_decomp);
      std::cout << "simplifying ... " << std::flush;
      Polygon_set_2 ps_decomp;
      ps_decomp.insert(sum);
      Arrangement_2& arr = ps_decomp.arrangement();
      simplify(arr);
      if (are_equal(ps_decomp, ps_conv)) {
        std::cout << "OK." << std::endl;
      }
      else {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_opt) {
      std::cout << "Using the optimal convex decomposition ... " << std::flush;
      Polygon_with_holes_2 sum = minkowski_sum_2(pgn1, pgn2, opt_decomp);
      std::cout << "simplifying ... " << std::flush;
      Polygon_set_2 ps_decomp;
      ps_decomp.insert(sum);
      Arrangement_2& arr = ps_decomp.arrangement();
      simplify(arr);
      if (are_equal(ps_decomp, ps_conv)) {
        std::cout << "OK." << std::endl;
      }
      else {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_hm) {
      std::cout << "Using the Hertel--Mehlhorn decomposition ... "
                << std::flush;
      Polygon_with_holes_2 sum = minkowski_sum_2(pgn1, pgn2, hm_approx_decomp);
      std::cout << "simplifying ... " << std::flush;
      Polygon_set_2 ps_decomp;
      ps_decomp.insert(sum);
      Arrangement_2& arr = ps_decomp.arrangement();
      simplify(arr);
      if (are_equal(ps_decomp, ps_conv)) {
        std::cout << "OK." << std::endl;
      }
      else {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_greene) {
      std::cout << "Using the Greene decomposition ... " << std::flush;
      Polygon_with_holes_2 sum = minkowski_sum_2(pgn1, pgn2, greene_decomp);
      std::cout << "simplifying ... " << std::flush;
      Polygon_set_2 ps_decomp;
      ps_decomp.insert(sum);
      Arrangement_2& arr = ps_decomp.arrangement();
      simplify(arr);
       if (are_equal(ps_decomp, ps_conv)) {
        std::cout << "OK." << std::endl;
      }
      else {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (use_vertical) {
      std::cout << "Using the vertical decomposition ... " << std::flush;
      Polygon_with_holes_2 sum = minkowski_sum_2(pgn1, pgn2, vertical_decomp);
      std::cout << "simplifying ... " << std::flush;
      Polygon_set_2 ps_decomp;
      ps_decomp.insert(sum);
      Arrangement_2& arr = ps_decomp.arrangement();
      simplify(arr);
      if (are_equal(ps_decomp, ps_conv)) {
        std::cout << "OK." << std::endl;
      }
      else {
        std::cout << "ERROR (different result)." << std::endl;
        return 1;
      }
    }

    if (((i + 2) < argc) && (argv[i+2][0] == '-')) i += 3;
    else i += 2;
  }

  return (0);
}
