// Checks that Kd_tree::build<Parallel_tag>() builds the same tree as
// Kd_tree::build<Sequential_tag>(), over the seven splitters, three traits and
// two point distributions, at several thread counts. Trees are compared by item
// count, node count, depth and a hash of the tree structure and of the points of
// every leaf; the order of the points inside a leaf may differ between the two
// builds.
//
// Optional argument: a point count, which additionally prints build times.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/use.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#  include <tbb/global_control.h>
#endif

// ---------------------------------------------------------------- signature

struct Signature {
  std::size_t num_items;
  std::size_t num_nodes;
  int depth;
  int dim;
  std::uint64_t structure;

  bool operator==(const Signature& o) const
  {
    return num_items == o.num_items && num_nodes == o.num_nodes
      && depth == o.depth && dim == o.dim && structure == o.structure;
  }
};

std::ostream& operator<<(std::ostream& s, const Signature& t)
{
  s << "items=" << t.num_items << " nodes=" << t.num_nodes
    << " depth=" << t.depth << " dim=" << t.dim
    << " structure=" << t.structure;
  return s;
}

inline void mix(std::uint64_t& h, std::uint64_t v)
{
  for (int i = 0; i != 8; ++i) {
    h ^= (v >> (8 * i)) & 0xffu;
    h *= 1099511628211ull;
  }
}

template <class FT>
std::uint64_t bits(const FT& x)
{
  const double v = CGAL::to_double(x);
  std::uint64_t b = 0;
  std::memcpy(&b, &v, sizeof(double));
  return b;
}

// Pre-order hash of the tree: the leaf/internal flag, size and points of every
// leaf, and the cutting dimension and value of every internal node. Node and
// item counts alone do not distinguish two trees of the same size. The hash of
// the points of a leaf does not depend on their order.
template <class Tree>
void hash_subtree(const Tree& tree, typename Tree::Node_const_handle n, std::uint64_t& h)
{
  if (n->is_leaf()) {
    typename Tree::Leaf_node_const_handle leaf
      = static_cast<typename Tree::Leaf_node_const_handle>(n);
    mix(h, 1);
    mix(h, leaf->size());
    typename Tree::Traits::Construct_cartesian_const_iterator_d construct_it
      = tree.traits().construct_cartesian_const_iterator_d_object();
    std::uint64_t points = 0;
    for (typename Tree::iterator it = leaf->begin(); it != leaf->end(); ++it) {
      std::uint64_t p = 14695981039346656037ull;
      typename Tree::Traits::Cartesian_const_iterator_d c = construct_it(*it);
      for (int i = 0; i != tree.dim(); ++i, ++c)
        mix(p, bits(*c));
      points += p;
    }
    mix(h, points);
  }
  else {
    typename Tree::Internal_node_const_handle node
      = static_cast<typename Tree::Internal_node_const_handle>(n);
    mix(h, 0);
    mix(h, static_cast<std::uint64_t>(node->cutting_dimension()));
    mix(h, bits(node->cutting_value()));
    hash_subtree(tree, node->lower(), h);
    hash_subtree(tree, node->upper(), h);
  }
}

// ------------------------------------------------------------------- driver

enum Mode { SEQUENTIAL, PARALLEL };

template <class Traits, class Splitter, class PointRange>
Signature build_and_sign(const PointRange& points, Mode mode,
                         double* seconds = nullptr)
{
  typedef CGAL::Kd_tree<Traits, Splitter> Tree;

  Tree tree(points.begin(), points.end(), Splitter(10));

  CGAL::Real_timer timer;
  timer.start();
  if (mode == PARALLEL)
    tree.template build<CGAL::Parallel_tag>();
  else
    tree.build();
  timer.stop();
  if (seconds != nullptr)
    *seconds = timer.time();

  const Tree& ctree = tree;
  Signature sig;
  sig.num_items = ctree.root()->num_items();
  sig.num_nodes = ctree.root()->num_nodes();
  sig.depth = ctree.root()->depth();
  sig.dim = ctree.dim();
  sig.structure = 14695981039346656037ull;
  hash_subtree(ctree, ctree.root(), sig.structure);
  return sig;
}

int failures = 0;

template <class Traits, class Splitter, class PointRange>
void check_one(const std::string& traits_name,
               const std::string& splitter_name,
               const std::string& input_name,
               const PointRange& points,
               const std::vector<int>& thread_counts)
{
  const Signature expected
    = build_and_sign<Traits, Splitter>(points, SEQUENTIAL);

  for (std::size_t i = 0; i != thread_counts.size(); ++i) {
    const int nb = thread_counts[i];
#ifdef CGAL_LINKED_WITH_TBB
    tbb::global_control control(tbb::global_control::max_allowed_parallelism,
                                static_cast<std::size_t>(nb));
#endif

    const Signature got = build_and_sign<Traits, Splitter>(points, PARALLEL);

    if (!(got == expected)) {
      ++failures;
      std::cerr << "MISMATCH  " << traits_name << " / " << splitter_name
                << " / " << input_name << " / " << nb << " threads\n"
                << "  sequential: " << expected << '\n'
                << "  parallel  : " << got << std::endl;
    }
  }
}

// Fair and Sliding_fair divide by a bounding box side length, which is zero on
// heavily duplicated points. Only a plain double yields inf there; Quotient and
// the exact kernels report a division by zero, so those two splitters are run on
// the uniform input only.
enum Splitter_set { ALL_SPLITTERS, WITHOUT_ASPECT_RATIO_SPLITTERS };

template <class Traits, class PointRange>
void check_splitters(const std::string& traits_name,
                     const std::string& input_name,
                     const PointRange& points,
                     const std::vector<int>& thread_counts,
                     Splitter_set which = ALL_SPLITTERS)
{
  check_one<Traits, CGAL::Median_of_max_spread<Traits> >
    (traits_name, "Median_of_max_spread", input_name, points, thread_counts);
  check_one<Traits, CGAL::Sliding_midpoint<Traits> >
    (traits_name, "Sliding_midpoint", input_name, points, thread_counts);
  check_one<Traits, CGAL::Median_of_rectangle<Traits> >
    (traits_name, "Median_of_rectangle", input_name, points, thread_counts);
  check_one<Traits, CGAL::Midpoint_of_max_spread<Traits> >
    (traits_name, "Midpoint_of_max_spread", input_name, points, thread_counts);
  check_one<Traits, CGAL::Midpoint_of_rectangle<Traits> >
    (traits_name, "Midpoint_of_rectangle", input_name, points, thread_counts);

  if (which == ALL_SPLITTERS) {
    check_one<Traits, CGAL::Fair<Traits> >
      (traits_name, "Fair", input_name, points, thread_counts);
    check_one<Traits, CGAL::Sliding_fair<Traits> >
      (traits_name, "Sliding_fair", input_name, points, thread_counts);
  }
}

// ------------------------------------------------------------------- inputs

// Snapping to a coarse grid makes many points coincide, so that subtrees reach
// zero tight spread and become leaves. A grid cell has to hold more than
// bucket_size points, and a few cells more than the size from which build()
// forks, or no such subtree is ever formed.
inline double snap(double v, double step)
{
  return step * std::floor(v / step);
}

template <class K>
std::vector<typename K::Point_3>
points_3(std::size_t n, bool duplicate_heavy, CGAL::Random& rnd, double step)
{
  typedef typename K::Point_3 Point_3;
  std::vector<Point_3> points;
  points.reserve(n);

  CGAL::Random_points_in_cube_3<CGAL::Simple_cartesian<double>::Point_3>
    gen(1.0, rnd);

  for (std::size_t i = 0; i != n; ++i, ++gen) {
    double x = (*gen).x(), y = (*gen).y(), z = (*gen).z();
    if (duplicate_heavy) {
      x = snap(x, step); y = snap(y, step); z = snap(z, step);
    }
    points.push_back(Point_3(x, y, z));
  }
  return points;
}

template <class K>
std::vector<typename K::Point_d>
points_d(std::size_t n, int dim, bool duplicate_heavy, CGAL::Random& rnd,
         double step)
{
  typedef typename K::Point_d Point_d;
  std::vector<Point_d> points;
  points.reserve(n);

  std::vector<double> coord(static_cast<std::size_t>(dim));
  for (std::size_t i = 0; i != n; ++i) {
    for (int j = 0; j != dim; ++j) {
      const double v = rnd.get_double(-1.0, 1.0);
      coord[static_cast<std::size_t>(j)] = duplicate_heavy ? snap(v, step) : v;
    }
    points.push_back(Point_d(dim, coord.begin(), coord.end()));
  }
  return points;
}

// ------------------------------------------------------------------- timing

template <class Traits, class Splitter, class PointRange>
void report_times(const PointRange& points,
                  const std::vector<int>& thread_counts)
{
  // One discarded build first, so that the first timed run does not pay for
  // faulting in the point vector.
  build_and_sign<Traits, Splitter>(points, SEQUENTIAL);

  double seq = 0.;
  build_and_sign<Traits, Splitter>(points, SEQUENTIAL, &seq);
  std::cout << "  sequential           \t" << seq << " s" << std::endl;

  for (std::size_t i = 0; i != thread_counts.size(); ++i) {
#ifdef CGAL_LINKED_WITH_TBB
    tbb::global_control control(tbb::global_control::max_allowed_parallelism,
                                static_cast<std::size_t>(thread_counts[i]));
#endif
    double par = 0.;
    build_and_sign<Traits, Splitter>(points, PARALLEL, &par);
    std::cout << "  parallel, " << thread_counts[i] << " threads\t" << par
              << " s\tspeedup " << (par > 0. ? seq / par : 0.) << std::endl;
  }
}

// --------------------------------------------------------------------- main

int main(int argc, char** argv)
{
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_USE(argc);
  CGAL_USE(argv);
  std::cout << "NOTICE: TBB is not linked, Parallel_tag is unavailable."
            << std::endl;
  return 0;
#else

  typedef CGAL::Simple_cartesian<double> Sc;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
  typedef CGAL::Homogeneous_d<double> Kd;

  typedef CGAL::Search_traits_3<Sc> Traits_sc_3;
  typedef CGAL::Search_traits_3<Epeck> Traits_epeck_3;
  typedef Kd Traits_dynamic;  // a d-dimensional kernel is itself a SearchTraits

  std::vector<int> thread_counts;
  thread_counts.push_back(1);
  thread_counts.push_back(4);

  const std::size_t n = 2000;
  const int dim = 5;
  CGAL::Random rnd(42);

  {
    const std::vector<Sc::Point_3> uni = points_3<Sc>(n, false, rnd, 1.0);
    const std::vector<Sc::Point_3> dup = points_3<Sc>(n, true, rnd, 1.0);
    check_splitters<Traits_sc_3>("Search_traits_3<Simple_cartesian>",
                                 "uniform", uni, thread_counts);
    check_splitters<Traits_sc_3>("Search_traits_3<Simple_cartesian>",
                                 "duplicate-heavy", dup, thread_counts,
                                 WITHOUT_ASPECT_RATIO_SPLITTERS);
  }
  {
    const std::vector<Epeck::Point_3> uni = points_3<Epeck>(n, false, rnd, 1.0);
    const std::vector<Epeck::Point_3> dup = points_3<Epeck>(n, true, rnd, 1.0);
    check_splitters<Traits_epeck_3>("Search_traits_3<Epeck>",
                                    "uniform", uni, thread_counts);
    check_splitters<Traits_epeck_3>("Search_traits_3<Epeck>",
                                    "duplicate-heavy", dup, thread_counts,
                                    WITHOUT_ASPECT_RATIO_SPLITTERS);
  }
  {
    const std::vector<Kd::Point_d> uni = points_d<Kd>(n, dim, false, rnd, 1.0);
    const std::vector<Kd::Point_d> dup = points_d<Kd>(n, dim, true, rnd, 1.0);
    check_splitters<Traits_dynamic>("Homogeneous_d", "uniform",
                                    uni, thread_counts);
    check_splitters<Traits_dynamic>("Homogeneous_d", "duplicate-heavy",
                                    dup, thread_counts,
                                    WITHOUT_ASPECT_RATIO_SPLITTERS);
  }

  // build() computes the bounding boxes, the partition and the median of a
  // node of more than 1000 points in parallel. These cases reach all three,
  // on one splitter of each kind.
  {
    const std::vector<int> many(1, 32);
    const std::vector<Sc::Point_3> uni = points_3<Sc>(300000, false, rnd, 0.5);
    const std::vector<Sc::Point_3> dup = points_3<Sc>(300000, true, rnd, 0.5);
    const std::vector<Sc::Point_3> dup_big = points_3<Sc>(600000, true, rnd, 0.5);
    check_one<Traits_sc_3, CGAL::Sliding_midpoint<Traits_sc_3> >
      ("Search_traits_3<Simple_cartesian>", "Sliding_midpoint", "uniform, 300000",
       uni, many);
    check_one<Traits_sc_3, CGAL::Sliding_midpoint<Traits_sc_3> >
      ("Search_traits_3<Simple_cartesian>", "Sliding_midpoint", "duplicate-heavy, 300000",
       dup, many);
    check_one<Traits_sc_3, CGAL::Median_of_rectangle<Traits_sc_3> >
      ("Search_traits_3<Simple_cartesian>", "Median_of_rectangle", "duplicate-heavy, 600000",
       dup_big, many);
    const std::vector<Epeck::Point_3> exact = points_3<Epeck>(70000, true, rnd, 0.5);
    check_one<Traits_epeck_3, CGAL::Median_of_rectangle<Traits_epeck_3> >
      ("Search_traits_3<Epeck>", "Median_of_rectangle", "duplicate-heavy, 70000",
       exact, many);
    const std::vector<Kd::Point_d> dyn = points_d<Kd>(70000, dim, true, rnd, 1.0);
    check_one<Traits_dynamic, CGAL::Median_of_rectangle<Traits_dynamic> >
      ("Homogeneous_d", "Median_of_rectangle", "duplicate-heavy, 70000",
       dyn, many);
  }

  if (failures != 0) {
    std::cerr << failures << " mismatch(es)." << std::endl;
    return 1;
  }
  std::cout << "Parallel and sequential builds agree." << std::endl;

  if (argc > 1) {
    const std::size_t big = static_cast<std::size_t>(std::atol(argv[1]));
    std::cout << "\nBuild times, Search_traits_3<Simple_cartesian>, "
              << "Sliding_midpoint, " << big << " uniform points:" << std::endl;
    const std::vector<Sc::Point_3> pts = points_3<Sc>(big, false, rnd, 1.0);
    report_times<Traits_sc_3, CGAL::Sliding_midpoint<Traits_sc_3> >
      (pts, thread_counts);
  }

  return 0;
#endif
}
