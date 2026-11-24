#define CGAL_NO_DEPRECATION_WARNINGS 1 // No_intersection_tag is tested
#define CGAL_NO_CDT_2_WARNING 1        // Exact_intersections_tag is tested with Constrained_triangulation_2

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Polygon_2.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;

int nb_failures = 0;

void failure(const std::string& what)
{
  std::cerr << "FAILURE: " << what << std::endl;
}

template <class K>
struct Intersection_tag_test
{
  using Point = typename K::Point_2;
  using Polyline = std::vector<Point>;     // closed if its extremities are equal
  using Polylines = std::vector<Polyline>; // inserted in every order
  using Edge = std::pair<Point, Point>;

  struct Configuration
  {
    std::string name;
    std::vector<Polylines> polylines; // inserted one after the other
    std::vector<Edge> edges;          // the expected constrained edges
  };

  static Polyline pt(double x, double y) { return {Point(x, y)}; }
  static Polyline seg(double x0, double y0, double x1, double y1) { return {Point(x0, y0), Point(x1, y1)}; }
  static Polyline pl(std::initializer_list<double> xy)
  {
    Polyline p;
    for(auto it = xy.begin(); it != xy.end(); it += 2)
      p.emplace_back(*it, *(it + 1));
    return p;
  }
  static Edge ce(const Point& p, const Point& q) { return (q < p) ? Edge(q, p) : Edge(p, q); }
  static Edge ce(double x0, double y0, double x1, double y1) { return ce(Point(x0, y0), Point(x1, y1)); }

  // ------------------------------------------------------------------------------------------
  // Configurations
  // ------------------------------------------------------------------------------------------

  // (a), (b), (c): constraints that do not intersect, or only at a common extremity
  const std::vector<Configuration> no_intersection =
  {
    {"(a) disjoint", {{seg(0,0, 4,0), seg(0,2, 4,3)}}, {ce(0,0, 4,0), ce(0,2, 4,3)}},
    {"(a) disjoint, collinear", {{seg(0,0, 1,0), seg(2,0, 4,0)}}, {ce(0,0, 1,0), ce(2,0, 4,0)}},
    {"(a) isolated vertex", {{seg(0,0, 4,0), pt(2,1)}}, {ce(0,0, 4,0)}},

    {"(b) common extremity", {{seg(0,0, 2,2), seg(2,2, 4,0)}}, {ce(0,0, 2,2), ce(2,2, 4,0)}},
    {"(b) common extremity, collinear", {{seg(0,0, 2,0), seg(2,0, 4,0)}}, {ce(0,0, 2,0), ce(2,0, 4,0)}},
    {"(b) three constraints with a common extremity", {{seg(2,2, 0,0), seg(2,2, 4,0), seg(2,2, 2,4)}},
      {ce(0,0, 2,2), ce(2,2, 4,0), ce(2,2, 2,4)}},
    {"(b) vertex at an extremity", {{seg(0,0, 4,0), pt(0,0)}}, {ce(0,0, 4,0)}},
    {"(b) zero-length constraint", {{seg(1,1, 1,1)}}, {}},
    {"(b) polyline with repeated points", {{pl({0,0, 0,0, 2,0, 2,0, 2,2})}}, {ce(0,0, 2,0), ce(2,0, 2,2)}},
    {"(b) nested closed polylines", {{pl({0,0, 6,0, 6,6, 0,6, 0,0}), pl({2,2, 4,2, 3,4, 2,2})}},
      {ce(0,0, 6,0), ce(6,0, 6,6), ce(6,6, 0,6), ce(0,6, 0,0), ce(2,2, 4,2), ce(4,2, 3,4), ce(3,4, 2,2)}},

    {"(c) closed polyline pinch", {{pl({2,2, 4,4, 4,0, 2,2, 0,0, 0,4, 2,2})}},
      {ce(2,2, 4,4), ce(4,4, 4,0), ce(4,0, 2,2), ce(2,2, 0,0), ce(0,0, 0,4), ce(0,4, 2,2)}},
    {"(c) open polyline pinch", {{pl({0,0, 2,2, 4,0, 4,4, 2,2, 0,4})}},
      {ce(0,0, 2,2), ce(2,2, 4,0), ce(4,0, 4,4), ce(4,4, 2,2), ce(2,2, 0,4)}},
    {"(c) closed polylines sharing a vertex", {{pl({0,0, 2,2, 0,4, 0,0}), pl({2,2, 4,0, 4,4, 2,2})}},
      {ce(0,0, 2,2), ce(2,2, 0,4), ce(0,4, 0,0), ce(2,2, 4,0), ce(4,0, 4,4), ce(4,4, 2,2)}},
  };

  // (d), (e), (f), (g): intersections at existing vertices (T-junctions, vertices in constraints,
  // overlapping and identical constraints); constraints are split at these vertices
  const std::vector<Configuration> intersection_at_vertices =
  {
    {"(d) T-junction", {{seg(0,0, 4,0), seg(2,0, 2,2)}}, {ce(0,0, 2,0), ce(2,0, 4,0), ce(2,0, 2,2)}},
    {"(d) double T-junction", {{seg(0,0, 0,4), seg(4,0, 4,4), seg(0,2, 4,2)}},
      {ce(0,0, 0,2), ce(0,2, 0,4), ce(4,0, 4,2), ce(4,2, 4,4), ce(0,2, 4,2)}},
    {"(d) constraint through a common extremity", {{seg(0,0, 2,0), seg(2,0, 4,0), seg(2,-2, 2,2)}},
      {ce(0,0, 2,0), ce(2,0, 4,0), ce(2,-2, 2,0), ce(2,0, 2,2)}},
    {"(d) constraint through a polyline vertex", {{pl({0,0, 2,2, 4,0}), seg(2,4, 2,0)}},
      {ce(0,0, 2,2), ce(2,2, 4,0), ce(2,4, 2,2), ce(2,2, 2,0)}},
    {"(d) polyline touching itself", {{pl({0,0, 4,0, 4,2, 2,0})}},
      {ce(0,0, 2,0), ce(2,0, 4,0), ce(4,0, 4,2), ce(4,2, 2,0)}},
    {"(d) hole touching the outer boundary", {{pl({0,0, 6,0, 6,6, 0,6, 0,0}), pl({3,0, 4,2, 2,2, 3,0})}},
      {ce(0,0, 3,0), ce(3,0, 6,0), ce(6,0, 6,6), ce(6,6, 0,6), ce(0,6, 0,0), ce(3,0, 4,2), ce(4,2, 2,2), ce(2,2, 3,0)}},

    {"(e) vertex in a constraint", {{seg(0,0, 4,0), pt(2,0)}}, {ce(0,0, 2,0), ce(2,0, 4,0)}},
    {"(e) vertices in a diagonal constraint", {{seg(0,0,4,4), pt(1,1), pt(3,3)}},
      {ce(0,0,1,1), ce(1,1,3,3), ce(3,3,4,4)}},
    {"(e) zero-length constraint in a constraint", {{seg(0,0, 4,0), seg(2,0, 2,0)}},
      {ce(0,0, 2,0), ce(2,0, 4,0)}},

    {"(f) collinear overlap", {{seg(0,0, 3,0), seg(1,0, 4,0)}}, {ce(0,0, 1,0), ce(1,0, 3,0), ce(3,0, 4,0)}},
    {"(f) containment", {{seg(0,0, 4,4), seg(1,1, 3,3)}}, {ce(0,0, 1,1), ce(1,1, 3,3), ce(3,3, 4,4)}},
    {"(f) overlap with a common extremity", {{seg(0,0, 4,0), seg(0,0, 2,0)}}, {ce(0,0, 2,0), ce(2,0, 4,0)}},
    {"(f) polyline folding back", {{pl({0,0, 4,0, 2,0})}}, {ce(0,0, 2,0), ce(2,0, 4,0)}},
    {"(f) closed polylines sharing part of an edge", {{pl({0,0, 2,0, 2,2, 0,2, 0,0}), pl({2,1, 4,1, 4,3, 2,3, 2,1})}},
     {ce(0,0, 2,0), ce(2,0, 2,1), ce(2,1, 2,2), ce(2,2, 0,2), ce(0,2, 0,0), ce(2,1, 4,1), ce(4,1, 4,3), ce(4,3, 2,3), ce(2,3, 2,2)}},

    {"(g) opposite constraints", {{seg(0,0, 4,0), seg(4,0, 0,0)}}, {ce(0,0, 4,0)}},
    {"(g) identical constraint inserted later", {{seg(0,0, 4,0)}, {seg(0,0, 4,0)}}, {ce(0,0, 4,0)}},
    {"(g) polyline going back and forth", {{pl({0,0, 4,0, 0,0})}}, {ce(0,0, 4,0)}},
    {"(g) closed polylines sharing an edge", {{pl({0,0, 2,0, 2,2, 0,2, 0,0}), pl({2,0, 4,0, 4,2, 2,2, 2,0})}},
     {ce(0,0, 2,0), ce(2,0, 2,2), ce(2,2, 0,2), ce(0,2, 0,0), ce(2,0, 4,0), ce(4,0, 4,2), ce(4,2, 2,2)}},
  };

  // (h): crossings; the expected edges are those of the tags that construct the crossing points
  const std::vector<Configuration> crossing =
  {
    {"(h) crossing", {{seg(0,0, 4,4), seg(0,4, 4,0)}},
      {ce(0,0, 2,2), ce(2,2, 4,4), ce(0,4, 2,2), ce(2,2, 4,0)}},
    {"(h) crossing at a non-lattice point", {{seg(0,0, 3,1), seg(0,1, 3,0)}},
      {ce(0,0, 1.5,0.5), ce(1.5,0.5, 3,1), ce(0,1, 1.5,0.5), ce(1.5,0.5, 3,0)}},
    {"(h) three constraints crossing at a point", {{seg(0,0, 4,4), seg(0,4, 4,0), seg(2,0, 2,4)}},
      {ce(0,0, 2,2), ce(2,2, 4,4), ce(0,4, 2,2), ce(2,2, 4,0), ce(2,0, 2,2), ce(2,2, 2,4)}},
    {"(h) self-crossing closed polyline", {{pl({0,0, 4,4, 4,0, 0,4, 0,0})}},
      {ce(0,0, 2,2), ce(2,2, 4,4), ce(4,4, 4,0), ce(4,0, 2,2), ce(2,2, 0,4), ce(0,4, 0,0)}},
  };

  // (i): crossings at a point that is a vertex when the second constraint is inserted
  const std::vector<Edge> x = {ce(0,0, 2,2), ce(2,2, 4,4), ce(0,4, 2,2), ce(2,2, 4,0)};
  const std::vector<Edge> x_and_stem = {ce(0,0, 2,2), ce(2,2, 4,4), ce(0,4, 2,2), ce(2,2, 4,0), ce(2,2, 2,5)};
  const std::vector<Configuration> crossing_at_vertex =
  {
    {"(i) crossing at a vertex inserted before", {{pt(2,2)}, {seg(0,0, 4,4), seg(0,4, 4,0)}}, x},
    {"(i) crossing at a vertex inserted in between", {{seg(0,0,4,4)}, {pt(2,2)}, {seg(0,4,4,0)}}, x},
    {"(i) crossing at the extremity of a constraint", {{seg(2,2, 2,5)}, {seg(0,0, 4,4), seg(0,4, 4,0)}}, x_and_stem},
    {"(i) crossing at a T-junction", {{seg(0,0, 4,4), seg(2,2, 2,5)}, {seg(0,4, 4,0)}}, x_and_stem},
    {"(i) polyline crossing at one of its vertices", {{pl({0,0, 4,4, 4,0, 0,4, 0,2, 2,2})}},
     {ce(0,0, 2,2), ce(2,2, 4,4), ce(4,4, 4,0), ce(4,0, 2,2), ce(2,2, 0,4), ce(0,4, 0,2), ce(0,2, 2,2)}},
  };

  // (i): crossings at a point that becomes a vertex after the second constraint is inserted
  const std::vector<Configuration> crossing_at_later_vertex =
  {
    {"(i) crossing at a vertex inserted after", {{seg(0,0, 4,4), seg(0,4, 4,0)}, {pt(2,2)}}, x},
    {"(i) crossing at the extremity of a constraint inserted after", {{seg(0,0, 4,4), seg(0,4, 4,0)}, {seg(2,2, 2,5)}}, x_and_stem},
  };

  // ------------------------------------------------------------------------------------------
  // Checks
  // ------------------------------------------------------------------------------------------

  template <class Tr>
  struct is_ctp : std::false_type {};

  template <class Tr>
  struct is_ctp<CGAL::Constrained_triangulation_plus_2<Tr>> : std::true_type {};

  // next order of the items of the polylines
  static bool next_order(std::vector<Polylines>& polylines)
  {
    for(Polylines& s : polylines)
      if(std::next_permutation(s.begin(), s.end()))
        return true;
    return false;
  }

  // Inserts the polylines of `c` with `insert(t, step)`, in every order of each step, with or without an
  // additional vertex (configurations of dimension 1 are thus also tested in dimension 2).
  // Checks that an exception is thrown if `exception` is true, and the result otherwise.
  template <class Tr, class Insert>
  void check(const Configuration& c,
             bool exception,
             const std::string& function,
             Insert insert)
  {
    for(bool additional_vertex : {false, true})
    {
      std::vector<Polylines> polylines = c.polylines;
      for(Polylines& s : polylines)
        std::sort(s.begin(), s.end());

      do
      {
        const std::string what = c.name + " with " + function;
        Tr t;
        if(additional_vertex)
          t.insert(Point(13, 37));

        bool thrown = false;
        try
        {
          for(const Polylines& s : polylines)
            insert(t, s);
        }
        catch(const typename Tr::Intersection_of_constraints_exception&)
        {
          thrown = true;
        }

        if(thrown != exception)
        {
          failure((exception ? "no exception: " : "unexpected exception: ") + what);
          continue;
        }

        if(thrown)
          continue;

        std::set<Point> vertices, expected_vertices;
        std::set<Edge> edges, expected_edges(c.edges.begin(), c.edges.end());
        for(auto v : t.finite_vertex_handles())
          vertices.insert(v->point());

        if(additional_vertex)
          expected_vertices.insert(Point(13, 37));

        for(const Polylines& s : polylines)
          for(const Polyline& p : s)
            expected_vertices.insert(p.begin(), p.end());

        for(const Edge& ed : c.edges)
          expected_vertices.insert({ed.first, ed.second});

        for(auto ed : t.constrained_edges())
          edges.insert(ce(ed.first->vertex(t.cw(ed.second))->point(),
                          ed.first->vertex(t.ccw(ed.second))->point()));

        if(!t.is_valid())
          failure("invalid triangulation: " + what);

        if(vertices != expected_vertices)
          failure("wrong vertices: " + what);

        if(edges != expected_edges)
          failure("wrong constrained edges: " + what);

        if constexpr(is_ctp<Tr>::value)
        {
          std::set<Edge> subconstraints;
          for(const auto& sc : t.subconstraints())
            subconstraints.insert(ce(sc.first->point(), sc.second->point()));

          if(subconstraints != edges)
            failure("subconstraints are not the constrained edges: " + what);
        }
      } while(next_order(polylines));
    }
  }

  // ------------------------------------------------------------------------------------------
  // Specific tests for each constrained class, to use specific insertion functions

  template <class CT>
  void test_constrained_triangulation_2(const Configuration& c, bool exception)
  {
    check<CT>(c, exception, "insert()",
              [](CT& t, const Polylines& steps)
              {
                for(const Polyline& p : steps)
                {
                  if(p.size() == 1)
                    t.insert(p[0]);
                  else if(p.size() == 2)
                    t.insert_constraint(p[0], p[1]);
                  else
                    t.insert_constraint(p.begin(), p.end());
                }
              });

    check<CT>(c, exception, "push_back()",
              [](CT& t, const Polylines& steps)
              {
                for(const Polyline& p : steps)
                {
                  if(p.size() == 1)
                    t.push_back(p[0]);
                  else if(p.size() == 2)
                    t.push_back(std::make_pair(p[0], p[1]));
                  else if(p.front() == p.back())
                    t.insert_constraint(p.begin(), p.end() - 1, true); // closed polyline
                  else
                    t.insert_constraint(p.begin(), p.end(), false);
                }
              });

    check<CT>(c, exception, "vertex handles",
              [](CT& t, const Polylines& steps)
              {
                for(const Polyline& p : steps)
                {
                  if(p.size() == 1)
                  {
                    t.insert_constraint(p[0], p[0]); // a zero-length constraint is a point
                    continue;
                  }
                  std::vector<typename CT::Vertex_handle> v;
                  for(const Point& q : p)
                    v.push_back(t.insert(q));
                  for(std::size_t i=1; i<v.size(); ++i)
                    t.insert_constraint(v[i - 1], v[i]);
                }
              });
  }

  template <class CDT>
  void test_constrained_Delaunay_triangulation_2(const Configuration& c, bool exception)
  {
    test_constrained_triangulation_2<CDT>(c, exception);

    // additional stuff
    check<CDT>(c, exception, "insert_constraints()",
               [](CDT& t, const Polylines& steps)
               {
                 std::vector<Point> points;
                 std::vector<std::pair<Point, Point>> segments;
                 for(const Polyline& p : steps)
                 {
                   if(p.size() == 1)
                     points.push_back(p[0]);
                   for(std::size_t i=1; i<p.size(); ++i)
                     segments.emplace_back(p[i - 1], p[i]);
                 }
                 t.insert(points.begin(), points.end());
                 t.insert_constraints(segments.begin(), segments.end());
               });

    check<CDT>(c, exception, "insert_constraints() with indices",
               [](CDT& t, const Polylines& steps)
               {
                 std::vector<Point> points;
                 std::vector<std::pair<std::size_t, std::size_t>> indices;
                 for(const Polyline& p : steps)
                 {
                   for(std::size_t i=0; i<p.size(); ++i)
                   {
                     points.push_back(p[i]);
                     if(i > 0)
                       indices.emplace_back(points.size() - 2, points.size() - 1);
                   }
                 }
                 t.insert_constraints(points.begin(), points.end(), indices.begin(), indices.end());
               });
  }

  template <class CTP>
  void test_constrained_triangulation_plus_2(const Configuration& c, bool exception)
  {
    test_constrained_Delaunay_triangulation_2<CTP>(c, exception);

    // additional stuff
    check<CTP>(c, exception, "insert_constraint(range)",
               [](CTP& t, const Polylines& steps)
               {
                 for(const Polyline& p : steps)
                 {
                   if(p.size() == 1)
                     t.insert(p[0]);
                   else if(p.size() > 2 && p.front() == p.back())
                     t.insert_constraint(CGAL::Polygon_2<K>(p.begin(), p.end() - 1));
                   else
                     t.insert_constraint(p);
                 }
               });
  }

  template <class Itag>
  void check(const std::vector<Configuration>& configurations, bool exception)
  {
    using CT = CGAL::Constrained_triangulation_2<K, CGAL::Default, Itag>;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>;

    for(const Configuration& c : configurations)
    {
      test_constrained_triangulation_2<CT>(c, exception);
      test_constrained_Delaunay_triangulation_2<CDT>(c, exception);
      test_constrained_triangulation_plus_2<CGAL::Constrained_triangulation_plus_2<CT> >(c, exception);
      test_constrained_triangulation_plus_2<CGAL::Constrained_triangulation_plus_2<CDT> >(c, exception);
    }
  }

  template <class Itag> void accepted(const std::vector<Configuration>& c)
  {
    check<Itag>(c, false);
  }
  template <class Itag> void rejected(const std::vector<Configuration>& c)
  {
    check<Itag>(c, true);
  }

  // ------------------------------------------------------------------------------------------
  // One test per intersection tag
  // ------------------------------------------------------------------------------------------

  void test_No_constraint_intersection_tag()
  {
    using Itag = CGAL::No_constraint_intersection_tag;

    accepted<Itag>(no_intersection);
    rejected<Itag>(intersection_at_vertices);
    rejected<Itag>(crossing);
    rejected<Itag>(crossing_at_vertex);
    rejected<Itag>(crossing_at_later_vertex);

    // insert_unique_constraints() only removes the constraints duplicated in its range
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>;
    const std::vector<Edge> twice = {{Point(0,0), Point(4,0)}, {Point(4,0), Point(0,0)}};
    CDT t;
    try
    {
      t.insert_unique_constraints(twice.begin(), twice.end());
    }
    catch(const typename CDT::Intersection_of_constraints_exception&)
    {
      failure("insert_unique_constraints() did not remove the duplicates");
    }

    try
    {
      t.insert_unique_constraints(twice.begin(), twice.begin() + 1);
      failure("insert_unique_constraints(): missing exception for a constraint already in the triangulation");
    }
    catch(const typename CDT::Intersection_of_constraints_exception&)
    {
    }
  }

  void test_No_constraint_intersection_requiring_constructions_tag()
  {
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;

    accepted<Itag>(no_intersection);
    accepted<Itag>(intersection_at_vertices);
    rejected<Itag>(crossing);
    accepted<Itag>(crossing_at_vertex);
    rejected<Itag>(crossing_at_later_vertex);
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  void test_No_intersection_tag() // deprecated, equivalent to No_constraint_intersection_requiring_constructions_tag
  {
    using Itag = CGAL::No_intersection_tag;

    accepted<Itag>(no_intersection);
    accepted<Itag>(intersection_at_vertices);
    rejected<Itag>(crossing);
    accepted<Itag>(crossing_at_vertex);
    rejected<Itag>(crossing_at_later_vertex);
  }
#endif

  void test_Exact_predicates_tag()
  {
    using Itag = CGAL::Exact_predicates_tag;

    accepted<Itag>(no_intersection);
    accepted<Itag>(intersection_at_vertices);
    accepted<Itag>(crossing);
    accepted<Itag>(crossing_at_vertex);
    accepted<Itag>(crossing_at_later_vertex);
  }

  void test_Exact_intersections_tag()
  {
    using Itag = CGAL::Exact_intersections_tag;

    accepted<Itag>(no_intersection);
    accepted<Itag>(intersection_at_vertices);
    accepted<Itag>(crossing);
    accepted<Itag>(crossing_at_vertex);
    accepted<Itag>(crossing_at_later_vertex);
  }
};

void test_exception()
{
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<Epick, CGAL::Default, CGAL::No_constraint_intersection_tag>;
  using Point = Epick::Point_2;

  static_assert(std::is_base_of_v<std::exception, CDT::Intersection_of_constraints_exception>);
  static_assert(std::is_same_v<CGAL::Constrained_Delaunay_triangulation_2<Epick>::Intersection_tag,
                               CGAL::No_constraint_intersection_requiring_constructions_tag>); // default tag

  try
  {
    CDT t;
    t.insert_constraint(Point(0,0), Point(2,2));
    t.insert_constraint(Point(0,2), Point(2,0));
    failure("no exception for a crossing");
  }
  catch(const CDT::Intersection_of_constraints_exception& e)
  {
    if(std::string(e.what()).empty())
      failure("empty what()");
  }
}

int main()
{
  test_exception();

  Intersection_tag_test<Epick> epick;
  Intersection_tag_test<Epeck> epeck;

  std::cout << "No_constraint_intersection_tag" << std::endl;
  epick.test_No_constraint_intersection_tag();
  epeck.test_No_constraint_intersection_tag();

  std::cout << "No_constraint_intersection_requiring_constructions_tag" << std::endl;
  epick.test_No_constraint_intersection_requiring_constructions_tag();
  epeck.test_No_constraint_intersection_requiring_constructions_tag();

#ifndef CGAL_NO_DEPRECATED_CODE
  std::cout << "No_intersection_tag" << std::endl;
  epick.test_No_intersection_tag();
#endif

  std::cout << "Exact_predicates_tag" << std::endl;
  epick.test_Exact_predicates_tag();
  epeck.test_Exact_predicates_tag();

  std::cout << "Exact_intersections_tag" << std::endl;
  epeck.test_Exact_intersections_tag();

  if(nb_failures > 0) {
    std::cerr << nb_failures << " failures" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}
