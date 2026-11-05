#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

#include <sstream>
#include <set>
#include <CGAL/intersections.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Tr = CGAL::Triangulation_3<K>;
using P = K::Point_3;

const std::array<P, 2> segment = {
    P{-20.304600000000001, 30.456900000000001, 0},
    P{402.2842, 41.8782, 0},
};

const std::string triangulation_content =
R"_tr(2
18
812.18299999999999 152.28399999999999 0
-20.304600000000001 30.456900000000001 0
395.93919999999997 91.370449999999991 0
662.43700000000001 365.48200000000003 0
824.87300000000005 53.299500000000002 0
819.79700000000003 -119.289 0
402.2842 41.8782 0
-116.751 -441.62400000000002 0
399.74619999999999 -44.416049999999998 0
-279.18799999999999 -418.78199999999998 0
-784.26400000000001 109.137 0
-710.65999999999997 225.88800000000001 0
-402.28430000000003 69.796949999999995 0
11 28 0
187.81729999999999 60.913674999999998 0
189.7208 -6.9795749999999988 0
84.708100000000002 11.7386625 0
32.201750000000004 21.097781250000001 0
34
13 10 2
5 0 1
8 0 6
10 8 2
4 3 1
4 1 0
16 8 9
15 2 14
7 5 1
7 1 3
6 0 5
9 5 7
10 0 8
17 16 15
9 8 6
9 6 5
11 0 10
13 12 11
13 11 10
12 2 4
12 4 0
12 0 11
13 2 12
16 7 15
18 14 2
18 2 8
16 9 7
15 7 3
15 3 4
15 4 2
17 8 16
18 8 17
18 17 15
18 15 14
3 22 18
5 8 10
10 14 12
25 0 12
9 5 28
1 20 4
14 26 30
24 33 29
1 9 11
4 27 8
1 15 2
8 26 15
2 3 16
23 32 30
2 15 6
10 11 14
12 18 21
21 18 22
16 0 17
29 20 22
5 21 19
16 17 20
19 17 0
27 13 26
7 25 33
3 31 24
11 23 6
9 28 23
4 29 27
19 7 28
6 13 31
30 32 25
13 33 31
7 24 32
)_tr";

void check(const Tr& tr)
{
  assert(tr.is_valid());
  for(auto a = tr.finite_vertices_begin(); a != tr.finite_vertices_end(); ++a)
    for(auto b = tr.finite_vertices_begin(); b != tr.finite_vertices_end(); ++b)
    {
      if(a == b) continue;
      Tr::Cell_handle last;
      std::set<Tr::Cell_handle> seen;
      for(auto it = tr.segment_traverser_cells_begin(a, b);
          it != tr.segment_traverser_cells_end(); ++it)
      {
        last = it;
        assert(seen.insert(last).second);
        assert(seen.size() <= tr.number_of_facets());
        assert(!tr.is_infinite(last->vertex(0)) && !tr.is_infinite(last->vertex(1)) && !tr.is_infinite(last->vertex(2)));
        assert(CGAL::do_intersect(K::Segment_3(a->point(), b->point()),
          K::Triangle_3(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point())));
      }
      assert(last != Tr::Cell_handle() && last->has_vertex(b));
    }
  std::vector<P> points;
  for(const auto& v : tr.finite_vertex_handles())
    points.push_back(v->point());
  for(const auto& facet : tr.finite_facets())
  {
    const auto cell = facet.first;
    points.push_back(CGAL::centroid(cell->vertex(0)->point(), cell->vertex(1)->point(), cell->vertex(2)->point()));
    for(int i = 0; i < 3; ++i)
      points.push_back(CGAL::midpoint(cell->vertex(i)->point(), cell->vertex((i + 1) % 3)->point()));
  }
  assert(points.size() > tr.number_of_vertices());
  for(const P& a : points)
    for(const P& b : points)
    {
      if(a == b) continue;
      Tr::Cell_handle last;
      std::set<Tr::Cell_handle> seen;
      for(auto it = tr.segment_traverser_cells_begin(a, b);
          it != tr.segment_traverser_cells_end(); ++it)
      {
        last = it;
        assert(seen.insert(last).second);
        assert(seen.size() <= tr.number_of_facets());
        assert(CGAL::do_intersect(K::Segment_3(a, b),
          K::Triangle_3(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point())));
      }
      assert(last != Tr::Cell_handle());
      assert(K::Triangle_3(last->vertex(0)->point(), last->vertex(1)->point(), last->vertex(2)->point()).has_on(b));
    }

}

int main()
{
  Tr tr;
  std::istringstream input(triangulation_content);
  input >> tr;
  check(tr);
  for(int plane = 0; plane < 3; ++plane)
  {
    Tr other;
    for(const auto& v : tr.finite_vertex_handles())
    {
      const auto& p = v->point();
      other.insert(plane == 0 ? P(0, p.x(), p.y()) :
                   plane == 1 ? P(p.x(), 0, p.y()) : P(-p.x(), p.y(), 0));
    }
    check(other);
  }
}
