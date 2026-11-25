#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/intersections.h>

#include <iostream>
#include <vector>
#include <cstdlib>

#define TEST_ASSERT(cond) do { \
    if (!(cond)) { \
        std::cerr << "Assertion failed: " << #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        std::abort(); \
    } \
} while(0)


template <typename Tr>
void verify_traversal(const Tr& tr,
                      typename Tr::Vertex_handle va,
                      typename Tr::Vertex_handle vb)
{
    using Segment_3 = typename Tr::Segment;
    using Cell_handle = typename Tr::Cell_handle;
    using Kernel = typename Tr::Geom_traits;
    using Triangle_3 = typename Kernel::Triangle_3;
    using Point_3 = typename Kernel::Point_3;

    const Point_3& pa = va->point();
    const Point_3& pb = vb->point();
    Segment_3 seg(pa, pb);

    auto it = tr.segment_traverser_cells_begin(va, vb);
    auto end = tr.segment_traverser_cells_end();

    TEST_ASSERT(it != end);

    std::vector<Cell_handle> visited;
    for (; it != end; ++it) {
        visited.push_back(it);
    }

    TEST_ASSERT(!visited.empty());

    // 1. First cell must contain va
    int idx_va = -1;
    TEST_ASSERT(visited.front()->has_vertex(va, idx_va));

    // 2. Last cell must contain vb
    int idx_vb = -1;
    TEST_ASSERT(visited.back()->has_vertex(vb, idx_vb));

    // 3. Brute-force Oracle:
    // (a) Every visited finite cell must intersect seg in a Segment_3 of non-zero length.
    // Cells that only touch the segment at an isolated vertex or point yield Point_3 and must NOT be in visited.
    std::vector<Segment_3> isects;
    for (size_t i = 0; i < visited.size(); ++i) {
        Cell_handle c = visited[i];
        if (!c->has_vertex(tr.infinite_vertex())) {
            Triangle_3 tri(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point());
            auto isect = CGAL::intersection(tri, seg);
            TEST_ASSERT(isect.has_value());
            TEST_ASSERT(std::holds_alternative<Segment_3>(*isect));
            Segment_3 s = std::get<Segment_3>(*isect);
            TEST_ASSERT(s.source() != s.target());

            // Orient s from pa towards pb
            if (CGAL::squared_distance(pa, s.source()) > CGAL::squared_distance(pa, s.target())) {
                s = Segment_3(s.target(), s.source());
            }
            isects.push_back(s);
        }
    }
    TEST_ASSERT(!isects.empty());

    // (b) Monotonic chained coverage from pa to pb:
    // Tolerance 1e-10 accounts for floating-point construction roundoff in inexact kernels (EPICK),
    // whereas in exact kernels (EPECK) coordinates match exactly.
    // First segment must start at pa
    TEST_ASSERT(CGAL::squared_distance(isects.front().source(), pa) <= 1e-10);
    // Last segment must end at pb
    TEST_ASSERT(CGAL::squared_distance(isects.back().target(), pb) <= 1e-10);
    // Consecutive intersection segments must connect without gaps
    for (size_t i = 0; i + 1 < isects.size(); ++i) {
        double gap = CGAL::to_double(CGAL::squared_distance(isects[i].target(), isects[i+1].source()));
        TEST_ASSERT(gap <= 1e-10);
    }

    // (c) Topological adjacency: consecutive visited cells must share a simplex
    for (size_t i = 0; i + 1 < visited.size(); ++i) {
        Cell_handle c1 = visited[i];
        Cell_handle c2 = visited[i+1];
        bool share_simplex = false;
        for (int v1 = 0; v1 < 3; ++v1) {
            for (int v2 = 0; v2 < 3; ++v2) {
                if (!c1->has_vertex(tr.infinite_vertex()) && c1->vertex(v1) == c2->vertex(v2)) {
                    share_simplex = true;
                }
            }
        }
        TEST_ASSERT(share_simplex);
    }

    // (d) Completeness: every finite cell whose interior is crossed by seg MUST be visited
    for (auto c = tr.finite_cells_begin(); c != tr.finite_cells_end(); ++c) {
        Triangle_3 tri(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point());
        auto isect = CGAL::intersection(tri, seg);
        if (isect.has_value() && std::holds_alternative<Segment_3>(*isect)) {
            Segment_3 s = std::get<Segment_3>(*isect);
            if (s.source() != s.target()) {
                bool on_edge = false;
                for (int e = 0; e < 3; ++e) {
                    Segment_3 edge(c->vertex((e+1)%3)->point(), c->vertex((e+2)%3)->point());
                    if (CGAL::do_intersect(edge, s.source()) && CGAL::do_intersect(edge, s.target())) {
                        on_edge = true;
                        break;
                    }
                }
                if (!on_edge) {
                    bool found = false;
                    Cell_handle hc = c;
                    for (Cell_handle vc : visited) {
                        if (vc == hc) { found = true; break; }
                    }
                    TEST_ASSERT(found);
                }
            }
        }
    }
}

template <typename Tr>
void verify_traversal_points(const Tr& tr,
                             const typename Tr::Point& pa,
                             const typename Tr::Point& pb)
{
    using Segment_3 = typename Tr::Segment;
    using Cell_handle = typename Tr::Cell_handle;
    using Kernel = typename Tr::Geom_traits;
    using Triangle_3 = typename Kernel::Triangle_3;
    using Point_3 = typename Kernel::Point_3;

    Segment_3 seg(pa, pb);

    auto it = tr.segment_traverser_cells_begin(pa, pb);
    auto end = tr.segment_traverser_cells_end();

    TEST_ASSERT(it != end);

    std::vector<Cell_handle> visited;
    for (; it != end; ++it) {
        visited.push_back(it);
    }

    TEST_ASSERT(!visited.empty());

    // 1. Every finite cell visited must intersect the segment in a non-zero length Segment_3
    std::vector<Segment_3> isects;
    for (size_t i = 0; i < visited.size(); ++i) {
        Cell_handle c = visited[i];
        if (!c->has_vertex(tr.infinite_vertex())) {
            Triangle_3 tri(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point());
            auto isect = CGAL::intersection(tri, seg);
            TEST_ASSERT(isect.has_value());
            TEST_ASSERT(std::holds_alternative<Segment_3>(*isect));
            Segment_3 s = std::get<Segment_3>(*isect);
            TEST_ASSERT(s.source() != s.target());

            if (CGAL::squared_distance(pa, s.source()) > CGAL::squared_distance(pa, s.target())) {
                s = Segment_3(s.target(), s.source());
            }
            isects.push_back(s);
        } else {
            // Finite edge of infinite cell must intersect segment
            int inf;
            c->has_vertex(tr.infinite_vertex(), inf);
            Segment_3 finite_edge(c->vertex(tr.cw(inf))->point(), c->vertex(tr.ccw(inf))->point());
            TEST_ASSERT(CGAL::do_intersect(finite_edge, seg));
        }
    }

    // 2. Chaining across finite cells (tolerance 1e-10 accounts for EPICK floating-point roundoff)
    for (size_t i = 0; i + 1 < isects.size(); ++i) {
        double gap = CGAL::to_double(CGAL::squared_distance(isects[i].target(), isects[i+1].source()));
        TEST_ASSERT(gap <= 1e-10);
    }

    // 3. Consecutive cells must be topologically connected
    for (size_t i = 0; i + 1 < visited.size(); ++i) {
        Cell_handle c1 = visited[i];
        Cell_handle c2 = visited[i+1];
        bool share_simplex = false;
        for (int v1 = 0; v1 < 3; ++v1) {
            for (int v2 = 0; v2 < 3; ++v2) {
                if (c1->vertex(v1) == c2->vertex(v2)) {
                    share_simplex = true;
                }
            }
        }
        TEST_ASSERT(share_simplex);
    }
}

template <typename Kernel>
void test_grid_and_collinear()
{
    using Tr = CGAL::Triangulation_3<Kernel>;
    using Point_3 = typename Kernel::Point_3;

    Tr tr;
    // 4x4 grid in z = 0 with multiple collinear edges and triples
    std::vector<typename Tr::Vertex_handle> vhs;
    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            vhs.push_back(tr.insert(Point_3(x, y, 0)));
        }
    }
    TEST_ASSERT(tr.dimension() == 2);

    for (size_t i = 0; i < vhs.size(); ++i) {
        for (size_t j = i + 1; j < vhs.size(); ++j) {
            verify_traversal(tr, vhs[i], vhs[j]);
            verify_traversal(tr, vhs[j], vhs[i]);
        }
    }
    std::cout << "Grid and collinear tests passed!" << std::endl;
}

template <typename Kernel>
void test_random_2d(unsigned int seed)
{
    using Tr = CGAL::Delaunay_triangulation_3<Kernel>;
    using Point_3 = typename Kernel::Point_3;

    std::cout << "Testing random 2D point sets (seed = " << seed << ")..." << std::endl;
    Tr tr;
    CGAL::Random rng(seed);
    std::vector<typename Tr::Vertex_handle> vhs;
    for (int i = 0; i < 20; ++i) {
        double x = rng.get_double(-100.0, 100.0);
        double y = rng.get_double(-100.0, 100.0);
        vhs.push_back(tr.insert(Point_3(x, y, 0)));
    }
    TEST_ASSERT(tr.dimension() == 2);

    for (size_t i = 0; i < vhs.size(); ++i) {
        for (size_t j = i + 1; j < vhs.size(); ++j) {
            verify_traversal(tr, vhs[i], vhs[j]);
            verify_traversal(tr, vhs[j], vhs[i]);
        }
    }
    std::cout << "Random 2D tests passed for seed " << seed << "!" << std::endl;
}

template <typename Kernel>
void test_arbitrary_points()
{
    using Tr = CGAL::Delaunay_triangulation_3<Kernel>;
    using Point_3 = typename Kernel::Point_3;

    Tr tr;
    tr.insert(Point_3(0, 0, 0));
    tr.insert(Point_3(4, 0, 0));
    tr.insert(Point_3(4, 4, 0));
    tr.insert(Point_3(0, 4, 0));
    tr.insert(Point_3(2, 2, 0));
    TEST_ASSERT(tr.dimension() == 2);

    // 1. Non-vertex points strictly inside cells
    verify_traversal_points(tr, Point_3(1, 0.5, 0), Point_3(3, 3.5, 0));
    verify_traversal_points(tr, Point_3(3, 3.5, 0), Point_3(1, 0.5, 0));

    // 2. Segment starting inside and ending outside convex hull (traverses to infinite cell)
    verify_traversal_points(tr, Point_3(2, 2, 0), Point_3(6, 2, 0));
    verify_traversal_points(tr, Point_3(2, 2, 0), Point_3(2, 6, 0));
    verify_traversal_points(tr, Point_3(2, 2, 0), Point_3(-2, 2, 0));
    verify_traversal_points(tr, Point_3(2, 2, 0), Point_3(2, -2, 0));

    // 3. Segment starting outside and ending outside (crossing the triangulation)
    verify_traversal_points(tr, Point_3(-1, 2, 0), Point_3(5, 2, 0));
    verify_traversal_points(tr, Point_3(5, 2, 0), Point_3(-1, 2, 0));

    // 4. Points on boundary edges
    verify_traversal_points(tr, Point_3(1, 0, 0), Point_3(3, 4, 0));

    std::cout << "Arbitrary points and infinite cells tests passed!" << std::endl;
}

// Reproducer from issue #7130
template <typename Kernel>
void test_issue_7130_reproducer()
{
    using Tr = CGAL::Triangulation_3<Kernel>;
    using P = typename Kernel::Point_3;

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

    Tr tr;
    std::istringstream iss(triangulation_content);
    iss >> tr;
    TEST_ASSERT(tr.is_valid());
    TEST_ASSERT(tr.dimension() == 2);

    typename Tr::Vertex_handle va, vb;
    for (auto v = tr.finite_vertices_begin(); v != tr.finite_vertices_end(); ++v) {
        double vx = CGAL::to_double(v->point().x());
        double vy = CGAL::to_double(v->point().y());
        if (vx < -20.0 && vx > -21.0 && vy > 30.0 && vy < 31.0) va = v;
        if (vx > 402.0 && vx < 403.0 && vy > 41.0 && vy < 42.0) vb = v;
    }
    TEST_ASSERT(va != typename Tr::Vertex_handle() && vb != typename Tr::Vertex_handle());
    const P& pa = va->point();
    const P& pb = vb->point();

    auto it = tr.segment_traverser_cells_begin(va, vb);
    auto end = tr.segment_traverser_cells_end();

    std::vector<typename Tr::Cell_handle> visited;
    for (; it != end; ++it) {
        visited.push_back(it);
    }

    TEST_ASSERT(visited.size() == 5);
    int idx = -1;
    TEST_ASSERT(visited.front()->has_vertex(va, idx));
    TEST_ASSERT(visited.back()->has_vertex(vb, idx));

    typename Kernel::Segment_3 seg(pa, pb);
    for (auto c : visited) {
        TEST_ASSERT(!c->has_vertex(tr.infinite_vertex()));
        typename Kernel::Triangle_3 tri(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point());
        auto isect = CGAL::intersection(tri, seg);
        TEST_ASSERT(isect.has_value());
        TEST_ASSERT(std::holds_alternative<typename Kernel::Segment_3>(*isect));
        typename Kernel::Segment_3 iseg = std::get<typename Kernel::Segment_3>(*isect);
        TEST_ASSERT(iseg.source() != iseg.target());
    }

    for (size_t i = 0; i + 1 < visited.size(); ++i) {
        typename Tr::Cell_handle c1 = visited[i];
        typename Tr::Cell_handle c2 = visited[i+1];
        bool are_neighbors = (c1->neighbor(0) == c2 || c1->neighbor(1) == c2 || c1->neighbor(2) == c2);
        TEST_ASSERT(are_neighbors);
    }
    std::cout << "Issue #7130 reproducer test passed!" << std::endl;
}

int main()
{
    std::cout << "=== Testing EPICK ===" << std::endl;
    test_issue_7130_reproducer<CGAL::Exact_predicates_inexact_constructions_kernel>();
    test_grid_and_collinear<CGAL::Exact_predicates_inexact_constructions_kernel>();
    test_random_2d<CGAL::Exact_predicates_inexact_constructions_kernel>(42);
    test_random_2d<CGAL::Exact_predicates_inexact_constructions_kernel>(1337);
    test_arbitrary_points<CGAL::Exact_predicates_inexact_constructions_kernel>();

    std::cout << "=== Testing EPECK ===" << std::endl;
    test_issue_7130_reproducer<CGAL::Exact_predicates_exact_constructions_kernel>();
    test_grid_and_collinear<CGAL::Exact_predicates_exact_constructions_kernel>();
    test_random_2d<CGAL::Exact_predicates_exact_constructions_kernel>(42);
    test_random_2d<CGAL::Exact_predicates_exact_constructions_kernel>(1337);
    test_arbitrary_points<CGAL::Exact_predicates_exact_constructions_kernel>();

    std::cout << "ALL 2D COMPREHENSIVE TESTS PASSED SUCCESSFULLY!" << std::endl;
    return 0;
}
