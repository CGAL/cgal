#ifndef CGAL_NO_CDT_2_WARNING
#define CGAL_NO_CDT_2_WARNING
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/Polygon_repair/Even_odd_rule.h>
#include <CGAL/version.h>

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using Kernel = CGAL::Simple_cartesian<CGAL::Exact_rational>;
using Point = Kernel::Point_2;
using Ring = std::vector<Point>;
using Polygon = CGAL::Polygon_2<Kernel>;
using PolygonWithHoles = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

enum class Location {
    Outside,
    Inside,
    Boundary
};

std::string qstr(const CGAL::Exact_rational& q) {
    std::ostringstream os;
    os << q;
    return os.str();
}

std::string pstr(const Point& p) {
    return "(" + qstr(p.x()) + ", " + qstr(p.y()) + ")";
}

std::string lstr(Location loc) {
    switch (loc) {
        case Location::Outside:
            return "OUTSIDE";
        case Location::Inside:
            return "INSIDE";
        case Location::Boundary:
            return "BOUNDARY";
    }
    return "UNKNOWN";
}

Point p(int x, int y) {
    return Point(CGAL::Exact_rational(x), CGAL::Exact_rational(y));
}

Point pq(int xn, int xd, int yn, int yd) {
    CGAL::Exact_rational x = CGAL::Exact_rational(xn) / CGAL::Exact_rational(xd);
    CGAL::Exact_rational y = CGAL::Exact_rational(yn) / CGAL::Exact_rational(yd);
    return Point(x, y);
}

Ring ring_from_ints(std::initializer_list<std::pair<int, int>> pts) {
    Ring out;
    out.reserve(pts.size());
    for (const auto& pt : pts) {
        out.push_back(p(pt.first, pt.second));
    }
    return out;
}

bool on_segment(const Point& a, const Point& b, const Point& q) {
    if (CGAL::orientation(a, b, q) != CGAL::COLLINEAR) {
        return false;
    }
    const CGAL::Exact_rational min_x = (std::min)(a.x(), b.x());
    const CGAL::Exact_rational max_x = (std::max)(a.x(), b.x());
    const CGAL::Exact_rational min_y = (std::min)(a.y(), b.y());
    const CGAL::Exact_rational max_y = (std::max)(a.y(), b.y());
    return (min_x <= q.x() && q.x() <= max_x && min_y <= q.y() && q.y() <= max_y);
}

struct WindingResult {
    bool boundary = false;
    int wn = 0;
};

WindingResult winding_number_ring(const Ring& ring, const Point& q) {
    WindingResult res;
    if (ring.size() < 3) {
        return res;
    }

    const size_t n = ring.size();
    for (size_t i = 0; i < n; ++i) {
        const Point& a = ring[i];
        const Point& b = ring[(i + 1) % n];

        if (on_segment(a, b, q)) {
            res.boundary = true;
            return res;
        }

        if (a.y() <= q.y()) {
            if (b.y() > q.y() && CGAL::orientation(a, b, q) == CGAL::LEFT_TURN) {
                ++res.wn;
            }
        } else {
            if (b.y() <= q.y() && CGAL::orientation(a, b, q) == CGAL::RIGHT_TURN) {
                --res.wn;
            }
        }
    }
    return res;
}

Location evenodd_from_input_winding(const std::vector<Ring>& rings, const Point& q) {
    int total_wn = 0;
    for (const auto& r : rings) {
        WindingResult wr = winding_number_ring(r, q);
        if (wr.boundary) {
            return Location::Boundary;
        }
        total_wn += wr.wn;
    }
    const int parity = ((total_wn % 2) + 2) % 2;
    return parity ? Location::Inside : Location::Outside;
}

Ring clean_ring(const Ring& ring) {
    Ring out;
    out.reserve(ring.size());
    for (const Point& pt : ring) {
        if (out.empty() || out.back() != pt) {
            out.push_back(pt);
        }
    }
    if (out.size() >= 2 && out.front() == out.back()) {
        out.pop_back();
    }
    return out;
}

Multipolygon repaired_evenodd(const std::vector<Ring>& rings) {
    Multipolygon in;
    for (const Ring& r0 : rings) {
        Ring r = clean_ring(r0);
        if (r.size() < 3) {
            continue;
        }
        Polygon poly(r.begin(), r.end());
        in.add_polygon_with_holes(PolygonWithHoles(poly));
    }
    return CGAL::Polygon_repair::repair(in, CGAL::Polygon_repair::Even_odd_rule());
}

Location locate_in_repaired(const Multipolygon& mp, const Point& q) {
    for (const auto& pwh : mp.polygons_with_holes()) {
        const CGAL::Bounded_side outer = pwh.outer_boundary().bounded_side(q);
        if (outer == CGAL::ON_BOUNDARY) {
            return Location::Boundary;
        }
        if (outer == CGAL::ON_UNBOUNDED_SIDE) {
            continue;
        }

        bool inside_hole = false;
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            const CGAL::Bounded_side hs = hit->bounded_side(q);
            if (hs == CGAL::ON_BOUNDARY) {
                return Location::Boundary;
            }
            if (hs == CGAL::ON_BOUNDED_SIDE) {
                inside_hole = true;
                break;
            }
        }
        if (!inside_hole) {
            return Location::Inside;
        }
    }
    return Location::Outside;
}

bool print_check(const std::string& label, const std::vector<Ring>& rings, const Multipolygon& mp, const Point& q) {
    const Location from_input = evenodd_from_input_winding(rings, q);
    const Location from_repair = locate_in_repaired(mp, q);
    const bool same = (from_input == from_repair);
    std::cout << "  " << label << "  p=" << pstr(q)
              << "  input_winding=" << lstr(from_input)
              << "  repair=" << lstr(from_repair)
              << "  => " << (same ? "MATCH" : "MISMATCH") << "\n";
    return same;
}

bool run_two_square_validation() {
    std::cout << "=== Validation: two overlapping squares (evenodd) ===\n";
    const std::vector<Ring> rings = {
        ring_from_ints({{0, 0}, {10, 0}, {10, 10}, {0, 10}}),
        ring_from_ints({{5, 0}, {15, 0}, {15, 10}, {5, 10}})
    };

    const Multipolygon repaired = repaired_evenodd(rings);

    bool ok = true;
    ok &= print_check("square witness A", rings, repaired, p(2, 5));  // inside first square only
    ok &= print_check("square witness B", rings, repaired, p(7, 5));  // in overlap, so outside under evenodd
    std::cout << '\n';
    return ok;
}

bool run_case45_repro() {
    std::cout << "=== Repro: case 45 witness points ===\n";
    std::cout << "Problem triangle reported in repaired output:\n";
    std::cout << "  T = {(4491,1100), (4499,1135), (4100,1105)}\n";

    // Exact ring bag that triggers the failure in CGAL Polygon_repair EVENODD.
    const std::vector<Ring> rings = {
        ring_from_ints({
            {4505, 401}, {4522, 1145}, {4503, 1162}, {2280, 1129}, {4100, 1105},
            {3360, 1050}, {3302, 1107}, {1026, 1126}, {1026, 235}
        }),
        ring_from_ints({
            {4100, 1105}, {4499, 1135}, {4491, 1100}
        }),
        ring_from_ints({
            {4491, 1100}, {4501, 1100}, {4501, 866}, {1146, 462}, {1071, 1067}, {4469, 1000}
        }),
        ring_from_ints({
            {3291, 1118}, {3360, 1050}, {4512, 1136}
        })
    };

    const Multipolygon repaired = repaired_evenodd(rings);

    bool controls_ok = true;
    controls_ok &= print_check("control OUT point", rings, repaired, p(500, 500));
    controls_ok &= print_check("control IN point ", rings, repaired, p(4200, 800));

    // Strictly interior point of the problem triangle: centroid.
    const Point tri_centroid = pq(13090, 3, 3340, 3);
    const bool tri_match = print_check("triangle centroid", rings, repaired, tri_centroid);

    std::cout << '\n';
    const bool reproduced = controls_ok && !tri_match;
    std::cout << "Repro verdict: " << (reproduced ? "REPRODUCED" : "NOT REPRODUCED") << "\n";
    return reproduced;
}

}  // namespace

int main() {
    std::cout << "CGAL_VERSION=" << CGAL_VERSION_STR << "\n";
    std::cout << "CGAL_VERSION_NR=" << CGAL_VERSION_NR << "\n\n";
    std::cout << "CGAL_RELEASE_DATE=" << CGAL_RELEASE_DATE << "\n\n";

    const bool validation_ok = run_two_square_validation();
    const bool repro_ok = run_case45_repro();

    if (!validation_ok) {
        std::cerr << "Validation failed: methodology check did not pass.\n";
        return 1;
    }
    if(!repro_ok) {
        std::cerr << "Mismatch as reported in Issue 9325 at triangle witness point was not observed.\n";
        return 0;
    }
    return 1;
}
