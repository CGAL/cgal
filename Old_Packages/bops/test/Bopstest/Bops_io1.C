#include "Bops_types.h"
#include "Bops_io.h"
#include <iostream>

using namespace CGAL;
using namespace std;

static
bool
read_point(TestPoint &pt)
{
    int x, y, w;
    CGAL_STD::cin >> x >> y >> w;
    if (!CGAL_STD::cin.good())
	return false;
    pt = TestPoint(TestNum(x), TestNum(y), TestNum(w));
    return true;
}


static
bool
read_pgn(TestPolygon &pgn)
{
    int n, i;
    CGAL_STD::cin >> n;
    if (!CGAL_STD::cin.good())
	return false;
    if (n < 3) {
	CGAL_STD::cin.clear(ios::failbit);
	return false;
    }
    vector<TestPoint> points(n);
    for (i=0; i<n; i++) {
	if (!read_point(points[i]))
	    return false;
    }
    pgn = TestPolygon(points.begin(), points.end());
    if (pgn.is_clockwise_oriented()) {
	pgn.reverse_orientation();
    }
    return true;
}

static
void
write_point(const TestPoint &pt)
{
    CGAL_STD::cout << to_double(pt.x()) << ' ' << to_double(pt.y());
}

static
void
write_segment(const TestSegment &seg)
{
    write_point(seg.source());
    CGAL_STD::cout << "  ";
    write_point(seg.target());
}

static
void
write_pgn(const TestPolygon &pgn)
{
    TestPolygon::Vertex_const_iterator cur(pgn.vertices_begin());
    CGAL_STD::cout << pgn.size() << '\n';
    while (cur != pgn.vertices_end()) {
	write_point(*cur);
	++cur;
	CGAL_STD::cout <<'\n';
    }
}

static
void
write_object(const Object &obj)
{
    TestPolygon pgn;
    ResultPolygon pgn2;
    if (assign(pgn, obj)) {
    	if (pgn.size()<3) {
    	    CGAL_STD::cerr << "Polygon with "<<pgn.size()<<" vertices.\n";
    	}
	if (!pgn.is_simple() || pgn.size()<3 ) {
	    CGAL_STD::cout << "Degenerate polygon\n";
	} else {
	    if (pgn.is_clockwise_oriented())
		CGAL_STD::cout << "Polygon hole\n";
	    else
		CGAL_STD::cout << "Polygon filled\n";
	}
	write_pgn(pgn);
	return;
    }
    if (assign(pgn2, obj)) {
	CGAL_STD::cout << "Polygon\n";
	write_pgn(pgn);
	return;
    }
    TestSegment seg;
    if (assign(seg, obj)) {
	CGAL_STD::cout << "Segment\n";
	write_segment(seg);
	return;
    }
    TestPoint pt;
    if (assign(pt, obj)) {
	CGAL_STD::cout << "Point\n";
	write_point(pt);
	return;
    }
    CGAL_STD::cout << "Unknown";
}


bool read_input(TestPolygon &pgn1, TestPolygon &pgn2)
{
    if (!read_pgn(pgn1)) {
	return false;
    }
    if (!read_pgn(pgn2)) {
	return false;
    }
    return true;
}

void write_result(const list<Object> & result)
{
    CGAL_STD::cout << result.size() << '\n';
    list<Object>::const_iterator cur = result.begin();
    while (cur != result.end() ) {
	write_object(*cur);
	CGAL_STD::cout << '\n';
	++cur;
    }
    CGAL_STD::cout << flush;
}
