// 3D intersection tests.

#include <CGAL/Object.h>
#include <CGAL/Point_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/MP_Float.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <cassert>

#include "create_bbox_mesh.h"

const double epsilon = 0.001;

struct randomint {
  randomint() ;
  int	get() const { return sequence[cur]; }
  int next() { cur = (cur+1)%11; return get();}
private:
  int sequence[11];
  int cur;
};

inline randomint::randomint()
{
  cur = 0;
  sequence[0] = 19;
  sequence[1] = 5;
  sequence[2] = 17;
  sequence[3] = 13;
  sequence[4] = 29;
  sequence[5] = 2;
  sequence[6] = 23;
  sequence[7] = 31;
  sequence[8] = 3;
  sequence[9] = 37;
  sequence[10] = 11;
}

randomint ri;

inline double to_nt(int d)
{
    return double(d);
}

template < typename K >
struct Test {

  typedef CGAL::Point_3< K >          P;
  typedef CGAL::Segment_3< K >        S;
  typedef CGAL::Line_3< K >           L;
  typedef CGAL::Plane_3< K >          Pl;
  typedef CGAL::Triangle_3< K >       Tr;
  typedef CGAL::Ray_3< K >            R;
  typedef CGAL::Iso_cuboid_3< K >     Cub;
  typedef CGAL::Bbox_3                Bbox;


  template < typename Type >
  bool approx_equal_nt(const Type &t1, const Type &t2)
  {
	if (t1 == t2)
		return true;
	if (CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < epsilon)
		return true;
	std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << "\n";
	return false;
  }

  template < typename Type >
  bool approx_equal(const Type&t1, const Type&t2)
  {
	return t1 == t2;
	// we need approx equal to check approx kernels, but maybe we should only test with exact kernels
	// (approx kernels were useful before, when the text output was checked by diff ?)
	// idea : test containment with intervals ?  or use some "epsilon double"?
	// I need to convert the text output to exact rationals in the source...
	// Well, for now the current scheme works.
  }

  bool approx_equal(const P & p, const P & q)
  {
	return approx_equal_nt(p.x(), q.x()) &&
	       approx_equal_nt(p.y(), q.y()) &&
	       approx_equal_nt(p.z(), q.z());
  }

  bool approx_equal(const S & p, const S & q)
  {
	return approx_equal(p.source(), q.source()) && approx_equal(p.target(), q.target());
  }

  /*
  bool approx_equal(const Pol & p, const Pol & q)
  {
	if (p.size() != q.size())
		return false;
	for(typename Pol::const_iterator itp = p.begin(), itq = q.begin(); itp != p.end(); ++itp, ++itq)
		if (!approx_equal(*itp, *itq))
			return false;
	return true;
  }
  */

  template < typename O1, typename O2>
  void check_no_intersection(const O1& o1, const O2& o2)
  {
	assert(!CGAL::do_intersect(o1, o2));
	assert(!CGAL::do_intersect(o2, o1));
  }

  template < typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2)
  {
	Res tmp;
	assert(CGAL::do_intersect(o1, o2));
	assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
	assert(CGAL::do_intersect(o2, o1));
	assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
  }

  template < typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2, const Res& result, bool do_opposite = true)
  {
	Res tmp;
	assert(CGAL::do_intersect(o1, o2));
	assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
	assert(approx_equal(tmp, result));
	if (do_opposite) {
	  assert(CGAL::do_intersect(o2, o1));
	  assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
	  assert(approx_equal(tmp, result));
	}
  }


  P p(int x, int y, int z)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(z*w), to_nt(w));
  }

  Pl pl(int a, int b, int c, int d)
  {
    int w = ri.next();
    return Pl(to_nt(a*w), to_nt(b*w), to_nt(c*w), to_nt(d*w));
  }

  void Cub_Cub()
  {
    std::cout << "Iso_cuboid - Iso_cuboid\n";
    check_intersection     (Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, -98725, 43), p(17, 9025473, 47)),
		            Cub(p(-7, 6, 43), p(17, 71, 47)));
    check_no_intersection  (Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, 98725, 43), p(17, 9025473, 47)));
    check_no_intersection  (Cub(p(-73, 6, 1), p(-70, 71, 58)), Cub(p(8, -98725, 43), p(17, 9025473, 47)));
    check_no_intersection  (Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, -98725, -47), p(17, 9025473, -43)));
  }

  void L_Cub()
  {
    std::cout << "Line - Iso_cuboid\n";
    check_intersection     (L(p(-3, 1,-5), p(  -2, -5, -7)), Cub(p(  -7,     -8, -9), p(-1,       2, -4)),
		            S(P(-3.16667, 2, -4.66667), P(-1.5, -8, -8)));
    check_intersection     (L(p( 0, 0, 3), p(   1,  2,  3)), Cub(p(   1,      1,  1), p( 3,       5,  8)),
		            S(P(       1, 2,        3), P( 2.5,  5,  3)));
    check_intersection     (L(p( 1, 0, 0), p(   1,  2,  3)), Cub(p(   1,      1,  1), p( 3,       5,  8)),
		            S(P(       1, 1,      1.5), P(   1,  5, 7.5)));
    check_intersection     (L(p( 0, 2, 0), p(   1,  2,  3)), Cub(p(   2,      1,  1), p( 3,       5,  8)),
		            S(P(       2, 2,        6), P(2.66667, 2, 8)));
    check_no_intersection  (L(p( 0, 0, 0), p(   1,  0,  3)), Cub(p(   2,      1,  1), p( 3,       5,  8)));
    check_no_intersection  (L(p( 4, 0, 0), p(   4,  1,  3)), Cub(p(   2,      1,  1), p( 3,       5,  8)));
    check_intersection     (L(p( 0, 0, 0), p(   1,  2,  3)), Cub(p(   1,      1,  1), p( 3,       5,  8)),
		            S(P(       1, 2,        3), P(2.5, 5, 7.5)));
  }

  void Pl_L()
  {
    std::cout << "Plane - Line\n";
    check_intersection     (pl(1, 1, 1, 0), L(p(   1,      1,  1), p( 2,       3,  4)),
		            P(0.5, 0, -0.5));
    check_intersection<L>  (pl(0, 0, 1,-1), L(p(   1,      1,  1), p( 2,       3,  1)));
    check_no_intersection  (pl(0, 0, 1,-2), L(p(   1,      1,  1), p( 2,       3,  1)));
    check_intersection     (pl(1, 0, 1, 3), L(p(   1,      1,  1), p( 2,       3, -1)),
		            P(  6, 11, -9));
    check_intersection     (pl(1, 2, 4, 7), L(p(   1,      1,  1), p( 2,       3,  4)),
		            P( 0.176471, -0.647059, -1.47059));
  }

  void Pl_Pl()
  {
    std::cout << "Plane - Plane\n";
    check_intersection     (pl(0, 0, 1, 0), pl(0, 1, 0, 0), L(P(0, 0, 0), P(-85, 0, 0)), false);
    check_intersection     (pl(0, 1, 0, 0), pl(0, 0, 1, 0), L(P(-85, 0, 0), P(0, 0, 0)), false);
    check_intersection<Pl> (pl(0, 0, 1, 1), pl(0, 0, 3, 3));
    check_no_intersection  (pl(2, 1, 3, 4), pl(6, 3, 9, 3));
    check_intersection<Pl> (pl(2, 1, 3, 4), pl(6, 3, 9, 12));
    check_intersection     (pl(2, 3, 7, 5), pl(9, 7, 1, 3), L(P(2,-3, 0), P(-3908, 5182, -1105)), false);
    check_intersection     (pl(9, 7, 1, 3), pl(2, 3, 7, 5), L(P(-3908, 5182, -1105), P(2,-3, 0)), false);
  }

  void Pl_Pl_Pl()
  {
    std::cout << "Plane - Plane - Plane\n";
    Pl pl1(1,0,0,0);
    Pl pl2(0,1,0,0);
    Pl pl3(0,0,1,0);

    // Generic intersection.
    CGAL::Object o = CGAL::intersection(pl1, pl2, pl3);
    P p;
    assert(assign(p, o));
    assert(p == P(0,0,0));

    // Empty intersection.
    Pl pl4(1,0,0,1); // pl4 is // to pl1.

    CGAL::Object o2 = CGAL::intersection(pl1, pl2, pl4);
    assert(o2.is_empty());

    CGAL::Object o3 = CGAL::intersection(pl1, pl4, pl2);
    assert(o3.is_empty());

    // Intersection in a line.
    Pl pl5(1,1,0,0); // pl1, pl2, pl5 intersect in the line l.
    L l;

    CGAL::Object o4 = CGAL::intersection(pl1, pl2, pl5);
    assert(assign(l, o4));

    assert(l == L(P(0,0,0), P(0,0,1)));

    // Intersection in a plane.
    CGAL::Object o5 = CGAL::intersection(pl1, pl1, pl1);
    Pl pl;
    assert(assign(pl, o5));
    assert(pl == pl1);
  }

  void Pl_R()
  {
    std::cout << "Plane - Ray\n";
    check_no_intersection  (pl( 1, 1, 1,  0), R(p(1, 1, 1), p(2, 3, 4)));
    check_intersection     (pl(-3, 7, 5, -2), R(p(12, -3, 7), p(-1, 5, -3)), P(5.06667, 1.26667, 1.66667));
    check_intersection<R>  (pl( 0, 0, 1, -1), R(p( 1,  1, 1), p( 2, 3,  1)));
    check_no_intersection  (pl( 0, 0, 1, -2), R(p( 1,  1, 1), p( 2, 3,  1)));
    check_intersection     (pl( 1, 0, 1,  3), R(p( 1,  1, 1), p( 2, 3, -1)), P(6, 11, -9));
    check_no_intersection  (pl( 1, 2, 4,  7), R(p( 1,  1, 1), p( 2, 3,  4)));
    check_intersection     (pl( 0, 0, 1,  0), R(p( 1,  1,-1), p( 2, 3,  4)), P(1.2, 1.4, 0));
    check_intersection     (pl( 0, 0, 1,  0), R(p( 2,  3, 4), p( 1, 1, -1)), P(1.2, 1.4, 0));
    check_intersection     (pl( 0, 0, 1,  0), R(p( 7,  1, 0), p(83, 1, -4)), P(  7,   1, 0));
    check_intersection     (pl( 0, 0, 1,  0), R(p(12,  6,-4), p( 7,25,  0)), P(  7,  25, 0));
  }

  void Pl_S()
  {
    std::cout << "Plane - Segment\n";
    check_no_intersection  (pl( 1, 1, 1,  0), S(p(1, 1, 1), p(2, 3, 4)));
    check_intersection     (pl(-3, 7, 5, -2), S(p(12, -3, 7), p(-1, 5, -3)), P(5.06667, 1.26667, 1.66667));
    check_intersection<S>  (pl( 0, 0, 1, -1), S(p( 1,  1, 1), p( 2, 3,  1)));
    check_no_intersection  (pl( 0, 0, 1, -2), S(p( 1,  1, 1), p( 2, 3,  1)));
    check_no_intersection  (pl( 1, 0, 1,  3), S(p( 1,  1, 1), p( 2, 3, -1)));
    check_no_intersection  (pl( 1, 2, 4,  7), S(p( 1,  1, 1), p( 2, 3,  4)));
    check_intersection     (pl( 0, 0, 1,  0), S(p( 1,  1,-1), p( 2, 3,  4)), P(1.2, 1.4, 0));
    check_intersection     (pl( 0, 0, 1,  0), S(p( 2,  3, 4), p( 1, 1, -1)), P(1.2, 1.4, 0));
    check_intersection     (pl( 0, 0, 1,  0), S(p( 7,  1, 0), p(83, 1, -4)), P(  7,   1, 0));
    check_intersection     (pl( 0, 0, 1,  0), S(p(12,  6,-4), p( 7,25,  0)), P(  7,  25, 0));
  }

  void R_Cub()
  {
    std::cout << "Ray - Iso_cuboid\n";
    check_intersection     (R(p( -3,  1,  -5), p( -2,  -5,  -7)), Cub(p( -7,  -8,  -9), p( -1,  2,  -4)),
		            S(P(-3, 1, -5), P(-1.5, -8, -8)));
    check_intersection     (R(p(  0,  0,   3), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            S(P( 1, 2,  3), P( 2.5,  5,  3)));
    check_intersection     (R(p(  1,  0,   0), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            S(P( 1, 1,1.5), P(   1,  5,7.5)));
    check_intersection     (R(p(  0,  2,   0), p(  1,   2,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)),
		            S(P( 2, 2,  6), P(2.66667,2, 8)));
    check_no_intersection  (R(p(  0,  0,   0), p(  1,   0,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)));
    check_no_intersection  (R(p(  4,  0,   0), p(  4,   1,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)));
    check_intersection     (R(p(  0,  0,   0), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            S(P( 1, 2,  3), P(2.5, 5, 7.5)));
  }

  void S_Cub()
  {
    std::cout << "Segment - Iso_cuboid\n";
    check_intersection     (S(p( -3,  1,  -5), p( -2,  -5,  -7)), Cub(p( -7,  -8,  -9), p( -1,  2,  -4)),
		            S(P(-3, 1, -5), P(  -2, -5, -7)));
    check_intersection     (S(p(  0,  0,   3), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            P( 1,  2,  3));
    check_intersection     (S(p(  1,  0,   0), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            S(P( 1, 1, 1.5), P(   1,  2, 3)));
    check_no_intersection  (S(p(  0,  2,   0), p(  1,   2,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)));
    check_no_intersection  (S(p(  0,  0,   0), p(  1,   0,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)));
    check_no_intersection  (S(p(  4,  0,   0), p(  4,   1,   3)), Cub(p(  2,   1,   1), p(  3,  5,   8)));
    check_intersection     (S(p(  0,  0,   0), p(  1,   2,   3)), Cub(p(  1,   1,   1), p(  3,  5,   8)),
		            P( 1, 2,  3));
  }
  
  void Pl_Tr()
  {
    std::cout << "Plane - Triangle\n";
    check_intersection     ( Pl(P(0,0,0),P(12,0,0),P(0,11,0)),Tr(P(0,0,0),P(1,0,0),P(0,1,0)),Tr(P(0,0,0),P(1,0,0),P(0,1,0)),true);
    check_intersection     ( Pl(P(0,0,0),P(12,0,0),P(0,11,0)),Tr(P(0,0,0),P(1,0,1),P(0,1,1)),P(0,0,0),true);
    check_intersection     ( Pl(P(0,0,0),P(12,0,0),P(0,11,0)),Tr(P(0,0,0),P(1,0,0),P(0,1,1)),S(P(0,0,0),P(1,0,0)),true);
    check_intersection     ( Pl(P(0,0,0),P(12,0,0),P(0,11,0)),Tr(P(1,0,-1),P(-1,0,-1),P(0,0,1)),S(P(0.5,0,0),P(-0.5,0,0)),true);
  }
  
  void S_L()
  {
    std::cout << "Segment - Line\n";
    check_intersection     ( S(P(0,0,0),P(12,0,0)),L(P(0,0,0),P(1,0,0)),S(P(0,0,0),P(12,0,0)),true);
    check_intersection     ( S(P(0,0,0),P(12,0,0)),L(P(0,1,0),P(0,-1,0)),P(0,0,0),true);
    check_intersection     ( S(P(-12,0,0),P(12,0,0)),L(P(0,1,0),P(0,-1,0)),P(0,0,0),true);
  }
  
  void R_L()
  {
    std::cout << "Ray - Line\n";
    check_intersection     (L(P(0,0,0),P(1,0,0)),R(P(3,0,0),P(6,0,0)),R(P(3,0,0),P(6,0,0)));
    check_no_intersection  (L(P(0,0,0),P(1,0,0)),R(P(3,0,1),P(6,0,1)));
    check_intersection     (L(P(0,0,0),P(1,0,0)),R(P(3,0,0),P(6,4,0)),P(3,0,0));
    check_intersection     (L(P(0,0,0),P(1,0,0)),R(P(0,-2,0),P(0,-1,0)),P(0,0,0));
    check_no_intersection  (L(P(0,0,0),P(1,0,0)),R(P(0, 2,0),P(0,4,0)));
    check_intersection     (L(P(0,0,0),P(1,0,0)),R(P(5,-2,0),P(5,-1,0)),P(5,0,0));
    check_no_intersection  (L(P(0,0,0),P(1,0,0)),R(P(6, 2,0),P(5,4,0)));
  }

  void R_S()
  {
    std::cout << "Ray - Segment\n";
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(0,-2,0),P(0,-1,0)),P(0,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,-2,0),P(-3,-1,0)),P(-3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(3,-2,0),P(3,-1,0)),P(3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(0,0,0),P(0,-1,0)),P(0,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,0,0),P(-3,-1,0)),P(-3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(3,0,0),P(3,-1,0)),P(3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(0,-2,0),P(0,0,0)),P(0,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,-2,0),P(-3,0,0)),P(-3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(3,-2,0),P(3,0,0)),P(3,0,0));
    check_no_intersection  (S(P(-3,0,0),P(3,0,0)),R(P(0,-1,0),P(0,-2,0)));
    check_no_intersection  (S(P(-3,0,0),P(3,0,0)),R(P(-3,-1,0),P(-3,-2,0)));
    check_no_intersection  (S(P(-3,0,0),P(3,0,0)),R(P(3,-1,0),P(3,-2,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-4,0,0),P(-2,0,0)),S(P(-3,0,0),P(3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-4,0,0),P(7,0,0)),S(P(-3,0,0),P(3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,0,0),P(-2,0,0)),S(P(-3,0,0),P(3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,0,0),P(7,0,0)),S(P(-3,0,0),P(3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(0,0,0),P(-2,0,0)),S(P(0,0,0),P(-3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(0,0,0),P(2,0,0)),S(P(0,0,0),P(3,0,0)));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(-3,0,0),P(-4,0,0)),P(-3,0,0));
    check_intersection     (S(P(-3,0,0),P(3,0,0)),R(P(3,0,0),P(4,0,0)),P(3,0,0));
    check_no_intersection  (S(P(-3,0,0),P(3,0,0)),R(P(4,0,0),P(5,0,0)));
    check_no_intersection  (S(P(-3,0,0),P(3,0,0)),R(P(-5,0,0),P(-8,0,0)));
  }

  void R_R()
  {
    std::cout << "Ray - Ray\n";
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(0,0,0),P(-1,0,0)),P(0,0,0));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(0,0,0),P(0,1,0)),P(0,0,0));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(5,0,0),P(5,1,0)),P(5,0,0));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(-1,0,0),P(3,0,0)),R(P(0,0,0),P(1,0,0)));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(2,0,0),P(3,0,0)),R(P(2,0,0),P(6,0,0)));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(2,0,0),P(-3,0,0)),S(P(0,0,0),P(2,0,0)),false);
    check_intersection     (R(P(2,0,0),P(-3,0,0)),R(P(0,0,0),P(1,0,0)),S(P(2,0,0),P(0,0,0)),false);
    check_no_intersection  (R(P(0,0,0),P(1,0,0)),R(P(-2,0,0),P(-3,0,0)));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(1,0,0),P(1,1,0)),P(1,0,0));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(1,-1,0),P(1,1,0)),P(1,0,0));
    check_intersection     (R(P(0,0,0),P(1,0,0)),R(P(1,-2,0),P(1,-1,0)),P(1,0,0));
    check_no_intersection  (R(P(0,0,0),P(1,0,0)),R(P(1,-1,0),P(1,-2,0)));
    check_no_intersection  (R(P(0,0,0),P(1,0,0)),R(P(0,-1,0),P(0,-2,0)));
    check_no_intersection  (R(P(0,0,0),P(1,0,0)),R(P(0,1,0),P(0,2,0)));
  }
  
  void Bbox_L(){
    std::cout << "Bbox - Line\n";
    Bbox box(-1,-1,-1,1,1,1);
    
    //not in x,y,z slab
    check_no_intersection  (box,L(P(2,0,0),P(2,0,0.5)));
    check_no_intersection  (box,L(P(0,2,0),P(0,2,0.5)));
    check_no_intersection  (box,L(P(0,0,2),P(0.5,0,2)));
    //not in x,y,z slab
    check_no_intersection  (box,L(P(2,0,0.5),P(2,0,0)));
    check_no_intersection  (box,L(P(0,2,0.5),P(0,2,0)));
    check_no_intersection  (box,L(P(0.5,0,2),P(0,0,2)));
    //in each slab, time not matching
    //xz
    check_no_intersection  (box,L(P(-8,0,0),P(0,0, 8)));
    check_no_intersection  (box,L(P( 8,0,0),P(0,0, 8))); 
    check_no_intersection  (box,L(P(-8,0,0),P(0,0,-8)));
    check_no_intersection  (box,L(P( 8,0,0),P(0,0,-8)));
    //yz
    check_no_intersection  (box,L(P(0,-8,0),P(0,0, 8)));
    check_no_intersection  (box,L(P(0, 8,0),P(0,0, 8))); 
    check_no_intersection  (box,L(P(0,-8,0),P(0,0,-8)));
    check_no_intersection  (box,L(P(0, 8,0),P(0,0,-8)));
    //xy
    check_no_intersection  (box,L(P(0,-8,0),P(8,0,0)));
    check_no_intersection  (box,L(P(0, 8,0),P(8,0,0))); 
    check_no_intersection  (box,L(P(0,-8,0),P(-8,0,0)));
    check_no_intersection  (box,L(P(0, 8,0),P(-8,0,0)));
    //Intersecting
    //xz
    check_intersection<S>  (box,L(P(-0.5,0,0),P(0,0, 0.5)));
    check_intersection<S>  (box,L(P( 0.5,0,0),P(0,0, 0.5))); 
    check_intersection<S>  (box,L(P(-0.5,0,0),P(0,0,-0.5)));
    check_intersection<S>  (box,L(P( 0.5,0,0),P(0,0,-0.5)));
    //yz
    check_intersection<S>  (box,L(P(0,-0.5,0),P(0,0, 0.5)));
    check_intersection<S>  (box,L(P(0, 0.5,0),P(0,0, 0.5))); 
    check_intersection<S>  (box,L(P(0,-0.5,0),P(0,0,-0.5)));
    check_intersection<S>  (box,L(P(0, 0.5,0),P(0,0,-0.5)));
    //xy
    check_intersection<S>  (box,L(P(0,-0.5,0),P(0.5,0,0)));
    check_intersection<S>  (box,L(P(0, 0.5,0),P(0.5,0,0))); 
    check_intersection<S>  (box,L(P(0,-0.5,0),P(-0.5,0,0)));
    check_intersection<S>  (box,L(P(0, 0.5,0),P(-0.5,0,0)));
  }
  
  void Bbox_R(){
    std::cout << "Bbox - Ray\n";
    Bbox box(-1,-1,-1,1,1,1);
    
    //not in x,y,z slab
    check_no_intersection  (box,R(P(2,0,0),P(2,0,0.5)));
    check_no_intersection  (box,R(P(0,2,0),P(0,2,0.5)));
    check_no_intersection  (box,R(P(0,0,2),P(0.5,0,2)));
    //not in x,y,z slab
    check_no_intersection  (box,R(P(2,0,0.5),P(2,0,0)));
    check_no_intersection  (box,R(P(0,2,0.5),P(0,2,0)));
    check_no_intersection  (box,R(P(0.5,0,2),P(0,0,2)));
    //in each slab, time not matching
    //xz
    check_no_intersection  (box,R(P(-8,0,0),P(0,0, 8)));
    check_no_intersection  (box,R(P( 8,0,0),P(0,0, 8))); 
    check_no_intersection  (box,R(P(-8,0,0),P(0,0,-8)));
    check_no_intersection  (box,R(P( 8,0,0),P(0,0,-8)));
    //yz
    check_no_intersection  (box,R(P(0,-8,0),P(0,0, 8)));
    check_no_intersection  (box,R(P(0, 8,0),P(0,0, 8))); 
    check_no_intersection  (box,R(P(0,-8,0),P(0,0,-8)));
    check_no_intersection  (box,R(P(0, 8,0),P(0,0,-8)));
    //xy
    check_no_intersection  (box,R(P(0,-8,0),P(8,0,0)));
    check_no_intersection  (box,R(P(0, 8,0),P(8,0,0))); 
    check_no_intersection  (box,R(P(0,-8,0),P(-8,0,0)));
    check_no_intersection  (box,R(P(0, 8,0),P(-8,0,0)));
    //Intersecting
    //xz
    check_intersection<S>  (box,R(P(-0.5,0,0),P(0,0, 0.5)));
    check_intersection<S>  (box,R(P( 0.5,0,0),P(0,0, 0.5))); 
    check_intersection<S>  (box,R(P(-0.5,0,0),P(0,0,-0.5)));
    check_intersection<S>  (box,R(P( 0.5,0,0),P(0,0,-0.5)));
    //yz
    check_intersection<S>  (box,R(P(0,-0.5,0),P(0,0, 0.5)));
    check_intersection<S>  (box,R(P(0, 0.5,0),P(0,0, 0.5))); 
    check_intersection<S>  (box,R(P(0,-0.5,0),P(0,0,-0.5)));
    check_intersection<S>  (box,R(P(0, 0.5,0),P(0,0,-0.5)));
    //xy
    check_intersection<S>  (box,R(P(0,-0.5,0),P(0.5,0,0)));
    check_intersection<S>  (box,R(P(0, 0.5,0),P(0.5,0,0))); 
    check_intersection<S>  (box,R(P(0,-0.5,0),P(-0.5,0,0)));
    check_intersection<S>  (box,R(P(0, 0.5,0),P(-0.5,0,0)));
  }

  void Bbox_Tr() {
    std::cout << "Bbox - Triangle\n";

    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    Bbox unit_bbox(-1., -1., -1.
                   ,1.,  1.,  1.);
    const Polyhedron unit_bbox_poly = create_bbox_mesh<Polyhedron>(unit_bbox);
    const Tree tree(faces(unit_bbox_poly).first,
                    faces(unit_bbox_poly).second,
                    unit_bbox_poly);

    const Tr tr(P(-3. ,  0. ,   0.),
                P(-2. ,  0.1,   0.),
                P(-0.5,  3. ,   0.));

    const bool b = CGAL::do_intersect(unit_bbox, tr);
    assert(b == false);
    assert(tree.do_intersect(tr) == b);

    CGAL::Random_points_in_cube_3<P> r(10.);

    std::size_t bbox_does_intersect_counter = 0;
    std::size_t plane_does_intersect_counter = 0;
    std::size_t do_intersect_counter = 0;
    std::cerr << "Begin random testing...\n"
              << "  (each 'o' in the following line is 1000 tests)\n  ";
#if __OPTIMIZE__
    const std::size_t nb_of_tests = 100000;
#else
    const std::size_t nb_of_tests = 10000;
#endif
    for(std::size_t i = 0, end = nb_of_tests; i < end; ++i)
    {
      if(i % 1000 == 0) std::cerr << "o";
      const P p0(*r++);
      const P p1(*r++);
      const P p2(*r++);
      const Tr tr(p0, p1, p2);

      const bool b = do_intersect(unit_bbox, tr);
      if(b) ++do_intersect_counter;
      const bool b_tree = tree.do_intersect(tr);
      if(b != b_tree) {
        std::stringstream err_msg;
        err_msg.precision(17);
        CGAL::set_pretty_mode(err_msg);
        err_msg << "do_intersect(\n"
                << "             " << unit_bbox << "\n,\n"
                << "             " << tr
                << "             ) = " << std::boolalpha << b << "\n"
                << "but the same test with AABB tree gives: "
                << std::boolalpha << b_tree << "\n";
        CGAL_error_msg(err_msg.str().c_str());
        std::exit(EXIT_FAILURE);
      }
      if(CGAL::do_overlap(unit_bbox, tr.bbox())) {
        ++bbox_does_intersect_counter;
        if(CGAL::do_intersect(unit_bbox, tr.supporting_plane())) {
          ++plane_does_intersect_counter;
        }
      }
      //      if();
    } // end for-loop
    std::cerr << "\n";
    std::cerr << "                      Number of tests: "
              << nb_of_tests << "\n";
    std::cerr << "        Number of bbox-does-intersect: "
              << bbox_does_intersect_counter << "\n";
    std::cerr << "Number of bbox-and-plane-do-intersect: "
              << plane_does_intersect_counter << "\n";
    std::cerr << "              Number of intersections: "
              << do_intersect_counter << "\n";
  } // end function Bbox_Tr

  void run()
  {
    std::cout << "3D Intersection tests\n";
    Cub_Cub();
    L_Cub();
    Pl_L();
    Pl_Pl();
    Pl_Pl_Pl();
    Pl_R();
    Pl_S();
    R_Cub();
    S_Cub();
    Pl_Tr();
    S_L();
    R_L();
    R_S();
    R_R();
    Bbox_L();
    Bbox_R();
    Bbox_Tr();
  }

};

int main()
{
	Test< CGAL::Cartesian<double>   >().run();
	Test< CGAL::Homogeneous<CGAL::MP_Float> >().run();
	// TODO : test more kernels.
}
