#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/intersections_d.h>
#include <CGAL/double.h>
#include <CGAL/test_macros.h>
#include <CGAL/use.h>

#ifdef CGAL_USE_LEDA

#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
typedef leda_integer RT;
typedef leda_real    FT;

#elif defined CGAL_USE_GMP

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpz RT;
typedef CGAL::Gmpq FT;

#else

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
typedef CGAL::MP_Float     RT;
typedef CGAL::Quotient<RT> FT;

#endif

int main()
{ CGAL::IO::set_pretty_mode ( std::cerr );
  CGAL_TEST_START;
{
  typedef CGAL::Homogeneous_d<RT> Kernel;
  typedef CGAL::Point_d<Kernel>      Point;
  typedef CGAL::Vector_d<Kernel>     Vector;
  typedef CGAL::Direction_d<Kernel>  Direction;
  typedef CGAL::Hyperplane_d<Kernel> Hyperplane;
  typedef CGAL::Sphere_d<Kernel>     Sphere;
  CGAL_USE_TYPE(Sphere);
  typedef CGAL::Segment_d<Kernel>    Segment;
  typedef CGAL::Ray_d<Kernel>        Ray;
  typedef CGAL::Line_d<Kernel>       Line;

  {
    Direction dir(1,1,1), dir2(0,0,1);
    Point o(-1,-1,-1,1), org(3);
    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), p4(1,1,1,1);
    std::vector<Point> V = make_vector(p1,p2,p3);
    Hyperplane h1(p1,dir2);
    Hyperplane h2(V.begin(),V.end(),o,CGAL::ON_NEGATIVE_SIDE);

    /* intersections in 3-space */
    Point ip; Line il; Ray ir; Segment is;
    CGAL::Object res;
    Line l(org, p4);
    CGAL_TEST(CGAL::do_intersect(h2,l)){}
    CGAL_TEST(CGAL::do_intersect(l,h2)){}
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && l.has_on(ip)){}
    res = CGAL::intersection(l,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && l.has_on(ip)){}

    l = Line(org, p1);
    CGAL_TEST(CGAL::do_intersect(h2,l)){}
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(ip,res) && ip==p1){}

    l = Line(p1, p2);
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(il,res)&& il==l){}
    l = Line(p1, p2) + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(l,h2)){}

    Ray r(org, dir);
    CGAL_TEST(CGAL::do_intersect(h2,r)){}
    CGAL_TEST(CGAL::do_intersect(r,h2)){}
    res = CGAL::intersection(h2,r);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && r.has_on(ip)){}
    res = CGAL::intersection(r,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && r.has_on(ip)){}

    r = Ray(p4,dir);
    CGAL_TEST(!CGAL::do_intersect(r,h2)){}
    r = Ray(p1, Direction(-1,1,0));
    res = CGAL::intersection(h2,r);
    CGAL_TEST(CGAL::assign(ir,res) && ir==r){}
    r = r + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(r,h2)){}

    Segment s(org, p4);
    CGAL_TEST(CGAL::do_intersect(h2,s)){}
    CGAL_TEST(CGAL::do_intersect(s,h2)){}
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}
    res = CGAL::intersection(s,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}

    s = Segment(org, Point(-1,-1,-1,1));
    CGAL_TEST(!CGAL::do_intersect(h2,s)){}
    s = Segment(p1,p1);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && ip==p1){}
    s = Segment(org,p4);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}
    s = Segment(p1,p2);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(is,s)){}
    s = s + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(h2,s)){}
  }

  {
    Direction dir(1,1,1), dir2(0,0,1);
    Point o(-1,-1,-1,1), org(3);
    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), p4(1,1,1,1);
    Vector v(1,0,0,1);
    Point ip; Line il; Ray ir; Segment is;
    CGAL::Object res;

    // PURE LINES
    Line l1(p1,p2), l2(p2,p3);
    CGAL_TEST(CGAL::do_intersect(l1,l2)){}
    res = CGAL::intersection(l1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && l1.has_on(ip) && l2.has_on(ip)){}
    CGAL_TEST(!CGAL::do_intersect(l1,l2+v)){}
    CGAL_TEST(!CGAL::do_intersect(l1,l1+v)){}
    CGAL_TEST(CGAL::do_intersect(l1,l1.opposite())){}
    res = CGAL::intersection(l1,l1.opposite());
    CGAL_TEST(CGAL::assign(il,res) && CGAL::weak_equality(il,l1)){}

    // PURE RAYS
    Ray r1(p2,p1), r2(p2,p3);
    CGAL_TEST(CGAL::do_intersect(r1,r2)){}
    res = CGAL::intersection(r1,r2);
    CGAL_TEST(CGAL::assign(ip,res) && r1.has_on(ip) && r2.has_on(ip)){}
    CGAL_TEST(!CGAL::do_intersect(r1,r2+(p3-p2))){}
    CGAL_TEST(!CGAL::do_intersect(r1,r1+v)){}
    CGAL_TEST(!CGAL::do_intersect(r1,r2+v)){}
    CGAL_TEST(CGAL::do_intersect(r1,r2+(p2-p3))){}
    CGAL_TEST(CGAL::do_intersect(r1,r1.opposite())){}
    CGAL_TEST(CGAL::do_intersect(r1+-(r1.direction().vector()),
                                 r1.opposite())){}
    res = CGAL::intersection(r1,r1.opposite());
    CGAL_TEST(CGAL::assign(ip,res) && ip == r1.source()){}
    res = CGAL::intersection(r1+-(r1.direction().vector()),r1.opposite());
    CGAL_TEST(CGAL::assign(is,res)){}

    // PURE SEGMENTS
    Segment s1(p2,p1), s2(p2,p3), s3;
    CGAL_TEST(CGAL::do_intersect(s1,s2)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+s2.vector())){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+-2*s2.vector())){}
    CGAL_TEST(!CGAL::do_intersect(s1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+v)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s1+2*s1.vector())){}
    res = CGAL::intersection(s1,s2);
    CGAL_TEST(CGAL::assign(ip,res) &&
              s1.source()==ip && s2.source()==ip){}
    res = CGAL::intersection(s1,s1+(s1.vector()/2));
    CGAL_TEST(CGAL::assign(is,res) &&
              s1.has_on(is.source()) && s1.has_on(is.target())){}
    res = CGAL::intersection(s1,s1+s1.vector());
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip)){}
    s3 = s1 + s1.vector()/2;
    res = CGAL::intersection(s1,s3);
    CGAL_TEST(CGAL::assign(is,res) &&
              s1.has_on(is.source()) && s1.has_on(is.target()) &&
              s3.has_on(is.source()) && s3.has_on(is.target())){}

    // MIXED LINE RAY
    CGAL_TEST(CGAL::do_intersect(l1,r1)){}
    CGAL_TEST(CGAL::do_intersect(r1,l1)){}
    CGAL_TEST(!CGAL::do_intersect(l1,r1+v)){}
    CGAL_TEST(!CGAL::do_intersect(l2,r1+r1.direction().vector())){}
    res = CGAL::intersection(r1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && r1.has_on(ip) && l2.has_on(ip)){}
    res = CGAL::intersection(l1,r1);
    CGAL_TEST(CGAL::assign(ir,res) && r1==ir){}

    // MIXED LINE SEGMENT
    CGAL_TEST(CGAL::do_intersect(l1,s1)){}
    CGAL_TEST(CGAL::do_intersect(s1,l1)){}
    CGAL_TEST(!CGAL::do_intersect(l1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(l2,s1+s1.vector())){}
    res = CGAL::intersection(s1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip) && l2.has_on(ip)){}
    res = CGAL::intersection(l1,s1);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(s1,is)){}

    // MIXED RAY SEGMENT
    CGAL_TEST(CGAL::do_intersect(r1,s1)){}
    CGAL_TEST(CGAL::do_intersect(s1,r1)){}
    CGAL_TEST(!CGAL::do_intersect(r1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(r2,s1+s1.vector())){}
    res = CGAL::intersection(s1,r2);
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip) && r2.has_on(ip)){}
    res = CGAL::intersection(r1,s1);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(s1,is)){}
    s3 = s1 + -s1.vector()/2;
    res = CGAL::intersection(s3,r1);
    CGAL_TEST(CGAL::assign(is,res) && is.source()==r1.source() &&
              is.target()==s3.target()){}

  }
}
{
  typedef CGAL::Cartesian_d<FT>  Kernel;
  typedef CGAL::Point_d<Kernel>      Point;
  typedef CGAL::Vector_d<Kernel>     Vector;
  typedef CGAL::Direction_d<Kernel>  Direction;
  typedef CGAL::Hyperplane_d<Kernel> Hyperplane;
  typedef CGAL::Sphere_d<Kernel>     Sphere;
  CGAL_USE_TYPE(Sphere);
  typedef CGAL::Segment_d<Kernel>    Segment;
  typedef CGAL::Ray_d<Kernel>        Ray;
  typedef CGAL::Line_d<Kernel>       Line;

  {
    Direction dir(1,1,1), dir2(0,0,1);
    Point o(-1,-1,-1,1), org(3);
    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), p4(1,1,1,1);
    std::vector<Point> V = make_vector(p1,p2,p3);
    Hyperplane h1(p1,dir2);
    Hyperplane h2(V.begin(),V.end(),o,CGAL::ON_NEGATIVE_SIDE);

    /* intersections in 3-space */
    Point ip; Line il; Ray ir; Segment is;
    CGAL::Object res;
    Line l(org, p4);
    CGAL_TEST(CGAL::do_intersect(h2,l)){}
    CGAL_TEST(CGAL::do_intersect(l,h2)){}
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && l.has_on(ip)){}
    res = CGAL::intersection(l,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && l.has_on(ip)){}

    l = Line(org, p1);
    CGAL_TEST(CGAL::do_intersect(h2,l)){}
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(ip,res) && ip==p1){}

    l = Line(p1, p2);
    res = CGAL::intersection(h2,l);
    CGAL_TEST(CGAL::assign(il,res)&& il==l){}
    l = Line(p1, p2) + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(l,h2)){}

    Ray r(org, dir);
    CGAL_TEST(CGAL::do_intersect(h2,r)){}
    CGAL_TEST(CGAL::do_intersect(r,h2)){}
    res = CGAL::intersection(h2,r);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && r.has_on(ip)){}
    res = CGAL::intersection(r,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && r.has_on(ip)){}

    r = Ray(p4,dir);
    CGAL_TEST(!CGAL::do_intersect(r,h2)){}
    r = Ray(p1, Direction(-1,1,0));
    res = CGAL::intersection(h2,r);
    CGAL_TEST(CGAL::assign(ir,res) && ir==r){}
    r = r + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(r,h2)){}

    Segment s(org, p4);
    CGAL_TEST(CGAL::do_intersect(h2,s)){}
    CGAL_TEST(CGAL::do_intersect(s,h2)){}
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}
    res = CGAL::intersection(s,h2);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}

    s = Segment(org, Point(-1,-1,-1,1));
    CGAL_TEST(!CGAL::do_intersect(h2,s)){}
    s = Segment(p1,p1);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && ip==p1){}
    s = Segment(org,p4);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(ip,res) && h2.has_on(ip) && s.has_on(ip)){}
    s = Segment(p1,p2);
    res = CGAL::intersection(h2,s);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(is,s)){}
    s = s + Vector(1,1,1,1);
    CGAL_TEST(!CGAL::do_intersect(h2,s)){}
  }

  {
    Direction dir(1,1,1), dir2(0,0,1);
    Point o(-1,-1,-1,1), org(3);
    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), p4(1,1,1,1);
    Vector v(1,0,0,1);
    Point ip; Line il; Ray ir; Segment is;
    CGAL::Object res;

    // PURE LINES
    Line l1(p1,p2), l2(p2,p3);
    CGAL_TEST(CGAL::do_intersect(l1,l2)){}
    res = CGAL::intersection(l1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && l1.has_on(ip) && l2.has_on(ip)){}
    CGAL_TEST(!CGAL::do_intersect(l1,l2+v)){}
    CGAL_TEST(!CGAL::do_intersect(l1,l1+v)){}
    CGAL_TEST(CGAL::do_intersect(l1,l1.opposite())){}
    res = CGAL::intersection(l1,l1.opposite());
    CGAL_TEST(CGAL::assign(il,res) && CGAL::weak_equality(il,l1)){}

    // PURE RAYS
    Ray r1(p2,p1), r2(p2,p3);
    CGAL_TEST(CGAL::do_intersect(r1,r2)){}
    res = CGAL::intersection(r1,r2);
    CGAL_TEST(CGAL::assign(ip,res) && r1.has_on(ip) && r2.has_on(ip)){}
    CGAL_TEST(!CGAL::do_intersect(r1,r2+(p3-p2))){}
    CGAL_TEST(!CGAL::do_intersect(r1,r1+v)){}
    CGAL_TEST(!CGAL::do_intersect(r1,r2+v)){}
    CGAL_TEST(CGAL::do_intersect(r1,r2+(p2-p3))){}
    CGAL_TEST(CGAL::do_intersect(r1,r1.opposite())){}
    CGAL_TEST(CGAL::do_intersect(r1+-(r1.direction().vector()),
                                 r1.opposite())){}
    res = CGAL::intersection(r1,r1.opposite());
    CGAL_TEST(CGAL::assign(ip,res) && ip == r1.source()){}
    res = CGAL::intersection(r1+-(r1.direction().vector()),r1.opposite());
    CGAL_TEST(CGAL::assign(is,res)){}

    // PURE SEGMENTS
    Segment s1(p2,p1), s2(p2,p3), s3;
    CGAL_TEST(CGAL::do_intersect(s1,s2)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+s2.vector())){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+-2*s2.vector())){}
    CGAL_TEST(!CGAL::do_intersect(s1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s2+v)){}
    CGAL_TEST(!CGAL::do_intersect(s1,s1+2*s1.vector())){}
    res = CGAL::intersection(s1,s2);
    CGAL_TEST(CGAL::assign(ip,res) &&
              s1.source()==ip && s2.source()==ip){}
    res = CGAL::intersection(s1,s1+(s1.vector()/2));
    CGAL_TEST(CGAL::assign(is,res) &&
              s1.has_on(is.source()) && s1.has_on(is.target())){}
    res = CGAL::intersection(s1,s1+s1.vector());
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip)){}
    s3 = s1 + s1.vector()/2;
    res = CGAL::intersection(s1,s3);
    CGAL_TEST(CGAL::assign(is,res) &&
              s1.has_on(is.source()) && s1.has_on(is.target()) &&
              s3.has_on(is.source()) && s3.has_on(is.target())){}

    // MIXED LINE RAY
    CGAL_TEST(CGAL::do_intersect(l1,r1)){}
    CGAL_TEST(CGAL::do_intersect(r1,l1)){}
    CGAL_TEST(!CGAL::do_intersect(l1,r1+v)){}
    CGAL_TEST(!CGAL::do_intersect(l2,r1+r1.direction().vector())){}
    res = CGAL::intersection(r1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && r1.has_on(ip) && l2.has_on(ip)){}
    res = CGAL::intersection(l1,r1);
    CGAL_TEST(CGAL::assign(ir,res) && r1==ir){}

    // MIXED LINE SEGMENT
    CGAL_TEST(CGAL::do_intersect(l1,s1)){}
    CGAL_TEST(CGAL::do_intersect(s1,l1)){}
    CGAL_TEST(!CGAL::do_intersect(l1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(l2,s1+s1.vector())){}
    res = CGAL::intersection(s1,l2);
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip) && l2.has_on(ip)){}
    res = CGAL::intersection(l1,s1);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(s1,is)){}

    // MIXED RAY SEGMENT
    CGAL_TEST(CGAL::do_intersect(r1,s1)){}
    CGAL_TEST(CGAL::do_intersect(s1,r1)){}
    CGAL_TEST(!CGAL::do_intersect(r1,s1+v)){}
    CGAL_TEST(!CGAL::do_intersect(r2,s1+s1.vector())){}
    res = CGAL::intersection(s1,r2);
    CGAL_TEST(CGAL::assign(ip,res) && s1.has_on(ip) && r2.has_on(ip)){}
    res = CGAL::intersection(r1,s1);
    CGAL_TEST(CGAL::assign(is,res) && CGAL::weak_equality(s1,is)){}
    s3 = s1 + -s1.vector()/2;
    res = CGAL::intersection(s3,r1);
    CGAL_TEST(CGAL::assign(is,res) && is.source()==r1.source() &&
              is.target()==s3.target()){}

  }
}
  CGAL_TEST_END;
}

