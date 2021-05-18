#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/double.h>
#include <CGAL/test_macros.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
typedef leda_integer  RT_;
typedef leda_rational FT_;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>
typedef CGAL::Gmpz RT_;
typedef CGAL::Quotient<RT_> FT_;
#else
typedef double RT_;
typedef double FT_;
#endif
#endif

int main()
{ CGAL_KD_SETDTHREAD(2);
  CGAL::IO::set_pretty_mode ( std::cerr );
  CGAL_TEST_START;
{ // Homogeneous Kernel
  typedef CGAL::Homogeneous_d<RT_> Kernel;
  typedef Kernel::RT             RT;
  typedef Kernel::FT             FT;
  typedef Kernel::LA             LA;
  typedef CGAL::Point_d<Kernel>      Point;
  typedef CGAL::Vector_d<Kernel>     Vector;
  typedef CGAL::Direction_d<Kernel>  Direction;
  typedef CGAL::Hyperplane_d<Kernel> Hyperplane;
  typedef CGAL::Sphere_d<Kernel>     Sphere;
  typedef CGAL::Segment_d<Kernel>    Segment;
  typedef CGAL::Ray_d<Kernel>        Ray;
  typedef CGAL::Line_d<Kernel>       Line;
  typedef LA::Vector                 IVector;
  bool DOIO=true;

  {
    /* some construction test */
    int IV1[] = {1,0,1};
    int IV2[] = {0,1};
    IVector iv1(IV1,IV1+3), iv2(IV2,IV2+2);
    Point p0(2,CGAL::ORIGIN), p1(2,IV1,IV1+3),
          p2(2,iv2.begin(),iv2.end(),1),
          p3(2), p4(p1), p5(2,iv2.begin(),iv2.end(),1);
    CGAL_TEST(p0 == p3 && p1 == p4 && p2 == p5){} // op==
    CGAL_TEST(p0 != p1 && p0 != p2){} // op!=
    CGAL_TEST(p0 == CGAL::ORIGIN){}

    if (DOIO) CGAL_IO_TEST(p1,p4,CGAL::IO::ASCII); p4 = p1;

    CGAL_TEST(p4.dimension()==p1.dimension()){} // dimension()
    CGAL_TEST(p4.dimension()!=Point(3).dimension()){} // dimension()

    /* some input and access test */
    CGAL_TEST(p1.hx()==RT(1) && p1.hy()==RT(0) && p1.hw()!=RT(0)){}
    CGAL_TEST(p1.hx()==p1.homogeneous(0) &&
              p1.hy()==p1.homogeneous(1) &&
              p1.hw()==p1.homogeneous(2)){} // hx,hy,hw
    CGAL_TEST(p1.x()==FT(1) && p1.y()==FT(0)){} // x(),y()
    CGAL_TEST(p1.x()==p1[0] && p1.y()==p1[1]){} // x(),y()
    CGAL_TEST(p1.x()==p1.cartesian(0) && p1.y()==p1.cartesian(1)){} // x(),y()
    Point::Homogeneous_const_iterator hcit; int i;
    for (i=0, hcit = p1.homogeneous_begin();
         hcit != p1.homogeneous_end(); ++i, ++hcit)
      CGAL_TEST(*hcit == iv1[i]){}

    Point::Homogeneous_const_iterator hit;
    Point::Cartesian_const_iterator cit;

    Kernel::Cartesian_const_iterator_d cit_begin =
      Kernel().construct_cartesian_const_iterator_d_object()(p1);
    Kernel::Cartesian_const_iterator_d cit_end =
      Kernel().construct_cartesian_const_iterator_d_object()(p1, 1);

    for (i=0,hit=p1.homogeneous_begin(),cit=p1.cartesian_begin();
         i<p1.dimension(); ++hit,++cit,++i) {
      CGAL_TEST(p1.homogeneous(i)==*hit){}
      CGAL_TEST(p1.cartesian(i)==*cit){}
      CGAL_TEST(p1.cartesian(i)==*cit_begin){}
      cit_begin++;
    }
    CGAL_TEST(cit==p1.cartesian_end()){}
    CGAL_TEST(cit_begin==cit_end){}
    CGAL_TEST(p1.homogeneous(p1.dimension())==*hit){}
    CGAL_TEST(++hit == p1.homogeneous_end()){}

    /* the geometric operations and predicates */

    CGAL_TEST(CGAL::project_along_d_axis(
              CGAL::lift_to_paraboloid(p1)) == p1){} // lift and project
    Vector vi(7,8,1);
    CGAL_TEST(CGAL::midpoint(p1+vi,p1-vi)==p1){}
    CGAL_TEST(CGAL::squared_distance(p1,p1+vi)==vi.squared_length()){}

    std::vector<Point> A = make_vector(p0,p1,p2);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::POSITIVE){}

    p3 = Point(1,1,2);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_BOUNDED_SIDE){}
    p3 = Point(1,1,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_ORIENTED_BOUNDARY){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_BOUNDARY){}
    p3 = Point(2,2,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_NEGATIVE_SIDE){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_UNBOUNDED_SIDE){}

    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_UNBOUNDED_SIDE){}

    Point B[2] = { p1, p2 };
    CGAL_TEST(CGAL::affine_rank(B,B+2)==1){}
    CGAL_TEST(CGAL::contained_in_affine_hull(B,B+2,Point(1,1,2))){}
    CGAL_TEST(CGAL::contained_in_affine_hull(B,B+2,Point(1,1,2)+5*(p1-p2))){}
    CGAL_TEST(CGAL::contained_in_simplex(B,B+2,Point(1,1,2))){}
    CGAL_TEST(!CGAL::contained_in_simplex(B,B+2,Point(5,5,1))){}
    CGAL_TEST(!CGAL::contained_in_simplex(B,B+2,Point(1,1,2)+5*(p1-p2))){}
    CGAL_TEST(CGAL::contained_in_simplex(A.begin(),A.end(),Point(1,1,3))){}
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),Point(5,5,1))){}
    CGAL_TEST(CGAL::affinely_independent(B,B+2)){}
    CGAL_TEST(CGAL::affinely_independent(A.begin(),A.end())){}
    A.push_back(Point(1,1));
    CGAL_TEST(!CGAL::affinely_independent(A.begin(),A.end())){}
    CGAL_TEST(CGAL::compare_lexicographically(p0,p1)==CGAL::SMALLER){}
    CGAL_TEST(CGAL::lexicographically_smaller(p0,p1)){}
    CGAL_TEST(CGAL::lexicographically_smaller_or_equal(p0,p1)){}
    CGAL_TEST(CGAL::compare_lexicographically(p1,p0)==CGAL::LARGER){}
    CGAL_TEST(CGAL::compare_lexicographically(p1,p1)==CGAL::EQUAL){}
  }

  {
    /* some construction test */
    Point p0(3); // the origin
    Vector e1 = Vector(3,Vector::Base_vector(),0);
    int IV1[] = {0,1,0,1};
    int IV2[] = {0,0,1};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    // the first unit vector

    Point p1(p0 + e1),
          p2(3,iv1.begin(),iv1.end()),
          p3(3,iv2.begin(),iv2.end(),1),
          p4(3);
    CGAL_TEST(CGAL::compare_lexicographically(p0,p4)==CGAL::EQUAL){}
    CGAL_TEST(CGAL::compare_lexicographically(p0,p1)==CGAL::SMALLER){}
    CGAL_TEST(CGAL::compare_lexicographically(p3,p0)==CGAL::LARGER){}

    CGAL_TEST(p1.cartesian(0)==FT(1)){}
    CGAL_TEST(p1[1]==FT(0)){}
    CGAL_TEST(p1.homogeneous(0)==Kernel::RT(1)){}
    CGAL_TEST( (p1 - CGAL::ORIGIN) == e1){}
    CGAL_TEST(CGAL::squared_distance(Point(1,1,1,1),p0)==FT(3)){}

    CGAL_TEST(p1-p0==e1){}
    CGAL_TEST(p1+e1==Point(2,0,0,1)){}
    CGAL_TEST(p1-e1==p0){}

    Point p5(3);
    CGAL_TEST((p5+=e1)==p1){}
    CGAL_TEST((p5-=e1)==p0){}

    /* orientation, sphere position, simplex position */
    std::vector<Point> A = make_vector(p0,p1,p2,p4);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::ZERO){}
    A[3] = p3;
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::POSITIVE){}
    std::swap(A[2],A[3]);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::NEGATIVE){}

    p4 = Point(1,1,1,2);
    std::swap(A[2],A[3]);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_POSITIVE_SIDE){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_ORIENTED_BOUNDARY){}
    p4 = Point(2,2,2,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_NEGATIVE_SIDE){}

    std::swap(A[2],A[3]);
    p4 = Point(1,1,1,2);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_BOUNDED_SIDE){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_BOUNDARY){}
    p4 = Point(2,2,2,1);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_UNBOUNDED_SIDE){}
    std::swap(A[2],A[3]);

    p4 = Point(1,1,1,3);
    CGAL_TEST(CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}
    p4 = Point(1,1,1,2);
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}

    /* affine hull, independence and rank */
    std::vector<Point> B = make_vector(p1,p2,p3);
    p4 = CGAL::midpoint(p1,p2);
    CGAL_TEST(CGAL::contained_in_affine_hull(B.begin(),B.end(),p4)){}

    p4 = p0;
    CGAL_TEST(!CGAL::contained_in_affine_hull(B.begin(),B.end(),p0)){}

    A[0] = Point(3,3,3,1);
    CGAL_TEST(CGAL::affinely_independent(A.begin(),A.end())){}
    CGAL_TEST(CGAL::affine_rank(A.begin(),A.end())==3){}
    A[3] = CGAL::midpoint(CGAL::midpoint(A[0],A[1]),A[2]);
    CGAL_TEST(!CGAL::affinely_independent(A.begin(),A.end())){}
    // degenerate_simplex
    CGAL_TEST(CGAL::affine_rank(A.begin(),A.end())==2){}
    // full_simplex
  }


  {
    /* construction and access */
    int IV1[] = {1,2,3,1};
    int IV2[] = {3,2,1};
    IVector i1(IV1,IV1+4), i2(IV2,IV2+3);
    Vector a0(3),
           a1(3,i1.begin(),i1.end()),
           a2(3,i2.begin(),i2.end(),1),
           a3(6,4,2,2), a4(1,0,1), a5(1,1,6), a6(a1);

    CGAL_TEST((a1==a6 && a3==a2)){}
    CGAL_TEST((a0!=a1 && a0!=a4)){}

    if (DOIO) CGAL_IO_TEST(a1,a6,CGAL::IO::ASCII); a6 = a1;

    CGAL_TEST(a2.dimension()==a3.dimension()){}
    CGAL_TEST(a3.dimension()!=a4.dimension()){}

    CGAL_TEST(a4.hx()==RT(1) && a4.hy()==RT(0) && a4.hw()!=RT(0)){}
    CGAL_TEST(a4.x()==FT(1) && a4.y()==FT(0)){}
    CGAL_TEST(Kernel::make_FT(a3.homogeneous(1),a3.homogeneous(3))==FT(2)){}
    CGAL_TEST(a3.cartesian(1)==FT(2)){}
    CGAL_TEST(a3[0]==FT(3)){}
    CGAL_TEST(a3.direction()==Direction(3,2,1)){}

    Vector::Homogeneous_const_iterator hit;
    Vector::Cartesian_const_iterator cit; int i;
    for (i=0,hit=a1.homogeneous_begin(),cit=a1.cartesian_begin();
         i<a1.dimension(); ++hit,++cit,++i) {
      CGAL_TEST(a1.homogeneous(i)==*hit){}
      CGAL_TEST(a1.cartesian(i)==*cit){}
    }
    CGAL_TEST(cit==a1.cartesian_end()){}
    CGAL_TEST(a1.homogeneous(a1.dimension())==*hit){}
    CGAL_TEST(++hit == a1.homogeneous_end()){}

    /* arithmetical operations and compare */
    CGAL_TEST(Vector::cmp(a0,a1)==CGAL::SMALLER){}
    CGAL_TEST(Vector::cmp(a1,Vector(1,2,3,1))==CGAL::EQUAL){}
    CGAL_TEST(Vector::cmp(a1,Vector(-1,0,0,1))==CGAL::LARGER){}

    CGAL_TEST(2*a1==Vector(2,4,6,1)){}
    CGAL_TEST(Kernel::make_FT(2,3)*a4==Vector(2,0,3)){}

    CGAL_TEST(a3/2==Vector(3,2,1,2)){}
    CGAL_TEST(a3/Kernel::make_FT(2,3)==Vector(9,6,3,2)){}

    Vector a7(a1);
    a7*=2;
    a7*=Kernel::make_FT(3,2);
    CGAL_TEST(3*a1==a7){}
    a7/=Kernel::make_FT(3,4);
    a7/=4;
    CGAL_TEST(a1==a7){}

    CGAL_TEST(a4*a4==FT(1)){}
    CGAL_TEST(a4.squared_length()==FT(1)){}
    CGAL_TEST(a1+a2==Vector(4,4,4,1)){}
    CGAL_TEST(a1-a2==Vector(-2,0,2,1)){}

    Vector eins(1,1,1,1);
    a7 = a1;
    CGAL_TEST((a7+=eins)==a1+eins){}
    CGAL_TEST((a7-=eins)==a1){}
    CGAL_TEST(a1-2*a1==-a1){}
    CGAL_TEST(a0.is_zero()){}

    /* some linear algebra */

    std::vector<Vector> P12 = make_vector(a1,a2);
    CGAL_TEST(CGAL::contained_in_linear_hull(P12.begin(),P12.end(),a1+a2)){}
    Vector z_off = a1 + Vector(3,Vector::Base_vector(),2);
    CGAL_TEST(!CGAL::contained_in_linear_hull(P12.begin(),P12.end(),z_off)){}

    std::vector<Vector> NLI = make_vector(a1,a2,a1+a2),
                        LI = make_vector(a1,a2,z_off);
    CGAL_TEST(CGAL::linearly_independent(LI.begin(),LI.end())){}
    CGAL_TEST(!CGAL::linearly_independent(NLI.begin(),NLI.end())){}

    CGAL_TEST(CGAL::linear_rank(LI.begin(),LI.end())==3){}
    CGAL_TEST(CGAL::linear_rank(NLI.begin(),NLI.end())==2){}
    std::vector<Vector> LB(3);
    std::vector<Vector>::iterator last =
      CGAL::linear_base(NLI.begin(),NLI.end(),LB.begin());
    CGAL_TEST(CGAL::linearly_independent(LB.begin(),last)){}
  }


  {
    /* construction and access */
    int IV1[] = {1,5,-3,2};
    int IV2[] = {1,1,0};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    Direction d0(3), d1(3,iv1.begin(),iv1.end()),
              d2(3,iv2.begin(),iv2.end()),
              d31(1,1,1), d32(1,1),
              d4 = Direction(3,Direction::Base_direction(),2),
              d5(d0);
    CGAL_TEST(d0==d5){}
    CGAL_TEST(d2==Direction(2,2,0)){}
    CGAL_TEST(d31!=d32 && d2!=d1)

    if (DOIO) CGAL_IO_TEST(d1,d5,CGAL::IO::ASCII); d5=d1;
    CGAL_TEST(d2.dimension()==d1.dimension()){}
    CGAL_TEST(d31.dimension()!=d32.dimension()){}
    CGAL_TEST(d1.dx()==d1.delta(0)){}
    CGAL_TEST(d1.dy()==d1.delta(1)){}
    CGAL_TEST(d1.dz()==d1.delta(2)){}
    CGAL_TEST(d1.dy()==d1[1]){}
    CGAL_TEST(d4.delta(0)==RT(0) && d4.delta(2)==RT(1)){}
    CGAL_TEST(d4==-Direction(0,0,-1)){}
    CGAL_TEST(d4==Direction(0,0,-1).opposite()){}
    CGAL_TEST(Direction::cmp(d0,d2)==CGAL::SMALLER){}
    Direction::Delta_const_iterator it; int i;
    for (i=0,it=d1.deltas_begin(); i<d1.dimension(); ++it,++i) {
      CGAL_TEST(d1.delta(i)==*it){}
    }
    CGAL_TEST(it==d1.deltas_end()){}
  }

  {
    /* construction and access */
    int IV1[] = {1,2,3,4};
    int IV2[] = {-4,-3,-2};
    IVector vi1(IV1,IV1+4);
    IVector vi2(IV2,IV2+3);
    // two ivec inits

    Point p1(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),0)),
          p2(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),1)),
          p3(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),2));
    // one three point init

    Direction dir(1,1,1);
    Point o(-1,-1,-1,1), org(3);
    // one dir init
    Hyperplane h0(3,vi1.begin(),vi1.end());
    Hyperplane h1(3,vi2.begin(),vi2.end(),1);
    std::vector<Point> V = make_vector(p1,p2,p3);
    Hyperplane h2(V.begin(),V.end(),o,CGAL::ON_NEGATIVE_SIDE);
    V = make_vector(org,org,org);
    Hyperplane h3(V.begin(),V.end(),o,CGAL::ON_POSITIVE_SIDE);
    Hyperplane h4(p1,-dir);
    Hyperplane h5(1,2,3,4);
    Hyperplane h6(1,0,1);
    Hyperplane h7(h0);
    Hyperplane h8(3);

    CGAL_TEST((h5==h0 && h7==h0)){}
    CGAL_TEST((h1!=h0 && h7!=h6)){}
    CGAL_TEST((h2.dimension()==h5.dimension() &&
               h2.dimension()!=h6.dimension())){}
    CGAL_TEST((h5.coefficient(1)==RT(2) && h5[2]==RT(3))){}
    CGAL_TEST((vi1 == h0.coefficient_vector())){}
    CGAL_TEST(h2.orthogonal_vector()*(p2-p1) == FT(0) ){}
    CGAL_TEST(h2.orthogonal_direction()==dir){}
    CGAL_TEST(h2.value_at(o)<FT(0)){}
    CGAL_TEST(h2.value_at(p3)==FT(0)){}
    CGAL_TEST(h2.value_at(p3+dir.vector())>FT(0)){}
    CGAL_TEST(h2.oriented_side(o)==CGAL::ON_NEGATIVE_SIDE){}
    CGAL_TEST(h2.oriented_side(p3)==CGAL::ON_ORIENTED_BOUNDARY){}
    CGAL_TEST(h2.oriented_side(p3+dir.vector())==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(h3.oriented_side(o)==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(h4.has_on(p1)){}
    CGAL_TEST(h3.has_on(org)){}
    CGAL_TEST(!h5.has_on(p1)){}
    CGAL_TEST(h3.has_on_positive_side(o)){}

    if (DOIO) CGAL_IO_TEST(h0,h8,CGAL::IO::ASCII); h8=h0;

    Hyperplane::Coefficient_const_iterator it; int i;
    for (i=0,it=h5.coefficients_begin();
         i<=h5.dimension(); ++it,++i) {
      CGAL_TEST(h5.coefficient(i)==*it){}
    }
    CGAL_TEST(it==h5.coefficients_end()){}

    // all kind of |HyperplaneHd<integer>| inits

    /* compares and other operations */
    CGAL_TEST(Hyperplane::strong_cmp(h0,h1)>0){}
    CGAL_TEST(Hyperplane::strong_cmp(h1,h0)<0){}
    CGAL_TEST(Hyperplane::strong_cmp(h0,h5)==0){}
    CGAL_TEST(Hyperplane::strong_cmp(h4,h2)<0){}
    CGAL_TEST(Hyperplane::weak_cmp(h4,h2)==0){}
    CGAL_TEST(Hyperplane::weak_cmp(h0,h1)>0){}
    CGAL_TEST(Hyperplane::weak_cmp(h1,h0)<0){}
    CGAL_TEST(h4!=h2 && CGAL::weak_equality(h4,h2)){}
  }

  {
    Point p(-1,-1,-1,1);
    Point q(1,-1,-1,1);
    Point r(-1,1,-1,1);
    Point s(1,1,1,1);

    std::vector< Point > A = make_vector(p,q,r,s);
    Sphere S1(3,A.begin(),A.end()), S2(3), S3(S1), S4 = S1;

    CGAL_TEST(S1 == S3 && S1 == S4 && S1 != S2){}
    CGAL_TEST(S1.dimension()==3){}
    CGAL_TEST(S1.point(0)==p){}
    CGAL_TEST(S2.is_degenerate()){}
    CGAL_TEST(S1.orientation()==CGAL::POSITIVE){}
    CGAL_TEST(S2.orientation()==CGAL::ZERO){}
    CGAL_TEST(S1.is_legal() && S2.is_legal()){}

    A[2]=A[3];
    Sphere S5 = Sphere(3,A.begin(),A.end());
    CGAL_TEST((!S5.is_legal())){}

    Point origin(0,0,0,1);
    CGAL_TEST((S1.center()==origin)){}
    CGAL_TEST(S1.has_on_bounded_side(origin)){}
    CGAL_TEST(S1.has_on_positive_side(origin)){}
    CGAL_TEST(!S1.has_on_unbounded_side(origin)){}
    CGAL_TEST(S1.has_on_unbounded_side(Point(2,0,0,1))){}
    CGAL_TEST(S1.has_on_negative_side(Point(2,0,0,1))){}
    CGAL_TEST(S1.has_on_boundary(Point(1,-1,1,1))){}
    CGAL_TEST(S1.oriented_side(origin)==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(S1.bounded_side(origin)==CGAL::ON_BOUNDED_SIDE){}

    Vector v(1,0,0,1);
    std::vector< Point > B = make_vector(p+v,q+v,r+v,s+v);
    CGAL_TEST( (S1+v) == Sphere(3,B.begin(),B.end()) ){}
    CGAL_TEST( CGAL::weak_equality(S1.opposite(),S1) ){}
    if (DOIO) CGAL_IO_TEST(S1,S3,CGAL::IO::ASCII); S3 = S1;
  }

  {
    Point p1(-5,1), p2(5,1);
    Vector v(0,5);
    Segment s0, s1(p1,p2), s2(p1,v),
              s3(Point(-1,-1),Point(5,5)),
              s4(Point(1,1,2),Point(5,5,3)),
              s5(Point(0,0,1),Point(4,4,1)),
              s6(s1);

    CGAL_TEST(s6==s1&&s2!=s1&&s5!=s2){}
    CGAL_TEST(s5.dimension()==2 && s2.dimension()==s3.dimension()){}
    CGAL_TEST(s1.source()==p1 && s1.target()==p2){}
    CGAL_TEST(s1.vertex(0)==p1 && s1.vertex(1)==p2){}
    CGAL_TEST(s1.point(0)==p1 && s1.point(1)==p2){}
    CGAL_TEST(s1[0]==p1 && s1[1]==p2){}
    CGAL_TEST(s1.opposite()==Segment(p2,p1)){}
    CGAL_TEST((s1.min)()==p1 && (s1.max)()==p2){}
    CGAL_TEST(s2.vector()==v){}
    CGAL_TEST(s2.direction()==v.direction()){}
    CGAL_TEST(s1.supporting_line()==Line(p1,p2)){}
    CGAL_TEST((s1 + v)==Segment(Point(-5,6),Point(5,6))){}

    s6 = Segment(Point(0,0),Point(10,10));
    Point p3(10,0);
    CGAL_TEST(s6.squared_length()==FT(200)){}
    CGAL_TEST(s6.has_on(Point(5,5))){}
    CGAL_TEST(Segment(Point(1,1),Point(1,1)).is_degenerate()){}
    CGAL_TEST(CGAL::weak_equality(s1,s1.opposite())){}
    CGAL_TEST(CGAL::parallel(s6,s6+v)){}
    CGAL_TEST(CGAL::common_endpoint(s1,s2,p3)&&p3==p1){}
    if (DOIO) CGAL_IO_TEST(s1,s6,CGAL::IO::ASCII); s6 = s1;
  }

  {
    Point p1(2), p2(5,0);
    Segment s(Point(1,-1),Point(1,5));
    Direction dir(1,0);
    Ray r0, r1(p1,p2), r2(p1,dir), r3(s), r4(r1);
    Point p3(3), p4(1,1,1,1);
    Ray r6(p3,p4);

    CGAL_TEST(r1==r4&&r1==r2){}
    CGAL_TEST(r1!=r3){}
    CGAL_TEST(r0.dimension()==0&&r1.dimension()==2&&r6.dimension()==3){}
    CGAL_TEST(r1.source()==p1){}
    CGAL_TEST(r1.point(0)==p1&&r1.point(1)==p2){}
    CGAL_TEST(r1.opposite()==Ray(p1,-dir)){}
    CGAL_TEST(r2.direction()==dir){}
    CGAL_TEST(r1.supporting_line()==r2.supporting_line()){}
    Ray r7(p1+Vector(1,-1),dir);
    CGAL_TEST(r1+Vector(1,-1) == r7){}
    CGAL_TEST(r1.has_on(Point(10,0))&&!r1.has_on(Point(-10,0))){}
    CGAL_TEST(CGAL::parallel(r1,(r2+Vector(0,1)).opposite())){}
    if (DOIO) CGAL_IO_TEST(r1,r4,CGAL::IO::ASCII); r4 = r1;
  }

  {
    Point origin(2), p1(5,0), p2(0,5), p3(3,3,2);
    Direction dir(7,0);
    Vector v(0,2);
    Segment s(Point(1,-2),Point(1,2));
    Ray r(origin, p3);
    Line l1(origin,p1), l2(origin,dir), l3(s), l4(r), l5(l1), l6;
    CGAL_TEST(l1==l5&&l1==l2){}
    CGAL_TEST(l1!=l3&&l1!=l6){}
    CGAL_TEST(l1.dimension()==2&&l6.dimension()==0){}
    CGAL_TEST(l1.point(0)==origin&&l1.point(1)==p1){}
    CGAL_TEST(l1.opposite()==Line(p1,-dir)){}
    CGAL_TEST(l2.direction()==dir){}
    CGAL_TEST(l1+v==Line(origin+v,dir)){}
    CGAL_TEST(l1.projection(Point(1,1))==Point(1,0)){}
    CGAL_TEST(l1.has_on(Point(3,0))&&!l1.has_on(p2)){}
    CGAL_TEST(CGAL::weak_equality(l1,l1.opposite())&&l1!=l1.opposite()){}
    CGAL_TEST(CGAL::parallel(l2,Line(p2,dir))){}
    if (DOIO) CGAL_IO_TEST(l1,l5,CGAL::IO::ASCII); l5 = l1;
  }




}
{ // Cartesian Kernel
  typedef CGAL::Cartesian_d<FT_> Kernel;
  typedef Kernel::FT             FT;
  typedef Kernel::RT             RT;
  typedef Kernel::LA             LA;
  typedef CGAL::Point_d<Kernel>      Point;
  typedef CGAL::Vector_d<Kernel>     Vector;
  typedef CGAL::Direction_d<Kernel>  Direction;
  typedef CGAL::Hyperplane_d<Kernel> Hyperplane;
  typedef CGAL::Sphere_d<Kernel>     Sphere;
  typedef CGAL::Segment_d<Kernel>    Segment;
  typedef CGAL::Ray_d<Kernel>        Ray;
  typedef CGAL::Line_d<Kernel>       Line;
  typedef LA::Vector                 IVector;
  bool DOIO=false;

  {
    /* some construction test */
    int IV1[] = {1,0,1};
    int IV2[] = {0,1};
    IVector iv1(IV1,IV1+3), iv2(IV2,IV2+2);
    Point p0(2,CGAL::ORIGIN), p1(2,IV1,IV1+3),
          p2(2,iv2.begin(),iv2.end(),1),
          p3(2), p4(p1), p5(2,iv2.begin(),iv2.end(),1);
    CGAL_TEST(p0 == p3 && p1 == p4 && p2 == p5){} // op==
    CGAL_TEST(p0 != p1 && p0 != p2){} // op!=
    CGAL_TEST(p0 == CGAL::ORIGIN){}

    if (DOIO) CGAL_IO_TEST(p1,p4,CGAL::IO::ASCII); p4 = p1;

    CGAL_TEST(p4.dimension()==p1.dimension()){} // dimension()
    CGAL_TEST(p4.dimension()!=Point(3).dimension()){} // dimension()

    /* some input and access test */
    CGAL_TEST(p1.hx()==RT(1) && p1.hy()==RT(0) && p1.hw()!=RT(0)){}
    CGAL_TEST(p1.hx()==p1.homogeneous(0) &&
              p1.hy()==p1.homogeneous(1) &&
              p1.hw()==p1.homogeneous(2)){} // hx,hy,hw
    CGAL_TEST(p1.x()==FT(1) && p1.y()==FT(0)){} // x(),y()
    CGAL_TEST(p1.x()==p1[0] && p1.y()==p1[1]){} // x(),y()
    CGAL_TEST(p1.x()==p1.cartesian(0) && p1.y()==p1.cartesian(1)){} // x(),y()
    Point::Homogeneous_const_iterator hcit; int i;
    for (i=0, hcit = p1.homogeneous_begin();
         hcit != p1.homogeneous_end(); ++i, ++hcit)
      CGAL_TEST(*hcit == iv1[i]){}

    Point::Homogeneous_const_iterator hit;
    Point::Cartesian_const_iterator cit;

    Kernel::Cartesian_const_iterator_d cit_begin =
      Kernel().construct_cartesian_const_iterator_d_object()(p1);
    Kernel::Cartesian_const_iterator_d cit_end =
      Kernel().construct_cartesian_const_iterator_d_object()(p1, 1);

    for (i=0,hit=p1.homogeneous_begin(),cit=p1.cartesian_begin();
         i<p1.dimension(); ++hit,++cit,++i) {
      CGAL_TEST(p1.homogeneous(i)==*hit){}
      CGAL_TEST(p1.cartesian(i)==*cit){}
      CGAL_TEST(p1.cartesian(i)==*cit_begin){}
      cit_begin++;
    }
    CGAL_TEST(cit==p1.cartesian_end()){}
    CGAL_TEST(cit_begin==cit_end){}
    CGAL_TEST(p1.homogeneous(p1.dimension())==*hit){}
    CGAL_TEST(++hit == p1.homogeneous_end()){}

    /* the geometric operations and predicates */

    CGAL_TEST(CGAL::project_along_d_axis(
              CGAL::lift_to_paraboloid(p1)) == p1){} // lift and project
    Vector vi(7,8,1);
    CGAL_TEST(CGAL::midpoint(p1+vi,p1-vi)==p1){}
    CGAL_TEST(CGAL::squared_distance(p1,p1+vi)==vi.squared_length()){}

    std::vector<Point> A = make_vector(p0,p1,p2);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::POSITIVE){}

    p3 = Point(1,1,2);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_BOUNDED_SIDE){}
    p3 = Point(1,1,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_ORIENTED_BOUNDARY){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_BOUNDARY){}
    p3 = Point(2,2,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_NEGATIVE_SIDE){}
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_UNBOUNDED_SIDE){}

    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p3) ==
          CGAL::ON_UNBOUNDED_SIDE){}

    Point B[2] = { p1, p2 };
    CGAL_TEST(CGAL::affine_rank(B,B+2)==1){}
    CGAL_TEST(CGAL::contained_in_affine_hull(B,B+2,Point(1,1,2))){}
    CGAL_TEST(CGAL::contained_in_affine_hull(B,B+2,Point(1,1,2)+5*(p1-p2))){}
    CGAL_TEST(CGAL::contained_in_simplex(B,B+2,Point(1,1,2))){}
    CGAL_TEST(!CGAL::contained_in_simplex(B,B+2,Point(5,5,1))){}
    CGAL_TEST(!CGAL::contained_in_simplex(B,B+2,Point(1,1,2)+5*(p1-p2))){}
    CGAL_TEST(CGAL::contained_in_simplex(A.begin(),A.end(),Point(1,1,3))){}
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),Point(5,5,1))){}
    CGAL_TEST(CGAL::affinely_independent(B,B+2)){}
    CGAL_TEST(CGAL::affinely_independent(A.begin(),A.end())){}
    A.push_back(Point(1,1));
    CGAL_TEST(!CGAL::affinely_independent(A.begin(),A.end())){}
    CGAL_TEST(CGAL::compare_lexicographically(p0,p1)==CGAL::SMALLER){}
    CGAL_TEST(CGAL::lexicographically_smaller(p0,p1)){}
    CGAL_TEST(CGAL::lexicographically_smaller_or_equal(p0,p1)){}
    CGAL_TEST(CGAL::compare_lexicographically(p1,p0)==CGAL::LARGER){}
    CGAL_TEST(CGAL::compare_lexicographically(p1,p1)==CGAL::EQUAL){}
  }

  {
    /* some construction test */
    Point p0(3); // the origin
    Vector e1 = Vector(3,Vector::Base_vector(),0);
    int IV1[] = {0,1,0,1};
    int IV2[] = {0,0,1};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    // the first unit vector

    Point p1(p0 + e1),
          p2(3,iv1.begin(),iv1.end()),
          p3(3,iv2.begin(),iv2.end(),1),
          p4(3);
    CGAL_TEST(CGAL::compare_lexicographically(p0,p4)==CGAL::EQUAL){}
    CGAL_TEST(CGAL::compare_lexicographically(p0,p1)==CGAL::SMALLER){}
    CGAL_TEST(CGAL::compare_lexicographically(p3,p0)==CGAL::LARGER){}

    CGAL_TEST(p1.cartesian(0)==FT(1)){}
    CGAL_TEST(p1[1]==FT(0)){}
    CGAL_TEST(p1.homogeneous(0)==Kernel::RT(1)){}
    CGAL_TEST( (p1 - CGAL::ORIGIN) == e1){}
    CGAL_TEST(CGAL::squared_distance(Point(1,1,1,1),p0)==FT(3)){}

    CGAL_TEST(p1-p0==e1){}
    CGAL_TEST(p1+e1==Point(2,0,0,1)){}
    CGAL_TEST(p1-e1==p0){}

    Point p5(3);
    CGAL_TEST((p5+=e1)==p1){}
    CGAL_TEST((p5-=e1)==p0){}

    /* orientation, sphere position, simplex position */
    std::vector<Point> A = make_vector(p0,p1,p2,p4);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::ZERO){}
    A[3] = p3;
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::POSITIVE){}
    std::swap(A[2],A[3]);
    CGAL_TEST(CGAL::orientation(A.begin(),A.end())==CGAL::NEGATIVE){}

    p4 = Point(1,1,1,2);
    std::swap(A[2],A[3]);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_POSITIVE_SIDE){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_ORIENTED_BOUNDARY){}
    p4 = Point(2,2,2,1);
    CGAL_TEST(CGAL::side_of_oriented_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_NEGATIVE_SIDE){}

    std::swap(A[2],A[3]);
    p4 = Point(1,1,1,2);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_BOUNDED_SIDE){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_BOUNDARY){}
    p4 = Point(2,2,2,1);
    CGAL_TEST(CGAL::side_of_bounded_sphere(A.begin(),A.end(),p4)
          == CGAL::ON_UNBOUNDED_SIDE){}
    std::swap(A[2],A[3]);

    p4 = Point(1,1,1,3);
    CGAL_TEST(CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}
    p4 = Point(1,1,1,2);
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}
    p4 = Point(1,1,1,1);
    CGAL_TEST(!CGAL::contained_in_simplex(A.begin(),A.end(),p4)){}

    /* affine hull, independence and rank */
    std::vector<Point> B = make_vector(p1,p2,p3);
    p4 = CGAL::midpoint(p1,p2);
    CGAL_TEST(CGAL::contained_in_affine_hull(B.begin(),B.end(),p4)){}

    p4 = p0;
    CGAL_TEST(!CGAL::contained_in_affine_hull(B.begin(),B.end(),p0)){}

    A[0] = Point(3,3,3,1);
    CGAL_TEST(CGAL::affinely_independent(A.begin(),A.end())){}
    CGAL_TEST(CGAL::affine_rank(A.begin(),A.end())==3){}
    A[3] = CGAL::midpoint(CGAL::midpoint(A[0],A[1]),A[2]);
    CGAL_TEST(!CGAL::affinely_independent(A.begin(),A.end())){}
    // degenerate_simplex
    CGAL_TEST(CGAL::affine_rank(A.begin(),A.end())==2){}
    // full_simplex
  }


  {
    /* construction and access */
    int IV1[] = {1,2,3,1};
    int IV2[] = {3,2,1};
    IVector i1(IV1,IV1+4), i2(IV2,IV2+3);
    Vector a0(3),
           a1(3,i1.begin(),i1.end()),
           a2(3,i2.begin(),i2.end(),1),
           a3(6,4,2,2), a4(1,0,1), a5(1,1,6), a6(a1);

    CGAL_TEST((a1==a6 && a3==a2)){}
    CGAL_TEST((a0!=a1 && a0!=a4)){}

    if (DOIO) CGAL_IO_TEST(a1,a6,CGAL::IO::ASCII); a6 = a1;

    CGAL_TEST(a2.dimension()==a3.dimension()){}
    CGAL_TEST(a3.dimension()!=a4.dimension()){}

    CGAL_TEST(a4.hx()==RT(1) && a4.hy()==RT(0) && a4.hw()!=RT(0)){}
    CGAL_TEST(a4.x()==FT(1) && a4.y()==FT(0)){}
    CGAL_TEST(Kernel::make_FT(a3.homogeneous(1),a3.homogeneous(3))==FT(2)){}
    CGAL_TEST(a3.cartesian(1)==FT(2)){}
    CGAL_TEST(a3[0]==FT(3)){}
    CGAL_TEST(a3.direction()==Direction(3,2,1)){}

    Vector::Homogeneous_const_iterator hit;
    Vector::Cartesian_const_iterator cit; int i;
    for (i=0,hit=a1.homogeneous_begin(),cit=a1.cartesian_begin();
         i<a1.dimension(); ++hit,++cit,++i) {
      CGAL_TEST(a1.homogeneous(i)==*hit){}
      CGAL_TEST(a1.cartesian(i)==*cit){}
    }
    CGAL_TEST(cit==a1.cartesian_end()){}
    CGAL_TEST(a1.homogeneous(a1.dimension())==*hit){}
    CGAL_TEST(++hit == a1.homogeneous_end()){}

    /* arithmetical operations and compare */
    CGAL_TEST(Vector::cmp(a0,a1)==CGAL::SMALLER){}
    CGAL_TEST(Vector::cmp(a1,Vector(1,2,3,1))==CGAL::EQUAL){}
    CGAL_TEST(Vector::cmp(a1,Vector(-1,0,0,1))==CGAL::LARGER){}

    CGAL_TEST(2*a1==Vector(2,4,6,1)){}
    CGAL_TEST(Kernel::make_FT(2,3)*a4==Vector(2,0,3)){}

    CGAL_TEST(a3/2==Vector(3,2,1,2)){}
    CGAL_TEST(a3/Kernel::make_FT(2,3)==Vector(9,6,3,2)){}

    Vector a7(a1);
    a7*=2;
    a7*=Kernel::make_FT(3,2);
    CGAL_TEST(3*a1==a7){}
    a7/=Kernel::make_FT(3,4);
    a7/=4;
    CGAL_TEST(a1==a7){}

    CGAL_TEST(a4*a4==FT(1)){}
    CGAL_TEST(a4.squared_length()==FT(1)){}
    CGAL_TEST(a1+a2==Vector(4,4,4,1)){}
    CGAL_TEST(a1-a2==Vector(-2,0,2,1)){}

    Vector eins(1,1,1,1);
    a7 = a1;
    CGAL_TEST((a7+=eins)==a1+eins){}
    CGAL_TEST((a7-=eins)==a1){}
    CGAL_TEST(a1-2*a1==-a1){}
    CGAL_TEST(a0.is_zero()){}

    /* some linear algebra */

    std::vector<Vector> P12 = make_vector(a1,a2);
    CGAL_TEST(CGAL::contained_in_linear_hull(P12.begin(),P12.end(),a1+a2)){}
    Vector z_off = a1 + Vector(3,Vector::Base_vector(),2);
    CGAL_TEST(!CGAL::contained_in_linear_hull(P12.begin(),P12.end(),z_off)){}

    std::vector<Vector> NLI = make_vector(a1,a2,a1+a2),
                        LI = make_vector(a1,a2,z_off);
    CGAL_TEST(CGAL::linearly_independent(LI.begin(),LI.end())){}
    CGAL_TEST(!CGAL::linearly_independent(NLI.begin(),NLI.end())){}

    CGAL_TEST(CGAL::linear_rank(LI.begin(),LI.end())==3){}
    CGAL_TEST(CGAL::linear_rank(NLI.begin(),NLI.end())==2){}
    std::vector<Vector> LB(3);
    std::vector<Vector>::iterator last =
      CGAL::linear_base(NLI.begin(),NLI.end(),LB.begin());
    CGAL_TEST(CGAL::linearly_independent(LB.begin(),last)){}
  }


  {
    /* construction and access */
    int IV1[] = {1,5,-3,2};
    int IV2[] = {1,1,0};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    Direction d0(3), d1(3,iv1.begin(),iv1.end()),
              d2(3,iv2.begin(),iv2.end()),
              d31(1,1,1), d32(1,1),
              d4 = Direction(3,Direction::Base_direction(),2),
              d5(d0);
    CGAL_TEST(d0==d5){}
    CGAL_TEST(d2==Direction(2,2,0)){}
    CGAL_TEST(d31!=d32 && d2!=d1)

    if (DOIO) CGAL_IO_TEST(d1,d5,CGAL::IO::ASCII); d5=d1;
    CGAL_TEST(d2.dimension()==d1.dimension()){}
    CGAL_TEST(d31.dimension()!=d32.dimension()){}
    CGAL_TEST(d1.dx()==d1.delta(0)){}
    CGAL_TEST(d1.dy()==d1.delta(1)){}
    CGAL_TEST(d1.dz()==d1.delta(2)){}
    CGAL_TEST(d1.dy()==d1[1]){}
    CGAL_TEST(d4.delta(0)==RT(0) && d4.delta(2)==RT(1)){}
    CGAL_TEST(d4==-Direction(0,0,-1)){}
    CGAL_TEST(d4==Direction(0,0,-1).opposite()){}
    CGAL_TEST(Direction::cmp(d0,d2)==CGAL::SMALLER){}
    Direction::Delta_const_iterator it; int i;
    for (i=0,it=d1.deltas_begin(); i<d1.dimension(); ++it,++i) {
      CGAL_TEST(d1.delta(i)==*it){}
    }
    CGAL_TEST(it==d1.deltas_end()){}
  }

  {
    /* construction and access */
    int IV1[] = {1,2,3,4};
    int IV2[] = {-4,-3,-2};
    IVector vi1(IV1,IV1+4);
    IVector vi2(IV2,IV2+3);
    // two ivec inits

    Point p1(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),0)),
          p2(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),1)),
          p3(CGAL::ORIGIN + Vector(3,Vector::Base_vector(),2));
    // one three point init

    Direction dir(1,1,1);
    Point o(-1,-1,-1,1), org(3);
    // one dir init
    Hyperplane h0(3,vi1.begin(),vi1.end());
    Hyperplane h1(3,vi2.begin(),vi2.end(),1);
    std::vector<Point> V = make_vector(p1,p2,p3);
    Hyperplane h2(V.begin(),V.end(),o,CGAL::ON_NEGATIVE_SIDE);
    V = make_vector(org,org,org);
    Hyperplane h3(V.begin(),V.end(),o,CGAL::ON_POSITIVE_SIDE);
    Hyperplane h4(p1,-dir);
    Hyperplane h5(1,2,3,4);
    Hyperplane h6(1,0,1);
    Hyperplane h7(h0);
    Hyperplane h8(3);

    CGAL_TEST((h5==h0 && h7==h0)){}
    CGAL_TEST((h1!=h0 && h7!=h6)){}
    CGAL_TEST((h2.dimension()==h5.dimension() &&
               h2.dimension()!=h6.dimension())){}
    CGAL_TEST((h5.coefficient(1)==RT(2) && h5[2]==RT(3))){}
    CGAL_TEST((vi1 == h0.coefficient_vector())){}
    CGAL_TEST(h2.orthogonal_vector()*(p2-p1) == FT(0) ){}
    CGAL_TEST(h2.orthogonal_direction()==dir){}
    CGAL_TEST(h2.value_at(o)<FT(0)){}
    CGAL_TEST(h2.value_at(p3)==FT(0)){}
    CGAL_TEST(h2.value_at(p3+dir.vector())>FT(0)){}
    CGAL_TEST(h2.oriented_side(o)==CGAL::ON_NEGATIVE_SIDE){}
    CGAL_TEST(h2.oriented_side(p3)==CGAL::ON_ORIENTED_BOUNDARY){}
    CGAL_TEST(h2.oriented_side(p3+dir.vector())==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(h3.oriented_side(o)==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(h4.has_on(p1)){}
    CGAL_TEST(h3.has_on(org)){}
    CGAL_TEST(!h5.has_on(p1)){}
    CGAL_TEST(h3.has_on_positive_side(o)){}

    if (DOIO) CGAL_IO_TEST(h0,h8,CGAL::IO::ASCII); h8=h0;

    Hyperplane::Coefficient_const_iterator it; int i;
    for (i=0,it=h5.coefficients_begin();
         i<=h5.dimension(); ++it,++i) {
      CGAL_TEST(h5.coefficient(i)==*it){}
    }
    CGAL_TEST(it==h5.coefficients_end()){}

    // all kind of |HyperplaneHd<integer>| inits

    /* compares and other operations */
    CGAL_TEST(Hyperplane::strong_cmp(h0,h1)>0){}
    CGAL_TEST(Hyperplane::strong_cmp(h1,h0)<0){}
    CGAL_TEST(Hyperplane::strong_cmp(h0,h5)==0){}
    CGAL_TEST(Hyperplane::strong_cmp(h4,h2)<0){}
    CGAL_TEST(Hyperplane::weak_cmp(h4,h2)==0){}
    CGAL_TEST(Hyperplane::weak_cmp(h0,h1)>0){}
    CGAL_TEST(Hyperplane::weak_cmp(h1,h0)<0){}
    CGAL_TEST(h4!=h2 && CGAL::weak_equality(h4,h2)){}
  }

  {
    Point p(-1,-1,-1,1);
    Point q(1,-1,-1,1);
    Point r(-1,1,-1,1);
    Point s(1,1,1,1);

    std::vector< Point > A = make_vector(p,q,r,s);
    Sphere S1(3,A.begin(),A.end()), S2(3), S3(S1), S4 = S1;

    CGAL_TEST(S1 == S3 && S1 == S4 && S1 != S2){}
    CGAL_TEST(S1.dimension()==3){}
    CGAL_TEST(S1.point(0)==p){}
    CGAL_TEST(S2.is_degenerate()){}
    CGAL_TEST(S1.orientation()==CGAL::POSITIVE){}
    CGAL_TEST(S2.orientation()==CGAL::ZERO){}
    CGAL_TEST(S1.is_legal() && S2.is_legal()){}

    A[2]=A[3];
    Sphere S5 = Sphere(3,A.begin(),A.end());
    CGAL_TEST((!S5.is_legal())){}

    Point origin(0,0,0,1);
    CGAL_TEST((S1.center()==origin)){}
    CGAL_TEST(S1.has_on_bounded_side(origin)){}
    CGAL_TEST(S1.has_on_positive_side(origin)){}
    CGAL_TEST(!S1.has_on_unbounded_side(origin)){}
    CGAL_TEST(S1.has_on_unbounded_side(Point(2,0,0,1))){}
    CGAL_TEST(S1.has_on_negative_side(Point(2,0,0,1))){}
    CGAL_TEST(S1.has_on_boundary(Point(1,-1,1,1))){}
    CGAL_TEST(S1.oriented_side(origin)==CGAL::ON_POSITIVE_SIDE){}
    CGAL_TEST(S1.bounded_side(origin)==CGAL::ON_BOUNDED_SIDE){}

    Vector v(1,0,0,1);
    std::vector< Point > B = make_vector(p+v,q+v,r+v,s+v);
    CGAL_TEST( (S1+v) == Sphere(3,B.begin(),B.end()) ){}
    CGAL_TEST( CGAL::weak_equality(S1.opposite(),S1) ){}
    if (DOIO) CGAL_IO_TEST(S1,S3,CGAL::IO::ASCII); S3 = S1;
  }

  {
    Point p1(-5,1), p2(5,1);
    Vector v(0,5);
    Segment s0, s1(p1,p2), s2(p1,v),
              s3(Point(-1,-1),Point(5,5)),
              s4(Point(1,1,2),Point(5,5,3)),
              s5(Point(0,0,1),Point(4,4,1)),
              s6(s1);

    CGAL_TEST(s6==s1&&s2!=s1&&s5!=s2){}
    CGAL_TEST(s5.dimension()==2 && s2.dimension()==s3.dimension()){}
    CGAL_TEST(s1.source()==p1 && s1.target()==p2){}
    CGAL_TEST(s1.vertex(0)==p1 && s1.vertex(1)==p2){}
    CGAL_TEST(s1.point(0)==p1 && s1.point(1)==p2){}
    CGAL_TEST(s1[0]==p1 && s1[1]==p2){}
    CGAL_TEST(s1.opposite()==Segment(p2,p1)){}
    CGAL_TEST((s1.min)()==p1 && (s1.max)()==p2){}
    CGAL_TEST(s2.vector()==v){}
    CGAL_TEST(s2.direction()==v.direction()){}
    CGAL_TEST(s1.supporting_line()==Line(p1,p2)){}
    CGAL_TEST((s1 + v)==Segment(Point(-5,6),Point(5,6))){}

    s6 = Segment(Point(0,0),Point(10,10));
    Point p3(10,0);
    CGAL_TEST(s6.squared_length()==FT(200)){}
    CGAL_TEST(s6.has_on(Point(5,5))){}
    CGAL_TEST(Segment(Point(1,1),Point(1,1)).is_degenerate()){}
    CGAL_TEST(CGAL::weak_equality(s1,s1.opposite())){}
    CGAL_TEST(CGAL::parallel(s6,s6+v)){}
    CGAL_TEST(CGAL::common_endpoint(s1,s2,p3)&&p3==p1){}
    if (DOIO) CGAL_IO_TEST(s1,s6,CGAL::IO::ASCII); s6 = s1;
  }

  {
    Point p1(2), p2(5,0);
    Segment s(Point(1,-1),Point(1,5));
    Direction dir(1,0);
    Ray r0, r1(p1,p2), r2(p1,dir), r3(s), r4(r1);
    Point p3(3), p4(1,1,1,1);
    Ray r6(p3,p4);

    CGAL_TEST(r1==r4&&r1==r2){}
    CGAL_TEST(r1!=r3){}
    CGAL_TEST(r0.dimension()==0&&r1.dimension()==2&&r6.dimension()==3){}
    CGAL_TEST(r1.source()==p1){}
    CGAL_TEST(r1.point(0)==p1&&r1.point(1)==p2){}
    CGAL_TEST(r1.opposite()==Ray(p1,-dir)){}
    CGAL_TEST(r2.direction()==dir){}
    CGAL_TEST(r1.supporting_line()==r2.supporting_line()){}
    Ray r7(p1+Vector(1,-1),dir);
    CGAL_TEST(r1+Vector(1,-1) == r7){}
    CGAL_TEST(r1.has_on(Point(10,0))&&!r1.has_on(Point(-10,0))){}
    CGAL_TEST(CGAL::parallel(r1,(r2+Vector(0,1)).opposite())){}
    if (DOIO) CGAL_IO_TEST(r1,r4,CGAL::IO::ASCII); r4 = r1;
  }

  {
    Point origin(2), p1(5,0), p2(0,5), p3(3,3,2);
    Direction dir(7,0);
    Vector v(0,2);
    Segment s(Point(1,-2),Point(1,2));
    Ray r(origin, p3);
    Line l1(origin,p1), l2(origin,dir), l3(s), l4(r), l5(l1), l6;
    CGAL_TEST(l1==l5&&l1==l2){}
    CGAL_TEST(l1!=l3&&l1!=l6){}
    CGAL_TEST(l1.dimension()==2&&l6.dimension()==0){}
    CGAL_TEST(l1.point(0)==origin&&l1.point(1)==p1){}
    CGAL_TEST(l1.opposite()==Line(p1,-dir)){}
    CGAL_TEST(l2.direction()==dir){}
    CGAL_TEST(l1+v==Line(origin+v,dir)){}
    CGAL_TEST(l1.projection(Point(1,1))==Point(1,0)){}
    CGAL_TEST(l1.has_on(Point(3,0))&&!l1.has_on(p2)){}
    CGAL_TEST(CGAL::weak_equality(l1,l1.opposite())&&l1!=l1.opposite()){}
    CGAL_TEST(CGAL::parallel(l2,Line(p2,dir))){}
    if (DOIO) CGAL_IO_TEST(l1,l5,CGAL::IO::ASCII); l5 = l1;
  }




}
  CGAL_TEST_END;
}

