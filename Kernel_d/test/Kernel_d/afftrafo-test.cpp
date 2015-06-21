#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/simple_transformations_d.h>
#include <CGAL/test_macros.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>


typedef CGAL::Exact_integer RT_;
typedef CGAL::Exact_rational FT_;


int main()
{ 
  CGAL_TEST_START;
{ 
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
  typedef LA::Vector                  IVector;
  typedef LA::Matrix                  IMatrix;
  typedef CGAL::Aff_transformation_d<Kernel> Afftrafo;

  { 
    Vector e1 = Vector(3,Vector::Base_vector(),0);
    int IV1[]={1,0,1,1};
    int IV2[]={0,0,1};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    Point p0(3), p1(p0 + e1), p2(3,iv1.begin(),iv1.end()), 
          p3(3,iv2.begin(),iv2.end(),1), p4(3);
    std::vector<Point> V = 
      make_vector(Point(-1,-1,-1,1),Point(1,-1,-1,1),
                  Point(-1,1,-1,1),Point(-1,-1,1,1));
    Sphere S(3,V.begin(),V.end());
    Line l1(p0,p1), l2(p0,Point(1,1,0,1)); 
    Ray r(p0,p1); Segment s(p0,p1);
    CGAL_TEST(CGAL::reflect(p1,p0)==Point(-1,0,0,1)){}
    CGAL_TEST(CGAL::reflect(p2,l1)==Point(1,0,-1,1)){}
    CGAL_TEST(CGAL::reflect(S,p0)==S.opposite()){}
    CGAL_TEST(CGAL::reflect(S,l1)==S){}
    CGAL_TEST(CGAL::reflect(l1,p0)==l1.opposite()){}
    CGAL_TEST(CGAL::reflect(l1,l2)==Line(p0,Point(0,1,0,1))){}
    CGAL_TEST(CGAL::reflect(r,p0)==r.opposite()){}
    CGAL_TEST(CGAL::reflect(r,l2)==Ray(p0,Point(0,1,0,1))){}
    CGAL_TEST(CGAL::reflect(s,p0)==Segment(p0,Point(-1,0,0,1))){}
    CGAL_TEST(CGAL::reflect(s,l2)==Segment(p0,Point(0,1,0,1))){}
  }

  { 
    int IV[]={2,3,1};
    IVector iv(IV,IV+3); 
    IMatrix M1(3),M2(3); 
    M1(0,0) = -1; M1(1,1) = -1; M1(2,2) = 1; // rot pi
    M2(0,0) = 0; M2(0,1) = -1; M2(0,2) = 1;
    M2(1,0) = 1; M2(1,1) = 0; M2(1,2) = 1;  M2(2,2) = 1;
    Vector rv(5,0);
    Direction rd(1,1);
   
    Afftrafo at0(2); // whole zero in 2d
    Afftrafo at1(M1); // rot(pi) in 2d
    Afftrafo at2(CGAL::SCALING,iv.begin(),iv.end()); // scale by (2,3) in 2d
    Afftrafo at3(CGAL::TRANSLATION,rv); // translate by (5,0) in 2d
    Afftrafo at4(at1); // rot(pi) in 2d
    Afftrafo at5(M2); // rot(-pi/2) + trans(1,1) in 2d
    Afftrafo at6(2,CGAL::IDENTITY); // id trafo in 2d
    Afftrafo at9(2,CGAL::SCALING,3,1);  // scale(3/1) in 2d
    Afftrafo at10(2,CGAL::ROTATION,RT(1),RT(0),RT(1)); // rot(pi/2) in 2d
    Afftrafo at11(2,CGAL::ROTATION,rd,1,1000); // d2_rot approx

    CGAL_TEST(at1==at4){}
    CGAL_TEST(at1!=at5){}
    CGAL_TEST(at1.dimension()==2){}
    CGAL_TEST(at1.matrix() == M1){}

    CGAL_IO_TEST(at4,at0,CGAL::IO::ASCII);

    Afftrafo at12 = at2 * at2.inverse();
    Point p(1,1); 
    CGAL_TEST(p.transform(at12)==p){} 
    CGAL_TEST(p.transform(at1)==Point(-1,-1)){}
    CGAL_TEST(p.transform(at2)==Point(2,3)){}
    CGAL_TEST(p.transform(at3)==Point(6,1)){}
    CGAL_TEST(p.transform(at4)==Point(-1,-1)){}
    CGAL_TEST(p.transform(at5)==Point(0,2)){}
    CGAL_TEST(p.transform(at6)==p){}
    CGAL_TEST(p.transform(at9)==Point(3,3)){}
    CGAL_TEST(p.transform(at10)==Point(-1,1)){}

    Vector v_approx = 
      Point(1,0).transform(at11) - CGAL::ORIGIN;
    /* we have v_approx = (x,y) building an angle alpha with
       the x-axis then sin alpha = y/length(v_approx) and we
       want to know that sin alpha <= num/den of at11. 
       we started with a direction rd(1,1) with angle pi/4
       and thus sin pi/4 = 1/ sqrt(2)
       Thus we check that the squared difference of the sine
       values is bounded by sqr(num/den):
       
       sqr(y/sqrt((sqr(x) + sqr(y))) - 1/sqrt(2)) < sqr(num/den)

           sqr(a - b) < sqr(c)
       <=> sqr(a) + sqr(b) - sqr(c) < 2ab
       <=> sqr( sqr(a) + sqr(b) - sqr(c) ) < 4 sqr(a) sqr(b)
           where sqr(b) = 1/2
       <=> sqr( sqr(a) + sqr(b) - sqr(c) ) < 2 sqr(a)
       
       in at11 we had num = 1 and den = 1000
    */
    FT sqr_sin_alpha = CGAL_NTS square(v_approx.y())/
                       v_approx.squared_length();
    // sqr(a)
    FT sqr_sin_piforth(1,2);
    // sqr(b)
    FT sqr_bound = Kernel::make_FT(RT(1),CGAL_NTS square(RT(1000)));  
    // sqr(c)
    FT sqr_sum = sqr_sin_alpha + sqr_sin_piforth - sqr_bound;
    CGAL_TEST( (sqr_sum*sqr_sum) < (FT(2)*sqr_sin_alpha) ){}

    // and analogously:
    FT sqr_cos_alpha = CGAL_NTS square(v_approx.x())/
                       v_approx.squared_length();
    FT sqr_cos_piforth(1,2);
    sqr_sum = sqr_cos_alpha + sqr_cos_piforth - sqr_bound;
    CGAL_TEST((sqr_sum*sqr_sum) < FT(2)*sqr_cos_alpha){}

    Vector v(2,2);
    CGAL_TEST(v.transform(at9)==3*v){}
    CGAL_TEST(v.transform(at3)==v){}
    CGAL_TEST(v.transform(at1)==-v){}

    Direction dir(5,1); 
    CGAL_TEST(dir.transform(at9)==dir){}
    CGAL_TEST(dir.transform(at3)==dir){}
    CGAL_TEST(dir.transform(at1)==-dir){}

    // 3 dim hyperplane:
    Afftrafo at20 = Afftrafo(CGAL::TRANSLATION,Vector(5,0,0,1));
    Afftrafo at21 = Afftrafo(3,CGAL::ROTATION,RT(1),RT(0),RT(1),0,2);
    // pi/2-rotation in x-z-plane in 3-space

    Afftrafo at22 = at21*at20; // composition of both

    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), o(3);
    CGAL_TEST(p1.transform(at22)==p1.transform(at20).transform(at21)){}

    std::vector<Point> U = make_vector(p1,p2,p3);
    Hyperplane hyp(U.begin(),U.end(),o,CGAL::ON_NEGATIVE_SIDE);

    Point pp1 = p1.transform(at22),
          pp2 = p2.transform(at22),
          pp3 = p3.transform(at22),
          o_trafo = o.transform(at22);

    std::vector<Point> V = make_vector(pp1,pp2,pp3);
    Hyperplane hyp_trafo(V.begin(),V.end(),o_trafo,
                         CGAL::ON_NEGATIVE_SIDE);

    CGAL_TEST(hyp.transform(at22)==hyp_trafo){}

    Segment s(p1,p2);
    CGAL_TEST(s.transform(at22)==Segment(pp1,pp2)){}
    Ray r(p1,p2);
    CGAL_TEST(r.transform(at22)==Ray(pp1,pp2)){}
    Line l(p1,p2);
    CGAL_TEST(l.transform(at22)==Line(pp1,pp2)){}
    std::vector<Point> W = make_vector(p1,p2,p3,o);
    Sphere S(3,W.begin(),W.end());
    W = make_vector(pp1,pp2,pp3,o_trafo);
    CGAL_TEST(S.transform(at22)==Sphere(3,W.begin(),W.end())){}
  }


}
{
  typedef CGAL::Cartesian_d<FT_> Kernel;
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
  typedef LA::Vector                  IVector;
  typedef LA::Matrix                  IMatrix;
  typedef CGAL::Aff_transformation_d<Kernel> Afftrafo;


  { 
    Vector e1 = Vector(3,Vector::Base_vector(),0);
    int IV1[]={1,0,1,1};
    int IV2[]={0,0,1};
    IVector iv1(IV1,IV1+4), iv2(IV2,IV2+3);
    Point p0(3), p1(p0 + e1), p2(3,iv1.begin(),iv1.end()), 
          p3(3,iv2.begin(),iv2.end(),1), p4(3);
    std::vector<Point> V = 
      make_vector(Point(-1,-1,-1,1),Point(1,-1,-1,1),
                  Point(-1,1,-1,1),Point(-1,-1,1,1));
    Sphere S(3,V.begin(),V.end());
    Line l1(p0,p1), l2(p0,Point(1,1,0,1)); 
    Ray r(p0,p1); Segment s(p0,p1);
    CGAL_TEST(CGAL::reflect(p1,p0)==Point(-1,0,0,1)){}
    CGAL_TEST(CGAL::reflect(p2,l1)==Point(1,0,-1,1)){}
    CGAL_TEST(CGAL::reflect(S,p0)==S.opposite()){}
    CGAL_TEST(CGAL::reflect(S,l1)==S){}
    CGAL_TEST(CGAL::reflect(l1,p0)==l1.opposite()){}
    CGAL_TEST(CGAL::reflect(l1,l2)==Line(p0,Point(0,1,0,1))){}
    CGAL_TEST(CGAL::reflect(r,p0)==r.opposite()){}
    CGAL_TEST(CGAL::reflect(r,l2)==Ray(p0,Point(0,1,0,1))){}
    CGAL_TEST(CGAL::reflect(s,p0)==Segment(p0,Point(-1,0,0,1))){}
    CGAL_TEST(CGAL::reflect(s,l2)==Segment(p0,Point(0,1,0,1))){}
  }

  { 
    int IV[]={2,3,1};
    IVector iv(IV,IV+3); 
    IMatrix M1(3),M2(3); 
    M1(0,0) = -1; M1(1,1) = -1; M1(2,2) = 1; // rot pi
    M2(0,0) = 0; M2(0,1) = -1; M2(0,2) = 1;
    M2(1,0) = 1; M2(1,1) = 0; M2(1,2) = 1;  M2(2,2) = 1;
    Vector rv(5,0);
    Direction rd(1,1);
   
    Afftrafo at0(2); // whole zero in 2d
    Afftrafo at1(M1); // rot(pi) in 2d
    Afftrafo at2(CGAL::SCALING,iv.begin(),iv.end()); // scale by (2,3) in 2d
    Afftrafo at3(CGAL::TRANSLATION,rv); // translate by (5,0) in 2d
    Afftrafo at4(at1); // rot(pi) in 2d
    Afftrafo at5(M2); // rot(-pi/2) + trans(1,1) in 2d
    Afftrafo at6(2,CGAL::IDENTITY); // id trafo in 2d
    Afftrafo at9(2,CGAL::SCALING,3,1);  // scale(3/1) in 2d
    Afftrafo at10(2,CGAL::ROTATION,RT(1),RT(0),RT(1)); // rot(pi/2) in 2d
    Afftrafo at11(2,CGAL::ROTATION,rd,1,1000); // d2_rot approx

    CGAL_TEST(at1==at4){}
    CGAL_TEST(at1!=at5){}
    CGAL_TEST(at1.dimension()==2){}
    CGAL_TEST(at1.matrix() == M1){}

    CGAL_IO_TEST(at4,at0,CGAL::IO::ASCII);

    Afftrafo at12 = at2 * at2.inverse();
    Point p(1,1); 
    CGAL_TEST(p.transform(at12)==p){} 
    CGAL_TEST(p.transform(at1)==Point(-1,-1)){}
    CGAL_TEST(p.transform(at2)==Point(2,3)){}
    CGAL_TEST(p.transform(at3)==Point(6,1)){}
    CGAL_TEST(p.transform(at4)==Point(-1,-1)){}
    CGAL_TEST(p.transform(at5)==Point(0,2)){}
    CGAL_TEST(p.transform(at6)==p){}
    CGAL_TEST(p.transform(at9)==Point(3,3)){}
    CGAL_TEST(p.transform(at10)==Point(-1,1)){}

    Vector v_approx = 
      Point(1,0).transform(at11) - CGAL::ORIGIN;
    /* we have v_approx = (x,y) building an angle alpha with
       the x-axis then sin alpha = y/length(v_approx) and we
       want to know that sin alpha <= num/den of at11. 
       we started with a direction rd(1,1) with angle pi/4
       and thus sin pi/4 = 1/ sqrt(2)
       Thus we check that the squared difference of the sine
       values is bounded by sqr(num/den):
       
       sqr(y/sqrt((sqr(x) + sqr(y))) - 1/sqrt(2)) < sqr(num/den)

           sqr(a - b) < sqr(c)
       <=> sqr(a) + sqr(b) - sqr(c) < 2ab
       <=> sqr( sqr(a) + sqr(b) - sqr(c) ) < 4 sqr(a) sqr(b)
           where sqr(b) = 1/2
       <=> sqr( sqr(a) + sqr(b) - sqr(c) ) < 2 sqr(a)
       
       in at11 we had num = 1 and den = 1000
    */
    FT sqr_sin_alpha = CGAL_NTS square(v_approx.y())/
                       v_approx.squared_length();
    // sqr(a)
    FT sqr_sin_piforth(1,2);
    // sqr(b)
    FT sqr_bound = Kernel::make_FT(RT(1),CGAL_NTS square(RT(1000)));  
    // sqr(c)
    FT sqr_sum = sqr_sin_alpha + sqr_sin_piforth - sqr_bound;
    CGAL_TEST( (sqr_sum*sqr_sum) < (FT(2)*sqr_sin_alpha) ){}

    // and analogously:
    FT sqr_cos_alpha = CGAL_NTS square(v_approx.x())/
                       v_approx.squared_length();
    FT sqr_cos_piforth(1,2);
    sqr_sum = sqr_cos_alpha + sqr_cos_piforth - sqr_bound;
    CGAL_TEST((sqr_sum*sqr_sum) < FT(2)*sqr_cos_alpha){}

    Vector v(2,2);
    CGAL_TEST(v.transform(at9)==3*v){}
    CGAL_TEST(v.transform(at3)==v){}
    CGAL_TEST(v.transform(at1)==-v){}

    Direction dir(5,1); 
    CGAL_TEST(dir.transform(at9)==dir){}
    CGAL_TEST(dir.transform(at3)==dir){}
    CGAL_TEST(dir.transform(at1)==-dir){}

    // 3 dim hyperplane:
    Afftrafo at20 = Afftrafo(CGAL::TRANSLATION,Vector(5,0,0,1));
    Afftrafo at21 = Afftrafo(3,CGAL::ROTATION,RT(1),RT(0),RT(1),0,2);
    // pi/2-rotation in x-z-plane in 3-space

    Afftrafo at22 = at21*at20; // composition of both

    Point p1(1,0,0,1), p2(0,1,0,1), p3(0,0,1,1), o(3);
    CGAL_TEST(p1.transform(at22)==p1.transform(at20).transform(at21)){}

    std::vector<Point> U = make_vector(p1,p2,p3);
    Hyperplane hyp(U.begin(),U.end(),o,CGAL::ON_NEGATIVE_SIDE);

    Point pp1 = p1.transform(at22),
          pp2 = p2.transform(at22),
          pp3 = p3.transform(at22),
          o_trafo = o.transform(at22);

    std::vector<Point> V = make_vector(pp1,pp2,pp3);
    Hyperplane hyp_trafo(V.begin(),V.end(),o_trafo,
                         CGAL::ON_NEGATIVE_SIDE);

    CGAL_TEST(hyp.transform(at22)==hyp_trafo){}

    Segment s(p1,p2);
    CGAL_TEST(s.transform(at22)==Segment(pp1,pp2)){}
    Ray r(p1,p2);
    CGAL_TEST(r.transform(at22)==Ray(pp1,pp2)){}
    Line l(p1,p2);
    CGAL_TEST(l.transform(at22)==Line(pp1,pp2)){}
    std::vector<Point> W = make_vector(p1,p2,p3,o);
    Sphere S(3,W.begin(),W.end());
    W = make_vector(pp1,pp2,pp3,o_trafo);
    CGAL_TEST(S.transform(at22)==Sphere(3,W.begin(),W.end())){}
  }


}
  CGAL_TEST_END;
}

