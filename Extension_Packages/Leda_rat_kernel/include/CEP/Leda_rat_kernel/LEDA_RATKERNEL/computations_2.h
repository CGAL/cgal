#ifndef CEP_LEDA_RAT_COMPUTATIONS_2_H
#define CEP_LEDA_RAT_COMPUTATIONS_2_H

// LEDA rational kernel 2d computation classes ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

#include <LEDA/rat_point.h>
#include <LEDA/rat_vector.h>
#include <LEDA/rat_line.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_triangle.h>
#include <LEDA/rat_rectangle.h>
#include <LEDA/rat_ray.h>
#include <LEDA/rat_circle.h>


#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE

template<class HELP_KERNEL>
class CGAL_compute_leda_rat_squared_distance_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rational           result_type;
  
  typedef typename HELP_KERNEL::Compute_squared_distance_2  Compute_squared_distance_2;

  template<class T1, class T2>
  leda_rational operator()(const T1& obj1, const T2& obj2) const
  {
     leda_to_cgal_2 conv;
     Compute_squared_distance_2 compute;
     
     CGAL::Quotient<leda_integer> result = compute(conv(obj1), conv(obj2));
     return leda_rational(result.numerator(), result.denominator());  
  }

};

#ifndef CGAL_NO_DEPRECATED_CODE
template<class HELP_KERNEL>
class CGAL_compute_leda_rat_y_at_x_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rational           result_type;
  
  typedef typename HELP_KERNEL::Compute_y_at_x_2  Compute_y_at_x_2;

  leda_rational operator()(const leda_rat_line& l, const leda_rational& x) const
  {
     leda_to_cgal_2 conv;
     Compute_y_at_x_2 compute;
     
     CGAL::Quotient<leda_integer> result = compute(conv(l), CGAL::Quotient<leda_integer>(x.numerator(), x.denominator()) );
     return leda_rational(result.numerator(), result.denominator());  
  }
};
#endif

/*

class Compute_leda_rat_squared_distance_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rational           result_type;

  // point + other object ...

  leda_rational operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
    return p1.sqr_dist(p2);
  }
  
  leda_rational operator()(const leda_rat_point& p, const leda_rat_line& l) const
  {
    return l.sqr_dist(p);
  }
  
  leda_rational operator()(const leda_rat_line& l, const leda_rat_point& p) const
  {
    return l.sqr_dist(p);
  } 
  
  // optimize this later (?)
  leda_rational operator()(const leda_rat_point& p, const leda_rat_ray& r) const
  {   
    // two cases: p lies in on of the halfplanes defined by the rotation of the supporting line
    // of r ...
    
    leda_rat_point src  = r.point1();
    
    leda_rat_ray r2 = r.rotate90(src);
    
    int ori = LEDA_NAMESPACE_NAME::orientation(r2,p);
    
    if (ori == 1) return src.sqr_dist(p);
    
    // distance to line is the correct distance ...
    leda_rat_line l(src,r.point2());
    return l.sqr_dist(p);
  }
  
  leda_rational operator()(const leda_rat_ray& r, const leda_rat_point& p) const
  { 
    return this->operator()(p,r);
  }  
  
  leda_rational operator()(const leda_rat_point& p, const leda_rat_segment& s) const
  {
    return s.sqr_dist(p);
  }
  
  leda_rational operator()(const leda_rat_segment& s, const leda_rat_point& p) const
  {
    return s.sqr_dist(p);  
  }  
  
  leda_rational operator()(const leda_rat_point& p, const leda_rat_triangle& t) const
  {
    // inside or outside ???
    // if the point is inside, squared dist == 0
    // else take the shortest distance to a segment
    LEDA_NAMESPACE_NAME::region_kind rk = t.region_of(p);
    
    if (! (rk == LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) return 0;
    
    // point is outside; find lowest dist to a segment of t ...
    leda_rat_point p1 = t.point1();
    leda_rat_point p2 = t.point2();
    leda_rat_point p3 = t.point3();        
    
    if (t.is_degenerate()){
      if (p1 == p2) {
        return leda_rat_segment(p1,p3).sqr_dist(p);
      }
      else {
        if (p1 == p3) {
	  return leda_rat_segment(p1,p2).sqr_dist(p);
	}
	else { // p2==p3
	  return leda_rat_segment(p1,p2).sqr_dist(p);
	}
      }
    }
    // non - degenerate triangle ...
    int ori = t.orientation();
    
    if (ori == +1) { // ccw ...
      int o1 = LEDA_NAMESPACE_NAME::orientation(p1,p2,p);
      
      if (o1 != 1) { return leda_rat_segment(p1,p2).sqr_dist(p); }
      
      int o2 = LEDA_NAMESPACE_NAME::orientation(p2,p3,p);
      
      if (o2 != 1) { return leda_rat_segment(p2,p3).sqr_dist(p); }
      
      return leda_rat_segment(p3,p1).sqr_dist(p);
    }
    
    // clockwise ...
    int o1 = LEDA_NAMESPACE_NAME::orientation(p2,p1,p);
      
    if (o1 != 1) { return leda_rat_segment(p1,p2).sqr_dist(p); }
      
    int o2 = LEDA_NAMESPACE_NAME::orientation(p3,p2,p);
      
    if (o2 != 1) { return leda_rat_segment(p2,p3).sqr_dist(p); }
      
    return leda_rat_segment(p3,p1).sqr_dist(p);    
  }  
  
  leda_rational operator()(const leda_rat_triangle& t, const leda_rat_point& p) const
  {
    return this->operator()(p,t);
  }
  
  // line ... 
  leda_rational operator()(const leda_rat_line& l1, const leda_rat_line& l2) const
  {
    if (LEDA_NAMESPACE_NAME::cmp_signed_dist(l1.point1(),l1.point2(),l2.point1(),l2.point2()) == 0) {
      // parallel lines ...
      return l1.sqr_dist(l2.point1());
    }
    return 0; // the lines intersect ...
  }
  
  leda_rational operator()(const leda_rat_line& l, const leda_rat_ray& r) const
  {  
    int ori1 = l.side_of(r.point1());
    int ori2 = l.side_of(r.point2());  
    
    if (ori1!=ori2 || ori1*ori2==0) return 0; 
    
    //distance comparison
    int cm = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(),l.point2(),r.point1(),r.point2());
    
    if (cm == 1) return 0; // source has larger dist from line ...
    
    return l.sqr_dist(r.point1());
  }
  
  leda_rational operator()(const leda_rat_ray& r, const leda_rat_line& l) const
  { 
    return this->operator()(l,r);      
  }    
  
  leda_rational operator()(const leda_rat_line& l, const leda_rat_segment& s) const
  {  
    bool inter = l.intersection(s);
    
    if (inter) return 0;
    int cm = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(),l.point2(),s.start(),s.end());
    
    if (cm >= 0){ // s.end() has smaller or equal distance ...
      return l.sqr_dist(s.end());
    }
    
    return l.sqr_dist(s.start());
  }
  
  leda_rational operator()(const leda_rat_segment& s, const leda_rat_line& l) const
  {
    return this->operator()(l,s);   
  }
 
  leda_rational operator()(const leda_rat_line& l, const leda_rat_triangle& t) const
  {
    bool inter = t.intersection(l);
    if (inter) return 0;
    
    // no intersection; compare dists of the three points of t to l
    int cm = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(),l.point2(),t.point1(),t.point2());
    
    if (cm >= 0) { //point2 has smaller or equal distance to l ...
      cm = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(),l.point2(),t.point2(),t.point3());

      if (cm >= 0) { //point3 has smaller or equal distance to l ...
        return l.sqr_dist(t.point3());
      }    
    
      return l.sqr_dist(t.point2());    
    }
    
    //point1 has smaller or equal distance (than point2) to l ...
    cm = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(),l.point2(),t.point1(),t.point3());

    if (cm >= 0) { //point3 has smaller or equal distance to l ...
      return l.sqr_dist(t.point3());
    }    
    
    return l.sqr_dist(t.point1());
  }

  leda_rational operator()(const leda_rat_triangle& t, const leda_rat_line& l) const
  {
    return this->operator()(l,t);   
  }    
  
  // ray ...
  leda_rational operator()(const leda_rat_ray& r1, const leda_rat_ray& r2) const
  {
    // intersection ??
    leda_rat_point pi;
    
    bool b = r1.intersection(r2, pi);
    
    if (b) return 0;
    
    // no intersection ...
    leda_rat_point s1 = r1.source();
    leda_rat_point s2 = r2.source();
    
    // take shorter distance: s1 from r2 or s2 from r1 ...
    leda_rational d1 = this->operator()(s1,r2);
    leda_rational d2 = this->operator()(s2,r1);
    
    if (d1 < d2) return d1;
    return d2;
  }
  
  leda_rational operator()(const leda_rat_ray& r, const leda_rat_segment& s) const
  {
    // intersection ??
    leda_rat_point pi;
    
    bool b = r.intersection(s, pi);
    
    if (b) return 0;  
    
    // no intersection ...
    leda_rat_point s1 = s.start();
    leda_rat_point s2 = s.end();
    leda_rat_point rs = r.source();
    
    // take shorter distance: s1 from r or s2 from r or rs from s ...
    leda_rational d1 = this->operator()(s1,r);
    leda_rational d2 = this->operator()(s2,r); 
    leda_rational d3 = this->operator()(rs,s);  
    
    if (d1 < d2) {
      if (d1 < d3) return d1;
      return d3;
    }
    
    if (d2 < d3) return d1;
    return d3;               
  }
  
  leda_rational operator()(const leda_rat_segment& s, const leda_rat_ray& r) const
  {
    return this->operator()(r,s);
  } 
  
  leda_rational operator()(const leda_rat_ray& r, const leda_rat_triangle& t) const
  {
    // does the ray intersect one of the triangle sides ??
    leda_rat_point pi;
    leda_rat_segment s1(t.point1(), t.point2());
    
    bool b = r.intersection(s1, pi);
    if (b) return 0;
    
    leda_rat_segment s2(t.point2(), t.point3());
    
    b = r.intersection(s2, pi);
    if (b) return 0;    
    
    leda_rat_segment s3(t.point3(), t.point1());
    
    b = r.intersection(s3, pi);
    if (b) return 0;    
    
    // no intersection ...
    leda_rational d1 = this->operator()(t.point1(),r);
    leda_rational d2 = this->operator()(t.point2(),r); 
    if (d2 < d1) d1 = d2;
    leda_rational d3 = this->operator()(t.point3(),r);
    if (d3 < d1) d1 = d3;
    leda_rational d4 = this->operator()(r.source(),t);  
    if (d4 < d1) return d4;
    return d1;
  }    
  
  leda_rational operator()(const leda_rat_triangle& t, const leda_rat_ray& r) const
  {
    return this->operator()(r,t);
  }    
  
  // segment ...
  leda_rational operator()(const leda_rat_segment& s1, const leda_rat_segment& s2) const
  {
    bool b = s1.intersection(s2);
    if (b) return 0;
    
    leda_rat_point p1 = s1.source();
    leda_rat_point p2 = s1.target();
    leda_rat_point p3 = s2.source();
    leda_rat_point p4 = s2.target();    
    
    // no intersection ...
    leda_rational d1 = s1.sqr_dist(p3);
    leda_rational d2 = s1.sqr_dist(p4); 
    if (d2 < d1) d1 = d2;
    leda_rational d3 = s2.sqr_dist(p1);
    if (d3 < d1) d1 = d3;
    leda_rational d4 = s2.sqr_dist(p2);  
    if (d4 < d1) return d4;
    return d1;    
  }
  
  leda_rational operator()(const leda_rat_segment& s, const leda_rat_triangle& t) const
  {
    bool b = t.intersection(s);
    if (b) return 0;
    
    // no intersection ...
    leda_rat_segment side1(t.point1(),t.point2());
    leda_rat_segment side2(t.point2(),t.point3());        
    leda_rat_segment side3(t.point3(),t.point1());
    
    leda_rational d1 = this->operator()(s,side1);
    leda_rational d2 = this->operator()(s,side2);
    if (d2 < d1) d1=d2;
    leda_rational d3 = this->operator()(s,side3);  
    if (d3 < d1) return d3;
    return d1;                  
  } 
  
  leda_rational operator()(const leda_rat_triangle& t, const leda_rat_segment& s) const
  {
    return this->operator()(s,t);
  }            
  
  
  // triangle
  // a1-a2-a3 and b1-b2-b3 form non-intersecting triangles ...
  leda_rational sq_dist(const leda_rat_point& a1, const leda_rat_point& a2, const leda_rat_point& a3,
                        const leda_rat_point& b1, const leda_rat_point& b2, const leda_rat_point& b3) const
  {
    leda_rat_segment S1(a1,a2);
  
    leda_rational d1 = S1.sqr_dist(b1);
    leda_rational d2 = S1.sqr_dist(b2);
    if (d2 < d1) d1=d2;
    leda_rational d3 = S1.sqr_dist(b3); 
    if (d3 < d1) d1=d3;   
  
    leda_rat_segment S2(a2,a3);
  
    leda_rational d4 = S2.sqr_dist(b1);
    if (d4 < d1) d1=d4;
    leda_rational d5 = S2.sqr_dist(b2);
    if (d5 < d1) d1=d5;
    leda_rational d6 = S2.sqr_dist(b3);
    if (d6 < d1) d1=d6;  
  
    leda_rat_segment S3(a3,a1);
  
    leda_rational d7 = S3.sqr_dist(b1);
    if (d7 < d1) d1=d7;
    leda_rational d8 = S3.sqr_dist(b2);
    if (d8 < d1) d1=d8;
    leda_rational d9 = S3.sqr_dist(b3);     
    if (d9 < d1) d1=d9; 
  
    leda_rat_segment S4(b1,b2);
  
    leda_rational e1 = S4.sqr_dist(a1);
    if (e1 < d1) d1=e1;
    leda_rational e2 = S4.sqr_dist(a2);
    if (e2 < d1) d1=e2;
    leda_rational e3 = S4.sqr_dist(a3); 
    if (e3 < d1) d1=e3;   
  
    leda_rat_segment S5(b2,b3);
  
    leda_rational e4 = S5.sqr_dist(a1);
    if (e4 < d1) d1=e4;
    leda_rational e5 = S5.sqr_dist(a2);
    if (e5 < d1) d1=e5;
    leda_rational e6 = S5.sqr_dist(a3);
    if (e6 < d1) d1=e6;  
  
    leda_rat_segment S6(b3,b1);
  
    leda_rational e7 = S6.sqr_dist(a1);
    if (e7 < d1) d1=e7;
    leda_rational e8 = S6.sqr_dist(a2);
    if (e8 < d1) d1=e8;
    leda_rational e9 = S6.sqr_dist(a3);     
    if (e9 < d1) d1=e9; 
  
    return d1;  
  }  
  
  leda_rational operator()(const leda_rat_triangle& t1, const leda_rat_triangle& t2) const
  {
    int ori = t1.orientation();
 
    // ori == 0 : the triangle t1 is degenerated to a segment ...
    if (ori == 0) {
     leda_rat_segment seg1(t1.point1(),t1.point2()), seg2(t1.point2(),t1.point3());
 
     if (t2.intersection(seg1) || t2.intersection(seg2)) return 0;
   
     return this->sq_dist(t1.point1(),t1.point2(),t1.point3(),t2.point1(),t2.point2(),t2.point3());
    }

    leda_rat_point a1, a2, a3;
    leda_rat_point b1, b2, b3;

    // positive ori of t1
   if (ori == 1){
     a1 = t1.point1();
     a2 = t1.point2();
     a3 = t1.point3();
   }
   else {
     a3 = t1.point1();
     a2 = t1.point2();
     a1 = t1.point3(); 
   }
   // a1-a2-a3 are counterclockwise ...
 
   b1 = t2.point1(); b2 = t2.point2(); b3 = t2.point3();
 
   int ori1 = LEDA_NAMESPACE_NAME::orientation(a1,a2,b1);
   int ori2 = LEDA_NAMESPACE_NAME::orientation(a1,a2,b2);
   int ori3 = LEDA_NAMESPACE_NAME::orientation(a1,a2,b3);
 
   if (ori1==-1 && ori2==-1 && ori3==-1) { //t2 completely in other halfplane
    // outside ...
    return this->sq_dist(a1,a2,a3,b1,b2,b3);
   }
 
   int ori4 = LEDA_NAMESPACE_NAME::orientation(a2,a3,b1);
   int ori5 = LEDA_NAMESPACE_NAME::orientation(a2,a3,b2);
   int ori6 = LEDA_NAMESPACE_NAME::orientation(a2,a3,b3); 
 
   if (ori4==-1 && ori5==-1 && ori6==-1) { //t2 completely in other halfplane 
    // outside ...
    return this->sq_dist(a1,a2,a3,b1,b2,b3); 
   }

   int ori7 = LEDA_NAMESPACE_NAME::orientation(a3,a1,b1);
   int ori8 = LEDA_NAMESPACE_NAME::orientation(a3,a1,b2);
   int ori9 = LEDA_NAMESPACE_NAME::orientation(a3,a1,b3); 
 
   if (ori7==-1 && ori8==-1 && ori9==-1) { //t2 completely in other halfplane 
    // outside ...
    return this->sq_dist(a1,a2,a3,b1,b2,b3);
   } 
 
   // is at least one point of t2 inside ?
   if (ori1!=-1 && ori4!=-1 && ori7!=-1) return 0;
   if (ori2!=-1 && ori5!=-1 && ori8!=-1) return 0;
   if (ori3!=-1 && ori6!=-1 && ori9!=-1) return 0;
 
   // all points of t2 outside t1 ...
   // intersection tests with sides of t1 ...
   if (! (ori1==1 && ori2==1 && ori3==1)) { if (t2.intersection(leda_rat_segment(a1,a2))) return 0; }
   if (! (ori4==1 && ori5==1 && ori6==1)) { if (t2.intersection(leda_rat_segment(a2,a3))) return 0; }
   if (! (ori7==1 && ori8==1 && ori9==1)) { if (t2.intersection(leda_rat_segment(a3,a1))) return 0; }
 
   return this->sq_dist(a1,a2,a3,b1,b2,b3);
  }         
};

*/

template<class K>
class Compute_leda_rat_squared_length_2 {

  typedef typename K::Segment_2   Segment_2;
  typedef typename K::FT          FT;

public:
  typedef Arity_tag< 1 > Arity;
  typedef FT             result_type;

  FT operator()(const Segment_2& s) const
  {
    return s.sqr_length();  
  }
};

template<class K>
class Compute_leda_rat_squared_radius_2 {

  typedef typename K::Point_2     Point_2;
  typedef typename K::FT          FT;

public:
  typedef FT             result_type;
  typedef Arity_tag< 1 > Arity;  

#if defined(CGAL_COMPATIBLE_CIRCLES)
  FT operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& C) const
  {
    return C.sqr_radius();
  }
#else
  // use LEDA circles ... 

  FT operator()(const leda_rat_circle& C) const
  {
    return C.sqr_radius();
  }
#endif  
  
  FT operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
    leda_rat_circle C(p1,p2,p3);
    return C.sqr_radius();  
  }  
  
  //this is not needed for the general kernel traits, but for the
  //Alpha shape traits ...
  FT operator()(const Point_2& p1, const Point_2& p2) const
  {
    // compute squared radius of circle with diameter (p1,p2)
    Point_2 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
    
    return p1.sqr_dist(m);
  }
};

template<class K>
class Compute_leda_rat_area_2 {

  typedef typename K::Triangle_2  Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT          FT;

public:
  typedef Arity_tag< 1 >          Arity;
  typedef FT                      result_type;

  FT operator()(const Iso_rectangle_2& r) const
  {
    return r.area();
  }
  
  FT operator()(const Triangle_2& t) const
  {
    return t.area();  
  }  
};


CGAL_END_NAMESPACE

#endif




