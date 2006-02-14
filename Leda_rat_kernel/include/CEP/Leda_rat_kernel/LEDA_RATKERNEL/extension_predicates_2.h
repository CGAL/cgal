#ifndef CEP_LEDA_RAT_EXTENSION_PREDICATES_2_H
#define CEP_LEDA_RAT_EXTENSION_PREDICATES_2_H

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>
#include <CGAL/functional_base.h>

#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>
#include <LEDA/floatf.h>
#include <cmath>

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

// this file includes some special
// predicates not present in the CGAL kernel
// they are useful for instance for Planar_maps/Arrangements

CGAL_BEGIN_NAMESPACE

template<class K>
class Predicate_leda_rat_do_intersect_to_right_2 {

typedef typename K::Point_2    __My_Point_2;
typedef typename K::Segment_2  __My_Segment_2;

public:

bool operator()(const __My_Segment_2& s1, const __My_Segment_2& s2,
                const __My_Point_2& p) const
{
  if ( s1.is_trivial() ) { // s1 is a point ...
    // first compare the potential intersection point with p, this is cheaper ...
    return ( (__My_Point_2::cmp_xy(s1.start(), p) == 1) && s2.contains(s1.source()) );
  }
  if ( s2.is_trivial() ) { // s2 is a point ...
    // first compare the potential intersection point with p, this is cheaper ...
    return ( (__My_Point_2::cmp_xy(s2.start(), p) == 1) && s1.contains(s2.source()) );
  }
  
  __My_Point_2 a = s1.source();
  __My_Point_2 b = s1.target();
  __My_Point_2 c = s2.source();
  __My_Point_2 d = s2.target();  
  
  // add check ?
  // is p in or left from the x_range of both segments s1, s2 ???
  
  // do we have an intersection ? 
#if (__LEDA__ <= 431)  
  int o1 = ::orientation(s1,c); 
  int o2 = ::orientation(s1,d);
#else   
  int o1 = s1.orientation(c); 
  int o2 = s1.orientation(d);
#endif  
  
  if (o1 == o2 && o1 != 0) return false; // we have no intersection point ...
 
#if (__LEDA__ <= 431)
  int o3 = ::orientation(s2,a);
  int o4 = ::orientation(s2,b);
#else  
  int o3 = s2.orientation(a);
  int o4 = s2.orientation(b);
#endif
  
  if (o3 == o4 && o3 != 0) return false; // we have no intersection point ...  

  if ( o1 == 0 && o2 == 0 ) {
     bool b1 = s2.contains(a);
     if (b1) if (__My_Point_2::cmp_xy(a, p) == 1) return true;
     bool b2 = s2.contains(b);
     if (b2) if (__My_Point_2::cmp_xy(b, p) == 1) return true;
     bool b3 = s1.contains(c);
     if (b3) if (__My_Point_2::cmp_xy(c, p) == 1) return true;
     bool b4 = s1.contains(d);
     if (b4) if (__My_Point_2::cmp_xy(d, p) == 1) return true;
  }
  else {
   if (o1 != o2 && o3 != o4) { // segment intersection !!!
     
    // double evaluation stage ...    
    // first compute the expression Ed

     double  s1_dyD = s1.dyD();
     double  s1_dxD = s1.dxD();
     double  s2_dxD = s2.dxD();
     double  s2_dyD = s2.dyD();
     double  axd = a.XD();
     double  bxd = b.XD();
     double  cxd = c.XD();
     double  dxd = d.XD();
     double  ayd = a.YD();
     double  byd = b.YD();
     double  cyd = c.YD();
     double  dyd = d.YD(); 
     
     double  pxd = p.XD();
     double  pwd = p.WD();    
     
     // for index computation see LEDA book page 613 ff
     //
     // all indices of input values of the expression
     // get a start index 1 (this is a bit too pessimistic)
     // rules: ind(A +/- B) = 1+ max(ind(A),ind(B))*s
     //        ind(A  *  B) = 1+ ind(A)*s + ind(B)*s^2
     //        we don't have division and we have only "float integers" -> no underflow
     // with s = 1+machine_eps0
     //
     // ind(s1_dyD*s2_dxD) < 3.5
     // ind(s2_dyD*s1_dxD) < 3.5
     // ...
     // ind(c2d * s1_dxD)  < 7.5
     // ind(c1d * s2_dxD)  < 7.5
     // ind(xd) < 9
     //
     // ind(xd*pwd) < 12
     // ind(pxd*wd) < 9
     //   ->   ind(E) < 14
     
     double  wd   = s1_dyD*s2_dxD - s2_dyD*s1_dxD;  // ind < 5
     double  c1d  = bxd*ayd       - axd*byd;        // ind < 5
     double  c2d  = dxd*cyd       - cxd*dyd;        // ind < 5
     
     double  xd   = c2d * s1_dxD - c1d * s2_dxD;    // ind < 9
     
     // toggle signs if w is negative ...
     if (wd < 0) { xd = -xd; wd = -wd; }
     
     double  Ed   = xd*pwd - pxd*wd;  // ind < 14
     
     // we could add s static test (with upper bound for eps)
     // if (Ed  > +eps_static) return true;
     // if (Ed  < -eps_static) return false; 
     
     // fabs the values ...   
     s1_dyD = std::fabs(s1_dyD); s1_dxD = std::fabs(s1_dxD);
     s2_dyD = std::fabs(s2_dyD); s2_dxD = std::fabs(s2_dxD);
     axd = std::fabs(axd); bxd = std::fabs(bxd);
     cxd = std::fabs(cxd); dxd = std::fabs(dxd);
     ayd = std::fabs(ayd); byd = std::fabs(byd);
     cyd = std::fabs(cyd); dyd = std::fabs(dyd); 
     pxd = std::fabs(pxd); // pwd always positive ...
     
     // compute mes; do this by "fabsing" the possibly negative 
     // input values and by replacing - with +
              
     double  bwd   = s1_dyD*s2_dxD + s2_dyD*s1_dxD;
     double  bc1d  = bxd*ayd       + axd*byd;
     double  bc2d  = dxd*cyd       + cxd*dyd;     
     
     double  bound_xd   = bc2d * s1_dxD + bc1d * s2_dxD;
     
     // eps = ind*mes*eps0     
     double  eps   =  14*(bound_xd*pwd + pxd*bwd)* LEDA_NAMESPACE_NAME::eps0;
         
     // compare the result of the expression (Ed) with the error bound (eps)     
     if (Ed  > +eps) return true;
     if (Ed  < -eps) return false;
   
     // we could add double eval. of y-coords, but this will occur
     // not so often in typical applications, so we don't do it         
   
    // exact evaluation stage ...
    
    bool toggle = false;
   
    leda_integer w  = s1.dy()*s2.dx() - s2.dy()*s1.dx();
    leda_integer c1 = b.X()*a.Y()     - a.X()*b.Y();
    leda_integer c2 = d.X()*c.Y()     - c.X()*d.Y();
     
    leda_integer xi = c2* s1.dx() - c1* s2.dx();
    
    if (w < 0) { w= -w; xi= -xi; toggle = true; }

    // compare x - coord with point ...
    leda_integer Ei = xi*p.W() - p.X()*w;
     
    if (Ei > 0) return true; 
    if (Ei < 0) return false;
     
    // Ei == 0; use y instead ...
     
    leda_integer yi = c2* s1.dy() - c1* s2.dy();
    if (toggle) yi = -yi;
    
    if (p.Y()*w < yi*p.W()) return true; 
     
    // lex. comparison returns false ...
   }
  }
  
  return false;
}

};


template<class K>
class Predicate_leda_rat_do_intersect_to_left_2 {

typedef typename K::Point_2    __My_Point_2;
typedef typename K::Segment_2  __My_Segment_2;


public:

bool operator()(const __My_Segment_2& s1, const __My_Segment_2& s2,
                const __My_Point_2& p) const 
{
  if ( s1.is_trivial() ) { // s1 is a point ...
    // first compare the potential intersection point with p, this is cheaper ...
    return ( (__My_Point_2::cmp_xy(s1.start(), p) == -1) && s2.contains(s1.source()) );
  }
  if ( s2.is_trivial() ) { // s2 is a point ...
    // first compare the potential intersection point with p, this is cheaper ...
    return ( (__My_Point_2::cmp_xy(s2.start(), p) == -1) && s1.contains(s2.source()) );
  }
  
  __My_Point_2 a = s1.source();
  __My_Point_2 b = s1.target();
  __My_Point_2 c = s2.source();
  __My_Point_2 d = s2.target();  
  
  // add check ?
  // is p in or left from the x_range of both segments s1, s2 ???  
  // do we have an intersection ?

#if (__LEDA__ <= 431)  
  int o1 = ::orientation(s1,c); 
  int o2 = ::orientation(s1,d);
#else   
  int o1 = s1.orientation(c); 
  int o2 = s1.orientation(d);
#endif
  
  if (o1 == o2 && o1 != 0) return false; // we have no intersection point ...

#if (__LEDA__ <= 431)  
  int o3 = ::orientation(s2,a);
  int o4 = ::orientation(s2,b);
#else  
  int o3 = s2.orientation(a);
  int o4 = s2.orientation(b);
#endif
  
  if (o3 == o4 && o3 != 0) return false; // we have no intersection point ...  

  if ( o1 == 0 && o2 == 0 ) {
     bool b1 = s2.contains(a);
     if (b1) if (__My_Point_2::cmp_xy(a, p) == -1) return true;
     bool b2 = s2.contains(b);
     if (b2) if (__My_Point_2::cmp_xy(b, p) == -1) return true;
     bool b3 = s1.contains(c);
     if (b3) if (__My_Point_2::cmp_xy(c, p) == -1) return true;
     bool b4 = s1.contains(d);
     if (b4) if (__My_Point_2::cmp_xy(d, p) == -1) return true;
  }
  else {
   if (o1 != o2 && o3 != o4) { // segment intersection !!!
   
     // double evaluation stage ...    
     // first compute the expression Ed
     double  s1_dyD = s1.dyD();
     double  s1_dxD = s1.dxD();
     double  s2_dxD = s2.dxD();
     double  s2_dyD = s2.dyD();
     double  axd = a.XD();
     double  bxd = b.XD();
     double  cxd = c.XD();
     double  dxd = d.XD();
     double  ayd = a.YD();
     double  byd = b.YD();
     double  cyd = c.YD();
     double  dyd = d.YD(); 
     
     double  pxd = p.XD();
     double  pwd = p.WD();    
     
     // for index computation see LEDA book page 613 ff
     // and do_intersect_to_right predicate
     // (there's a short explanation of the index computation)
     
     double  wd   = s1_dyD*s2_dxD - s2_dyD*s1_dxD;  // ind < 5
     double  c1d  = bxd*ayd       - axd*byd;        // ind < 5
     double  c2d  = dxd*cyd       - cxd*dyd;        // ind < 5
     
     double  xd   = c2d * s1_dxD - c1d * s2_dxD;    // ind < 9
     
     // toggle signs if w is negative ...
     if (wd < 0) { xd = -xd; wd = -wd; }     
     
     double  Ed   = xd*pwd - pxd*wd;  // ind < 14
     
     // fabs the values ...   
     s1_dyD = std::fabs(s1_dyD); s1_dxD = std::fabs(s1_dxD);
     s2_dyD = std::fabs(s2_dyD); s2_dxD = std::fabs(s2_dxD);
     axd = std::fabs(axd); bxd = std::fabs(bxd);
     cxd = std::fabs(cxd); dxd = std::fabs(dxd);
     ayd = std::fabs(ayd); byd = std::fabs(byd);
     cyd = std::fabs(cyd); dyd = std::fabs(dyd); 
     pxd = std::fabs(pxd); // pwd always positive ...
              
     double  bwd   = s1_dyD*s2_dxD + s2_dyD*s1_dxD;
     double  bc1d  = bxd*ayd       + axd*byd;
     double  bc2d  = dxd*cyd       + cxd*dyd;     
     
     double  bound_xd   = bc2d * s1_dxD + bc1d * s2_dxD;
     
     // eps = ind*mes*eps0
     
     double  eps   =  14*(bound_xd*pwd + pxd*bwd)* LEDA_NAMESPACE_NAME::eps0;
     
     if (Ed  > +eps) return false;
     if (Ed  < -eps) return true;
     
     // we could add double eval. of y-coords, but this will occur
     // not so often in typical applications, so we don't do it         
   
     // exact evaluation stage ...
     bool toggle = false;
   
     leda_integer w  = s1.dy()*s2.dx() - s2.dy()*s1.dx();
     leda_integer c1 = b.X()*a.Y()     - a.X()*b.Y();
     leda_integer c2 = d.X()*c.Y()     - c.X()*d.Y();
     
     leda_integer xi = c2* s1.dx() - c1* s2.dx();
     
     if (w < 0) { w= -w; xi= -xi; toggle = true; }
     
     // compare x - coord with point ...
     leda_integer Ei = xi*p.W() - p.X()*w;
     
     if (Ei > 0) return false; 
     if (Ei < 0) return true;
     
     // Ei == 0; use y instead ...
     
     leda_integer yi = c2* s1.dy() - c1* s2.dy();
     if (toggle) yi = -yi;
     
     if (p.Y()*w > yi*p.W()) return true; 
   }
  }
  
  return false;
}

};


CGAL_END_NAMESPACE

#endif
