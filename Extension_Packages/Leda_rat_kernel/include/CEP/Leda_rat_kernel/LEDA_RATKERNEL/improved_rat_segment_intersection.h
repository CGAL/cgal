#ifndef CEP_LEDA_RAT_SEGMENT_INTERSECTION_2_H
#define CEP_LEDA_RAT_SEGMENT_INTERSECTION_2_H

#include <LEDA/rat_segment.h>
// improved version of rat_segment intersection ...
// returns -1 : no intersection
// returns  1 : intersection is a point
// returns  2 : intersection is a segment

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

int rat_segment_intersection(const leda_rat_segment& s,
                             const leda_rat_segment& t, 
			     leda_rat_segment& I)
{ 
  if ( s.is_trivial() )
  {   
    if ( t.contains(  s.source()) )
       { I = leda_rat_segment(s.source(),s.source()); return 1; }
    else
       return -1;
  }
  if ( t.is_trivial() ) 
  { 
    if ( s.contains(t.source()) )
       { I = leda_rat_segment(t.source(),t.source()); return 1; }
    else 
       return -1;
  }

/* 
  // add a check of x ranges ...
  // (this makes it slower ...)
  leda_rat_point minp, maxp;
  if (leda_rat_point::cmp_x(s.start(),s.end()) == -1) { minp=s.start(); maxp=s.end(); }
  else { maxp=s.start(); minp=s.end(); }
  
  int res1 = leda_rat_point::cmp_x(minp,t.start());
  int res2 = leda_rat_point::cmp_x(minp,t.end());
  
  if (res1 == 1 && res2 == 1) return -1; // no intersection possible ...
  
  int res3 = leda_rat_point::cmp_x(maxp,t.start());
  int res4 = leda_rat_point::cmp_x(maxp,t.end()); 
   
  if (res3 == -1 && res4 == -1) return -1; // no intersection possible ... 
*/
   
  int o1 = s.orientation(t.start()); 
  int o2 = s.orientation(t.end());
  
  //int o1 = LEDA_NAMESPACE_NAME::orientation(s.start(), s.end(), t.start()); 
  //int o2 = LEDA_NAMESPACE_NAME::orientation(s.start(), s.end(), t.end());  
  
  // two orientation tests were moved ...

  // special case : collinearity ...
  if ( o1 == 0 && o2 == 0 )
  { leda_rat_point sa = s.source(); 
    leda_rat_point sb = s.target();
    if ( LEDA_NAMESPACE_NAME::compare(sa,sb) > 0 )
       { leda_rat_point h = sa; sa = sb; sb = h; }

    leda_rat_point ta = t.source(); 
    leda_rat_point tb = t.target();

    if ( LEDA_NAMESPACE_NAME::compare (ta,tb) > 0 )
       { leda_rat_point h = ta; ta = tb; tb = h; }

    leda_rat_point a = sa;
    if (LEDA_NAMESPACE_NAME::compare(sa,ta) < 0) a = ta;

    leda_rat_point b = tb; 
    if (LEDA_NAMESPACE_NAME::compare(sb,tb) < 0) b = sb;

    if ( LEDA_NAMESPACE_NAME::compare(a,b) <= 0 )
    { I = leda_rat_segment(a,b);
      if (a == b) return 1;
      else return 2;
    }
    return -1;
  }
  
  if ( o1 != o2){ // was && ...
   // moved ...
   int o3 = t.orientation(s.start());
   int o4 = t.orientation(s.end());  
   //int o3 = LEDA_NAMESPACE_NAME::orientation(t.start(),t.end(),s.start());
   //int o4 = LEDA_NAMESPACE_NAME::orientation(t.start(),t.end(),s.end());     
  
   if ( o3 != o4 )
   { leda_integer w  = s.dy()*t.dx() - t.dy()*s.dx();
     leda_integer c1 = s.X2()*s.Y1() - s.X1()*s.Y2();
     leda_integer c2 = t.X2()*t.Y1() - t.X1()*t.Y2();

     leda_rat_point p(c2*s.dx()-c1*t.dx(), c2*s.dy()-c1*t.dy(), w);
     I = leda_rat_segment(p,p);
     return 1;
   }
  }
  return -1;
}

#endif
