#ifndef CEP_LEDA_RAT_INTERSECTIONS_2_H
#define CEP_LEDA_RAT_INTERSECTIONS_2_H

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
#include <LEDA/rat_geo_alg.h>

// undefine the LEDA vector and list definition
// (this was only present in older LEDA versions)
#if defined(vector)
#undef vector
#undef list
#endif

#include <vector>


#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE

class Assign_leda_rat_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

  template<class T>
  bool operator()(T& t, CGAL::Object& obj) const
  { 
    return CGAL::assign(t, obj);
  }
};

// computations of 2d intersections ...
// attention - we need the correct conversion for special segments ...

template<class HELP_KERNEL>
class CGAL_intersect_leda_rat_2 {
public:
  typedef Arity_tag< 2 >     Arity;
  typedef CGAL::Object       result_type;
  
  typedef typename HELP_KERNEL::Intersect_2  Intersect_2;
  
  // for segments (uses LEDA directly) ...
  CGAL::Object operator()(const leda_rat_segment& s1, const leda_rat_segment& s2) const
  {
    CGAL::Object obj;    
    leda_rat_segment result;    
    bool bi = s1.intersection(s2, result);    
    if (bi) {
      //result.normalize();
    
      if (result.start() == result.end()) {
        obj = CGAL::make_object(result.start()); // only a point ...
      }
      else obj = CGAL::make_object(result); // segment ...
    }    
    return obj;
  }     
  
  // this one uses the CGAL kernel ...
  
  template<class T1, class T2>
  CGAL::Object operator()(const T1& obj1, const T2& obj2) const
  {
     leda_to_cgal_2 conv;
     Intersect_2 inter;
     
     // result "stores" a CGAL object ...
     CGAL::Object result = inter(conv(obj1), conv(obj2));
     
     // convert it back ...
     // attention !
     // at the moment this might cause trouble for "special" types
     // like rat_segment_cache
     
     cgal_to_leda_2 conv_back;
     
     return conv_back(result);
  }  
};



/*
class Intersect_leda_rat_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef CGAL::Object       result_type;

  // ------------------------
  // operators for line :
  // ------------------------

  CGAL::Object operator()(const leda_rat_line& l1, const leda_rat_line& l2) const
  {
    // result is a line ?
    
    if (LEDA_NAMESPACE_NAME::equal_as_sets(l1,l2)){
      return CGAL::make_object(l1);
    }
    
    // result is a point ?
    leda_rat_point pi;
    bool inter = l1.intersection(l2,pi);
    
    if (inter) return CGAL::make_object(pi);
    
    // lines do not intersect ...
    CGAL::Object obj;
    return obj;
  }
  
  CGAL::Object operator()(const leda_rat_line& l, const leda_rat_segment& s) const
  {
    // result is segment ?
    
    leda_rat_point p1 = s.start();
    leda_rat_point p2 = s.end();
    
    if ((l.side_of(p1) == 0) && (l.side_of(p2) == 0)){
      return CGAL::make_object(s);
    }
    
    // result is a point ?
    leda_rat_point pi;
    bool inter = l.intersection(s,pi);    
    
    if (inter) return CGAL::make_object(pi);
    
    //  no intersection ...
    CGAL::Object obj;
    return obj;    
  }
  
  CGAL::Object operator()(const leda_rat_segment& s,const leda_rat_line& l) const
  {
    return this->operator()(l,s);   
  } 
  
  // warum dafuer nichts in LEDA ????
  
  CGAL::Object operator()(const leda_rat_line& l, const leda_rat_ray& r) const
  {
    leda_rat_line lhelp(r);
    
    // result is a ray ?
    
    if (LEDA_NAMESPACE_NAME::equal_as_sets(l,lhelp)){
      return CGAL::make_object(r);
    }
    
    // result is a point ?
    leda_rat_point pi;
    bool inter = l.intersection(lhelp,pi);
    
    if (inter) {
      // is the intersection really with the ray ???
      if (r.contains(pi)) return CGAL::make_object(pi);
    }
    
    // no intersection ...
    CGAL::Object obj;
    return obj;        
  }  
  
  CGAL::Object operator()(const leda_rat_ray& r,const leda_rat_line& l) const
  {
    return this->operator()(l,r);   
  }   
  
  // warum dafuer nichts in LEDA ????
  
  CGAL::Object operator()(const leda_rat_line& l, const leda_rat_triangle& t) const
  {
    // is there an intersection at all ???
    if (! t.intersection(l)) {
      // no intersection ...
      CGAL::Object obj;
      return obj;      
    }
  
    leda_rat_point p1 = t.point1();
    leda_rat_point p2 = t.point2();
    leda_rat_point p3 = t.point3();
    
    // special handling of degenerate triangles
    if (t.is_degenerate()){ // triangle is a point or segment !!!!
      if (p1==p2 && p2==p3) { // point ...
        return CGAL::make_object(p1);
      }
      // segment !!!
      // soetwas sollte im LEDA kernel sein ...
      leda_rat_point smallest, largest;
	
      //get point with smallest x - coordinate ...
      int res = leda_rat_point::cmp_x(p1,p2);
	
      if (res == 0){ // use y -coords ...
	   res = leda_rat_point::cmp_y(p1,p2);
	   if (res==-1) { smallest = p1; largest = p2; }
	   else { smallest = p2; largest = p1; }
	   
	   res = leda_rat_point::cmp_y(p3,smallest);
	   if (res == -1) smallest = p3;
	   else { // new largest ???
	     res = leda_rat_point::cmp_y(largest,p3);
	     if (res == -1) largest = p3;	   
	   }
      }
      else {
	   if (res==-1) { smallest = p1; largest = p2; }
	   else { smallest = p2; largest = p1; }
	   
	   res = leda_rat_point::cmp_x(p3,smallest);
	   if (res == -1) smallest = p3;
	   else { // new largest ???
	     res = leda_rat_point::cmp_x(largest,p3);
	     if (res == -1) largest = p3;	   
	   }      
      }
      // the triangle is degenerated to a segment
      // but what about the intersection, is it a segment
      // or a point ?????
      leda_rat_segment sd(smallest,largest);
      
      if (l.side_of(smallest)==0 && l.side_of(largest)==0) {
        return CGAL::make_object(sd);
      }

      leda_rat_point pinter;
      bool bi = l.intersection(sd,pinter);
      
      CGAL::Object obj;
      if (bi) {
        return CGAL::make_object(pinter);
      }
      // no intersection ...
      return obj;
               
    } // end of handling of degenerate triangles ...
    
    leda_rat_point i1,i2,i3;
    bool b1 = l.intersection(leda_rat_segment(p1,p2),i1);
    bool b2 = l.intersection(leda_rat_segment(p2,p3),i2);
    bool b3 = l.intersection(leda_rat_segment(p3,p1),i3);
    
    // result is a segment ...
    
    if (!b1 && !b2 && !b3){ // result on one side of the triangle ???
      int ori1 = l.side_of(p1);
      int ori2 = l.side_of(p2);
      int ori3 = l.side_of(p3);
      
      if (ori1==0 && ori2==0) return CGAL::make_object(leda_rat_segment(p1,p2));
      if (ori1==0 && ori3==0) return CGAL::make_object(leda_rat_segment(p1,p3));
      if (ori2==0 && ori3==0) return CGAL::make_object(leda_rat_segment(p2,p3)); 
      
      CGAL::Object obj;
      return obj;     
    }
    
    // intersection goes through a point...
    if (b1 && b2 && b3){
      if (i1==i2) return CGAL::make_object(leda_rat_segment(i1,i3));
      if (i1==i3) return CGAL::make_object(leda_rat_segment(i1,i2)); 
      return CGAL::make_object(leda_rat_segment(i2,i3));     
    }
    
    if (b1 && b2){
      return CGAL::make_object(leda_rat_segment(i1,i2));
    }
    
    if (b1 && b3){
      return CGAL::make_object(leda_rat_segment(i1,i3));    
    }
    
    if (b2 && b3){
      return CGAL::make_object(leda_rat_segment(i2,i3));
    }
    
    // result is a point ...
    if (b1) return CGAL::make_object(i1);
    if (b2) return CGAL::make_object(i2);
    return CGAL::make_object(i3);            
  }   
  
  CGAL::Object operator()(const leda_rat_triangle& t,const leda_rat_line& l) const
  {
    return this->operator()(l,t);   
  }  

  CGAL::Object operator()(const leda_rat_line& l, const leda_rat_rectangle& r) const
  {
    leda_rat_segment clip_seg;
    bool b = r.clip(l,clip_seg);
    
    if (b) return CGAL::make_object(clip_seg);
    CGAL::Object obj;
    return obj;
  }
  
  CGAL::Object operator()(const leda_rat_rectangle& r,const leda_rat_line& l) const
  {
    return this->operator()(l,r);   
  }      

  // ------------------------
  // operators for ray :
  // ------------------------
  
  CGAL::Object operator()(const leda_rat_ray& r1, const leda_rat_ray& r2) const
  {
    leda_rat_point pi;
    bool bi = r1.intersection(r2,pi);
  
    if (bi) { 
     // we have an intersection, but what kind of result ???
     // is it a ray, a segment or a point ??
     
     int ori1 = LEDA_NAMESPACE_NAME::orientation(r1, r2.point1());
     if (ori1 != 0) {  // point !
       return CGAL::make_object(pi);
     }
     else {
       int ori2 = LEDA_NAMESPACE_NAME::orientation(r1, r2.point2());
       if (ori2 != 0) {  // point !
         return CGAL::make_object(pi);
       }
       else {  
         // this is the case that both rays are in a single line ...
         // the intersection might be a segment or a ray or a point ...

	 int cmp1 = leda_rat_point::cmp_xy(r1.point1(),r1.point2());
	 int cmp2 = leda_rat_point::cmp_xy(r2.point1(),r2.point2());
	 
	 if (cmp1==cmp2){ // result is a ray ...
	    if (r1.contains(r2.point1())) 
	    { // ray r2 ...
	      return CGAL::make_object(r2);
	    }
	    else return CGAL::make_object(r1);
	 }
	 else { 
	   // result is a point or segment ...
	   if (r1.point1() == r2.point1()) { // point
	      return CGAL::make_object(r1.point1());
	   }
	   else { // segment ...
	      return CGAL::make_object(leda_rat_segment(r1.point1(),r2.point1()));
	   }
	 }
       }       
     }
    }
    // no intersection ...
    CGAL::Object obj;
    return obj;  
  }
  
  CGAL::Object operator()(const leda_rat_ray& r, const leda_rat_segment& s) const
  {
    leda_rat_point pi;
    bool bi = r.intersection(s,pi);
  
    if (bi) { 
     // we have an intersection, but what kind of result ???
     // is it  a segment or a point ??
     
     int ori1 = LEDA_NAMESPACE_NAME::orientation(r, s.start());
     if (ori1 != 0) {  // point !
       return CGAL::make_object(pi); 
     }
     else {
       int ori2 = LEDA_NAMESPACE_NAME::orientation(r, s.end());
       if (ori2 != 0) {  // point !
         return CGAL::make_object(pi);
       }
       else {  
         // this is the case that both the ray 
	 // and the segment are in a single line ...
         // the intersection might be a segment or a point ...
         int cm   = leda_rat_point::cmp_xy(r.source(), r.point2());
         int cmp1 = leda_rat_point::cmp_xy(r.source(), s.start());
	 int cmp2 = leda_rat_point::cmp_xy(r.source(), s.end());
	 
	 if ((cmp1==cm || cmp1==0) && (cmp2==cm || cmp2==0)){
	    // segment is completely on the ray ...
	    return CGAL::make_object(s);
	 }
	 else {
	    // result only partly on the ray ...
	    // attention : we might return null-length segments ...
	    if (cmp1==cm || cmp1==0){
	      return CGAL::make_object(leda_rat_segment(r.source(), s.start()));
	    }
	    else {
	      return CGAL::make_object(leda_rat_segment(r.source(), s.end()));
	    }
	 }
       }       
     }
    }
    // no intersection ...
    CGAL::Object obj;
    return obj;  
  }  

  CGAL::Object operator()(const leda_rat_segment& s, const leda_rat_ray& r) const
  {
    return this->operator()(r,s);
  }  
  
  CGAL::Object operator()(const leda_rat_ray& r, const leda_rat_triangle& t) const
  {
    leda_rat_line l(r.point1(), r.point2()); // get supporting line ...
    CGAL::Object  obj;   // result of intersection of supp. line with triangle ... 
    CGAL::Object nothing;
    obj = this->operator()(l,t);
    
    leda_rat_segment s_result;
    leda_rat_point   p_result;
    
    if (CGAL::assign(s_result,obj)){  // segment
      leda_rat_point i1 = s_result.start();
      leda_rat_point i2 = s_result.end();
  
      bool bw1 = r.contains(i1);
      bool bw2 = r.contains(i2);
  
      if (!bw1 && !bw2) return nothing; // intersection empty ...
      if (bw1 && bw2) { return obj; }   // the resulting segment is the intersection
  
      if (bw1) 
        s_result = leda_rat_segment(r.source(),i1);
      else 
        s_result = leda_rat_segment(r.source(),i2);    
    }
    else {
      if (CGAL::assign(p_result,obj)){  // point
        if (r.contains(p_result)) return obj;
	// intersection empty ...
	return nothing;
      }
      else { // empty ...
        return obj;
      }      
    }
    
    return obj;
  }  
  
  CGAL::Object operator()(const leda_rat_triangle& t, const leda_rat_ray& r) const
  {
    return this->operator()(r,t);
  }      
  
  CGAL::Object operator()(const leda_rat_ray& r, const leda_rat_rectangle& rect) const
  {
    leda_rat_segment seg;
    bool bi = rect.clip(r, seg);
    
    CGAL::Object obj;
    
    if (bi) obj = CGAL::make_object(seg);
    
    return obj;
  }   
  
  CGAL::Object operator()(const leda_rat_rectangle& rect, const leda_rat_ray& r) const
  {
    return this->operator()(r,rect);
  } 
  
  // ------------------------
  // operators for segment :
  // ------------------------    

  CGAL::Object operator()(const leda_rat_segment& s1, const leda_rat_segment& s2) const
  {
    CGAL::Object obj;
    
    leda_rat_segment result;
    
    bool bi = s1.intersection(s2, result);
    
    if (bi) {
      if (result.start() == result.end()) obj = CGAL::make_object(result.start());
      else obj = CGAL::make_object(result);
    }
    
    return obj;
  }  
  
  CGAL::Object operator()(const leda_rat_segment& s, const leda_rat_triangle& t) const
  {
   LEDA_NAMESPACE_NAME::region_kind reg1 = t.region_of(s.start());
   LEDA_NAMESPACE_NAME::region_kind reg2 = t.region_of(s.end());
   
   CGAL::Object obj, nothing;
   
   // both inside or on triangle
   if (reg1!=LEDA_NAMESPACE_NAME::UNBOUNDED_REGION && reg2!=LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) 
   { 
     return CGAL::make_object(s); 
   }
   
   // both start and end are outside
   if (reg1==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION && reg2==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)
   {
      leda_rat_line l(s.start(),s.end());
      
      obj = this->operator()(l,t);
      
      leda_rat_segment s_result;
      leda_rat_point   p_result;
    
      if (CGAL::assign(s_result,obj)){  // segment
        leda_rat_point i1 = s_result.start();
        leda_rat_point i2 = s_result.end();
  
        bool bw1 = s.contains(i1);
        bool bw2 = s.contains(i2);
  
        if (!bw1 && !bw2) return nothing; // intersection empty ...
        if (bw1 && bw2) { return obj; }   // the resulting segment is the intersection
  
        // these cases should not occur ...
        if (bw1) 
         return CGAL::make_object(leda_rat_segment(i1,i1));
        else 
         return CGAL::make_object(leda_rat_segment(i2,i2));    
      }
      else {
        if (CGAL::assign(p_result,obj)){  // point
          if (s.contains(p_result)) return obj;
	  // intersection empty ...
	  return nothing;
        }
        else { // empty ...
          return obj;
        }      
      }      
    }
        
    // one segment end inside, one outside ...
    leda_rat_point inside, outside;
    if (reg2==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) 
    { inside=s.start(); outside=s.end(); }
    else    
    { inside=s.end(); outside=s.start(); }
   
    // get the side of the triangle, where our
    // intersecting segment leaves it ...   
    // handle case, that one is outside, one on the triangle ...
    leda_rat_point p1 = t.point1();
    leda_rat_point p2 = t.point2();
    leda_rat_point p3 = t.point3();
   
    leda_rat_point pi,pinter;
   
    bool i1 = s.intersection(leda_rat_segment(p1,p2), pi);
   
    if (i1 && pi!=inside) { pinter = pi; }
    else {
     i1 = s.intersection(leda_rat_segment(p2,p3), pi);
     
     if (i1 && pi!=inside) { pinter = pi; }
     else {
     
       i1 = s.intersection(leda_rat_segment(p3,p1),pi);
       if (i1 && pi!=inside) { pinter = pi; }
       else { pinter = inside; } 
     }
    }
   
    return CGAL::make_object(leda_rat_segment(inside,pinter));
  }

  CGAL::Object operator()(const leda_rat_triangle& t, const leda_rat_segment& s) const
  {
    return this->operator()(s,t);
  }  

  CGAL::Object operator()(const leda_rat_segment& s, const leda_rat_rectangle& rect) const
  {
    leda_rat_segment seg;
    bool bi = rect.clip(s, seg);
    
    CGAL::Object obj;
    
    if (bi) obj = CGAL::make_object(seg);
    
    return obj;
  }   
  
  CGAL::Object operator()(const leda_rat_rectangle& rect, const leda_rat_segment& s) const
  {
    return this->operator()(s,rect);
  } 

  // -----------------------------
  // operators for triangle :
  // -----------------------------

  CGAL::Object operator()(const leda_rat_triangle& t1, const leda_rat_triangle& t2) const
  {
    CGAL::Object obj,nothing;
  
    // check bounding boxes ...
    leda_rat_point t1_xmax, t1_xmin, t1_ymax, t1_ymin;
    leda_rat_point t2_xmax, t2_xmin, t2_ymax, t2_ymin;
  
    leda_rat_point a1 = t1.point1();
    leda_rat_point a2 = t1.point2();
    leda_rat_point a3 = t1.point3();
  
    // x
    if (leda_rat_point::cmp_x(a1,a2)==-1) { t1_xmin = a1; t1_xmax =a2; }
    else { t1_xmin = a2; t1_xmax =a1; }
  
    if (leda_rat_point::cmp_x(t1_xmin,a3)==-1) {
     if (leda_rat_point::cmp_x(t1_xmax,a3)==-1) t1_xmax=a3;
    }
    else t1_xmin=a3;
  
    // y
    if (leda_rat_point::cmp_y(a1,a2)==-1) { t1_ymin = a1; t1_ymax =a2; }
    else { t1_ymin = a2; t1_ymax =a1; }
  
    if (leda_rat_point::cmp_y(t1_ymin,a3)==-1) {
     if (leda_rat_point::cmp_y(t1_ymax,a3)==-1) t1_ymax=a3;
    }
    else t1_ymin=a3;

    leda_rat_point b1 = t2.point1();
    leda_rat_point b2 = t2.point2();
    leda_rat_point b3 = t2.point3();
  
    // x
    if (leda_rat_point::cmp_x(b1,b2)==-1) { t2_xmin = b1; t2_xmax =b2; }
    else { t2_xmin = b2; t2_xmax =b1; }
  
    if (leda_rat_point::cmp_x(t2_xmin,b3)==-1) {
     if (leda_rat_point::cmp_x(t2_xmax,b3)==-1) t2_xmax=b3;
    }
    else t2_xmin=b3;
  
    // y
    if (leda_rat_point::cmp_y(b1,b2)==-1) { t2_ymin = b1; t2_ymax =b2; }
    else { t2_ymin = b2; t2_ymax =b1; }
  
    if (leda_rat_point::cmp_y(t2_ymin,b3)==-1) {
     if (leda_rat_point::cmp_y(t2_ymax,b3)==-1) t2_ymax=b3;
    }
    else t2_ymin=b3;

    // now we have the bounding boxes; compare them !
    if (leda_rat_point::cmp_x(t1_xmax,t2_xmin) == -1) return obj;
    if (leda_rat_point::cmp_x(t1_xmin,t2_xmax) == 1) return obj;
    if (leda_rat_point::cmp_y(t1_ymax,t2_ymin) == -1) return obj;
    if (leda_rat_point::cmp_y(t1_ymin,t2_ymax) == 1) return obj;  

    //bounding boxes intersect !
    LEDA_NAMESPACE_NAME::region_kind cn1 = t1.region_of(b1);
    LEDA_NAMESPACE_NAME::region_kind cn2 = t1.region_of(b2);
    LEDA_NAMESPACE_NAME::region_kind cn3 = t1.region_of(b3);

    // t2 completely contained in t1 ?
    if ((!(cn1==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && 
        (!(cn2==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && 
	(!(cn3==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION))){
	
      leda_rat_triangle t_result;	
      if (t2.orientation() == -1) { t_result = leda_rat_triangle(b1, b2, b3); }
      else { t_result = leda_rat_triangle(b3, b2, b1); }

      return CGAL::make_object(t_result);
    }

    LEDA_NAMESPACE_NAME::region_kind dn1 = t2.region_of(a1);
    LEDA_NAMESPACE_NAME::region_kind dn2 = t2.region_of(a2);
    LEDA_NAMESPACE_NAME::region_kind dn3 = t2.region_of(a3);

    // t1 completely contained in t2 ?
    if ((!(dn1==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && 
        (!(dn2==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && 
	(!(dn3==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION))){
	
      leda_rat_triangle t_result;	
      if (t1.orientation() == -1) { t_result = leda_rat_triangle(a1,a2,a3); }
      else { t_result = leda_rat_triangle(a3,a2,a1); }

      return CGAL::make_object(t_result);
    }

    // no triangle completely in the other ...
    leda_rat_segment test_seg;
    leda_rat_point   test_pt;
    leda_list<leda_rat_point> result_pts;
    
    CGAL::Object iobj1 = this->operator()(leda_rat_segment(a1,a2),t2);
    CGAL::Object iobj2 = this->operator()(leda_rat_segment(a2,a3),t2);
    CGAL::Object iobj3 = this->operator()(leda_rat_segment(a3,a1),t2);
    
    bool inter1 = iobj1.is_empty();
    bool inter2 = iobj2.is_empty();
    bool inter3 = iobj3.is_empty();

    if (!inter1 && !inter2 && !inter3) {
     // intersection empty ...
     return nothing;
    }  

    // build resulting intersection ...
    // attention - do not forget points of t2 that are corners of the intersection

    if (inter1) {
     if (CGAL::assign(test_seg, iobj1)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj1)) result_pts.append(test_pt); 
     }
    }
    if (inter2) {
     if (CGAL::assign(test_seg, iobj2)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj2)) result_pts.append(test_pt); 
     }
    }
    if (inter3) {
     if (CGAL::assign(test_seg, iobj3)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj3)) result_pts.append(test_pt); 
     }
    }
    
    if (dn1 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a1);
    if (dn2 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a2);
    if (dn3 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a3); 
    if (cn1 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(b1);
    if (cn2 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(b2);
    if (cn3 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(b3);     
    
    // we have collected the intersection points ...
    // replace this later ...
    leda_list<leda_rat_point> convex_obj = LEDA_NAMESPACE_NAME::CONVEX_HULL(result_pts);
    
    switch (convex_obj.size()) {
       case 1: { obj = CGAL::make_object(leda_rat_point(convex_obj.head())); break ; }
       case 2: { 
          obj = CGAL::make_object(leda_rat_segment(convex_obj.head(),convex_obj[convex_obj.get_item(1)])); 
	  break ; 
       }
       case 3: { 
          obj = CGAL::make_object(leda_rat_triangle(convex_obj.head(),convex_obj[convex_obj.get_item(1)],convex_obj[convex_obj.get_item(2)])); 
	  break ; 
       }
       default: { // use a std::vector to store the result ...
          std::vector<leda_rat_point> pvec(convex_obj.size());
	  leda_rat_point piter;
	  int i = 0;
	  forall(piter,convex_obj) pvec[i++] = piter;
	  obj = CGAL::make_object(pvec);
       }
    }
    
    return obj;  
  }   

  CGAL::Object operator()(const leda_rat_triangle& Tr, const leda_rat_rectangle& Re) const
  {
    CGAL::Object obj,nothing;
  
    // check bounding boxes ...
    leda_rat_point r_xmax, r_xmin, r_ymax, r_ymin;
    leda_rat_point t_xmax, t_xmin, t_ymax, t_ymin;
  
    leda_rat_point a1 = Re.lower_left();
    leda_rat_point a2 = Re.lower_right();
    leda_rat_point a3 = Re.upper_right();
    leda_rat_point a4 = Re.upper_left();

    r_xmax = a2; r_xmin = a1; r_ymax = a3; r_ymin = a1;
  
    leda_rat_point b1 = Tr.point1();
    leda_rat_point b2 = Tr.point2();
    leda_rat_point b3 = Tr.point3();
  
    // x
    if (leda_rat_point::cmp_x(b1,b2)==-1) { t_xmin = b1; t_xmax =b2; }
    else { t_xmin = b2; t_xmax =b1; }
  
    if (leda_rat_point::cmp_x(t_xmin,b3)==-1) {
     if (leda_rat_point::cmp_x(t_xmax,b3)==-1) t_xmax=b3;
    }
    else t_xmin=b3;
  
    // y
    if (leda_rat_point::cmp_y(b1,b2)==-1) { t_ymin = b1; t_ymax =b2; }
    else { t_ymin = b2; t_ymax =b1; }
  
    if (leda_rat_point::cmp_y(t_ymin,b3)==-1) {
     if (leda_rat_point::cmp_y(t_ymax,b3)==-1) t_ymax=b3;
    }
    else t_ymin=b3;

    // now we have the bounding boxes; compare them !
    if (leda_rat_point::cmp_x(r_xmax,t_xmin) == -1) return nothing;
    if (leda_rat_point::cmp_x(r_xmin,t_xmax) == 1) return nothing;
    if (leda_rat_point::cmp_y(r_ymax,t_ymin) == -1) return nothing;
    if (leda_rat_point::cmp_y(r_ymin,t_ymax) == 1) return nothing;  

    //bounding boxes intersect !
    LEDA_NAMESPACE_NAME::region_kind cn1 = Tr.region_of(a1);
    LEDA_NAMESPACE_NAME::region_kind cn2 = Tr.region_of(a2);
    LEDA_NAMESPACE_NAME::region_kind cn3 = Tr.region_of(a3);
    LEDA_NAMESPACE_NAME::region_kind cn4 = Tr.region_of(a4);
    
    if ((!(cn1==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && (!(cn2==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && 
        (!(cn3==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION)) && (!(cn4==LEDA_NAMESPACE_NAME::UNBOUNDED_REGION))){
       return CGAL::make_object(Re);
    }

    bool dn1 = Re.inside_or_contains(b1);
    bool dn2 = Re.inside_or_contains(b2);
    bool dn3 = Re.inside_or_contains(b3);

    // t1 completely contained in t2 ?
    if (dn1 && dn2 && dn3){
      leda_rat_triangle t_result;
      if (Tr.orientation() == -1) { t_result = leda_rat_triangle(b1,b2,b3); }
      else { t_result = leda_rat_triangle(b3,b2,b1); }

      return CGAL::make_object(t_result);
    }

    // no object completely in the other one ...    
    leda_rat_segment test_seg;
    leda_rat_point   test_pt;
    leda_list<leda_rat_point> result_pts;
    
    CGAL::Object iobj1 = this->operator()(leda_rat_segment(a1,a2),Tr);
    CGAL::Object iobj2 = this->operator()(leda_rat_segment(a2,a3),Tr);
    CGAL::Object iobj3 = this->operator()(leda_rat_segment(a3,a4),Tr);
    CGAL::Object iobj4 = this->operator()(leda_rat_segment(a4,a1),Tr);
    
    bool inter1 = iobj1.is_empty();
    bool inter2 = iobj2.is_empty();
    bool inter3 = iobj3.is_empty();
    bool inter4 = iobj4.is_empty();

    if (!inter1 && !inter2 && !inter3 && !inter4) {
     // intersection empty ...
     return nothing;
    }  
    
    // build resulting intersection ...
    if (inter1) {
     if (CGAL::assign(test_seg, iobj1)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj1)) result_pts.append(test_pt); 
     }
    }
    if (inter2) {
     if (CGAL::assign(test_seg, iobj2)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj2)) result_pts.append(test_pt); 
     }
    }
    if (inter3) {
     if (CGAL::assign(test_seg, iobj3)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj3)) result_pts.append(test_pt); 
     }
    }
    if (inter4) {
     if (CGAL::assign(test_seg, iobj4)){
       if (test_seg.start() != test_seg.end()) result_pts.append(test_seg.start());
       result_pts.append(test_seg.end());     
     }
     else {
       if (CGAL::assign(test_pt, iobj4)) result_pts.append(test_pt); 
     }
    }
    
    if (dn1) result_pts.append(b1);
    if (dn2) result_pts.append(b2);
    if (dn3) result_pts.append(b3); 
    if (cn1 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a1);
    if (cn2 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a2);
    if (cn3 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a3); 
    if (cn4 != LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) result_pts.append(a4);        
    
    // we have collected the intersection points ...
    // replace this later ...
    leda_list<leda_rat_point> convex_obj = LEDA_NAMESPACE_NAME::CONVEX_HULL(result_pts);
    
    switch (convex_obj.size()) {
       case 1: { obj = CGAL::make_object(leda_rat_point(convex_obj.head())); break ; }
       case 2: { 
          obj = CGAL::make_object(leda_rat_segment(convex_obj.head(),convex_obj[convex_obj.get_item(1)])); 
	  break ; 
       }
       case 3: { 
          obj = CGAL::make_object(leda_rat_triangle(convex_obj.head(),convex_obj[convex_obj.get_item(1)],convex_obj[convex_obj.get_item(2)])); 
	  break ; 
       }
       default: { // use a std::vector to store the result ...
          std::vector<leda_rat_point> pvec(convex_obj.size());
	  leda_rat_point piter;
	  int i = 0;
	  forall(piter,convex_obj) pvec[i++] = piter;
	  obj = CGAL::make_object(pvec);
       }
    }
    
    return obj;          
  }
  
  CGAL::Object operator()(const leda_rat_rectangle& Re, const leda_rat_triangle& Tr) const
  {
    return this->operator()(Tr,Re);
  }
      
  
  // -----------------------------
  // operators for iso-rectangle :
  // -----------------------------   

  CGAL::Object operator()(const leda_rat_rectangle& r1, const leda_rat_rectangle& r2) const
  {
    CGAL::Object obj;
    
    leda_list<leda_rat_rectangle> rects = r1.intersection(r2);
    
    if (! rects.empty()){
      leda_rat_rectangle rs = rects.head();
      
      if (rs.is_point()) {
        obj = CGAL::make_object(rs.lower_left());
      }
      else {
        if (rs.is_segment()) {
	  obj = CGAL::make_object(leda_rat_segment(rs.lower_left(),rs.upper_right()));
	}
	else
         obj = CGAL::make_object(rs);
      }
    }
    
    return obj;
  } 
    
};
*/

CGAL_END_NAMESPACE

#endif




