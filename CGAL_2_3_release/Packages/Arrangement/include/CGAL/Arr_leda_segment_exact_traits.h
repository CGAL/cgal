// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-25 $
// release_date  : $CGAL_Date: 2000/07/07 $
//
// file          : include/CGAL/Arr_leda_segment_exact_traits.h
// package       : arr (1.39)
// maintainer    : Eyal Flato <flato@post.tau.ac.il>
// author(s)     : Iddo Hanniel
//                 Eyal Flato <flato@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_LEDA_SEGMENT_EXACT_TRAITS
#define CGAL_ARR_LEDA_SEGMENT_EXACT_TRAITS

#include <list>

#include <CGAL/Pm_leda_segment_exact_traits.h>

CGAL_BEGIN_NAMESPACE

#define CGAL_XT_SINGLE_POINT 1
#define	CGAL_XT_ORIGINAL_POINT 2

class Arr_leda_segment_exact_traits 
        : public Pm_leda_segment_exact_traits
{
public:
        Arr_leda_segment_exact_traits() 
                : Pm_leda_segment_exact_traits() {}

public:
  typedef Pm_leda_segment_exact_traits Base;
  
  typedef Base::Curve_status           Curve_status;
  typedef Base::Curve_point_status     Curve_point_status;
  
  typedef Base::X_curve                X_curve;
  typedef X_curve                      Curve;

  bool is_x_monotone(const Curve& cv) {return true;}
  //segments are x_monotone:
  void make_x_monotone(const Curve& cv, std::list<Curve>& l) {} 

  X_curve curve_flip(const X_curve& cv) const {
	return cv.reversal();
  }

  void curve_split(const X_curve& cv, X_curve& c1, X_curve& c2, 
                   const Point& split_pt) const
  {
    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_get_point_status(cv,split_pt)==ON_CURVE);
    // does not suit pmwx 
	//CGAL_precondition(curve_source(cv)!=split_pt);
    //CGAL_precondition(curve_target(cv)!=split_pt);
    
    c1=X_curve(cv.source(),split_pt);
    c2=X_curve(split_pt,cv.target());
  }


public:

  //returns true iff the intersection is strictly right of pt
  bool do_intersect_to_right(const X_curve& c1, const X_curve& c2,
                             const Point& pt) const 
  {
    return intersection_base(c1, c2, pt, true, false, dummy_pnt1, dummy_pnt2, 
			     dummy_int);

    // Following implementation was commented out during to the 
    // introduction of intersection_base by Eyal to speed up the traits class.
    /*	if (!c1.intersection(c2))
	return false;
	X_curve xcv;
	bool res = c1.intersection(c2, xcv);
	if (!res) return false;
    
	if (lexicographically_xy_larger(xcv.source(),pt) || 
        lexicographically_xy_larger(xcv.target(),pt))
	return true;
    
	return false;
    */
  }

  
  bool nearest_intersection_to_right(const X_curve& c1,
                                     const X_curve& c2,
                                     const Point& pt,
                                     Point& p1,
                                     Point& p2) const 
  {
    bool res = intersection_base(c1, c2, pt, true, true, p1, p2, dummy_int);
	if ((res) && (dummy_int & CGAL_XT_SINGLE_POINT))
		p2 = p1;
	return res;

    // Following implementation was commented out during to the 
    // introduction of intersection_base by Eyal to speed up the traits class.
/*    X_curve xcv;
    bool res = c1.intersection(c2, xcv);
    if (!res) return false;

    if (lexicographically_xy_larger(xcv.source(),xcv.target()))
      xcv=curve_flip(xcv);
    if (lexicographically_xy_larger(xcv.target(),pt)) {
      p2=point_normalize(xcv.target());
      if (lexicographically_xy_larger(xcv.source(),pt))
        p1=point_normalize(xcv.source());
      else
        p1=pt;
      
      return true;
    }

    return false; */
  }

#ifndef CGAL_PMWX_TRAITS_HAVE_INTERSECT_TO_LEFT

  X_curve curve_reflect_in_x_and_y (const X_curve& cv) const
  {
    X_curve reflected_cv( point_reflect_in_x_and_y( cv.source()),
			  point_reflect_in_x_and_y( cv.target()));
    return reflected_cv;
  }
      

  Point point_reflect_in_x_and_y (const Point& pt) const
  {
    Point reflected_pt( -pt.xcoord(), -pt.ycoord());
    return reflected_pt;
  }
      
#else
  bool do_intersect_to_left(const X_curve& c1, const X_curve& c2,
			    const Point& pt) const 
  {
    return intersection_base(c1, c2, pt, false, false, dummy_pnt1, dummy_pnt2,
			     dummy_int);
    /*	if (!c1.intersection(c2))
	return false;
	X_curve xcv;
	bool res = c1.intersection(c2, xcv);
	if (!res) return false;
		
	if (lexicographically_xy_smaller(xcv.source(),pt) || 
	lexicographically_xy_smaller(xcv.target(),pt))
	return true;
		
	return false;*/
  }

  bool nearest_intersection_to_left(const X_curve& c1,
                                     const X_curve& c2,
                                     const Point& pt,
                                     Point& p1,
                                     Point& p2) const 
  {
    bool res = intersection_base(c1, c2, pt, false, true, p1, p2, dummy_int);
	if ((res) && (dummy_int & CGAL_XT_SINGLE_POINT))
		p2 = p1;
	return res;
    /*X_curve xcv;
    bool res = c1.intersection(c2, xcv);
    if (!res) return false;

    if (lexicographically_xy_smaller(xcv.source(),xcv.target()))
      xcv=curve_flip(xcv);
    if (lexicographically_xy_smaller(xcv.target(),pt)) {
      p2=point_normalize(xcv.target());
      if (lexicographically_xy_smaller(xcv.source(),pt))
        p1=point_normalize(xcv.source());
      else
        p1=pt;
      
      return true;
    }

    return false;*/
  }

#endif

  bool curves_overlap(const X_curve& ca, const X_curve& cb) const {
    X_curve xcv;
    //    bool res = 
    ca.intersection(cb, xcv);
    return !(xcv.is_trivial());
  }


  // returns values in p1 and p2 only if return_intersection is true
  // if (xsect_type | CGAL_XT_SINGLE_POINT) then only p1 is returned
  bool intersection_base(const X_curve& c1, const X_curve& c2,
			 const Point& pt, bool right, bool return_intersection,
			 Point &p1, Point &p2, int &xsect_type) const 
  {
    xsect_type = 0;
    if ( c1.is_trivial() )
      { 
	if ( c2.contains(c1.source()) )
	  { 
	    if (right)
	      {
		if ( lexicographically_xy_larger(c1.source(),pt) )
		  {
		    // intersection is c1.source()
		    xsect_type = CGAL_XT_SINGLE_POINT | CGAL_XT_ORIGINAL_POINT;
		    if (return_intersection)
		      {
			p1 = c1.source();
			//p2 = p1;
		      }	
		    return true; 
		  }
	      }
	    else
	      {
		if ( lexicographically_xy_smaller(c1.source(),pt) )
		  {
		    // intersection is c1.source()
		    xsect_type = CGAL_XT_SINGLE_POINT | CGAL_XT_ORIGINAL_POINT;
		    if (return_intersection)
		      {
			p1 = c1.source();
			//p2 = p1;
		      }	
		    return true; 
		  }
	      }
	  }
	else
	  {
	    return false;
	  }
      }
	  
    if ( c2.is_trivial() )
      { 
	if ( c1.contains(c2.source()) )
	  { 
	    if (right)
	      {
		if ( lexicographically_xy_larger(c2.source(),pt) )
		  {
		    // intersection is c2.source()
		    xsect_type = CGAL_XT_SINGLE_POINT | CGAL_XT_ORIGINAL_POINT;
		    if (return_intersection)
		      {
			p1 = c2.source();
			//p2 = p1;
		      }	
		    return true; 
		  }
	      }
	    else
	      {
		if ( lexicographically_xy_smaller(c2.source(),pt) )
		  {
		    // intersection is c2.source()
		    xsect_type = CGAL_XT_SINGLE_POINT | CGAL_XT_ORIGINAL_POINT;
		    if (return_intersection)
		      {
			p1 = c2.source();
			//p2 = p1;
		      }	
		    return true; 
		  }
	      }
	  }
	else
	  {
	    return false;
	  }
      }
	  
    int o1 = ::orientation(c1,c2.start()); 
    int o2 = ::orientation(c1,c2.end());
	  
    if ( o1 == 0 && o2 == 0 )
      { 
	leda_rat_point sa = c1.source(); 
	leda_rat_point sb = c1.target();
	if ( ::compare (sa,sb) > 0 )
	  { 
	    leda_rat_point h = sa; 
	    sa = sb; 
	    sb = h; 
	  }
		  
	leda_rat_point ta = c2.source(); 
	leda_rat_point tb = c2.target();
		  
	if ( ::compare (ta,tb) > 0 )
	  { 
	    leda_rat_point h = ta; 
	    ta = tb; 
	    tb = h; 
	  }
		  
	leda_rat_point a = sa;
	if (::compare(sa,ta) < 0) 
	  a = ta;
		  
	leda_rat_point b = tb; 
	if (::compare(sb,tb) < 0) 
	  b = sb;
		  
	if ( ::compare(a,b) <= 0 )
	  { 
	    // a is left-low to b
	    if (right)
	      {
		//intersection (not to the right) is rat_segment(a,b);
		if (lexicographically_xy_larger(b,pt))
		  {
		    xsect_type = 0;
		    if (return_intersection)
		      {
			//if (b_right) 
			p2 = point_normalize(b);
			if (lexicographically_xy_larger(a,pt))
			  p1 = point_normalize(a);
			else
			  p1 = pt;
		      }	
		    return true;
		  }
	      }
	    else
	      {
		//intersection (not to the right) is rat_segment(a,b);
		if (lexicographically_xy_smaller(a,pt))
		  {
		    xsect_type = 0;
		    if (return_intersection)
		      {
			p2 = point_normalize(a);
			if (lexicographically_xy_smaller(b,pt))
			  p1 = point_normalize(b);
			else
			  p1 = pt;
		      }	
		    return true;
		  }
	      }
	  }
	return false;
      }

    int o3 = ::orientation(c2,c1.start());
    int o4 = ::orientation(c2,c1.end());
	  
    if ( o1 != o2 && o3 != o4 )
      { 
	leda_integer w  = c1.dy()*c2.dx() - c2.dy()*c1.dx();
	leda_integer m1 = c1.X2()*c1.Y1() - c1.X1()*c1.Y2();
	leda_integer m2 = c2.X2()*c2.Y1() - c2.X1()*c2.Y2();
		  
	leda_rat_point p(m2*c1.dx()-m1*c2.dx(), m2*c1.dy()-m1*c2.dy(), w);
	if (right)
	  {
	    if (lexicographically_xy_larger(p,pt) )
	      {
		//intersection is rat_segment(p,p);
		if (return_intersection)
		  {
		    xsect_type = CGAL_XT_SINGLE_POINT;
		    p1 = point_normalize(p);
		    //p2 = p1;
		  }	
		return true;
	      }
	  }
	else
	  {
	    if (lexicographically_xy_smaller(p,pt) )
	      {
		//intersection is rat_segment(p,p);
		if (return_intersection)
		  {
		    xsect_type = CGAL_XT_SINGLE_POINT;
		    p1 = point_normalize(p);
		    //p2 = p1;
		  }	
		return true;
	      }
	  }
      }

    return false;
  }




private:
  Point point_normalize(const Point &pt) const {

    leda_integer g, x, y, w;
    x = pt.X();
    y = pt.Y();
    w = pt.W();
    if (x.iszero() &&  y.iszero()) {
      //g = w;
      return Point(x,y,leda_integer(1));
    }
    else {
      g = gcd(x, y);
      g = gcd(g, w);

      return Point(x/g,y/g,w/g);
    }

  }

  // Dummies  
mutable leda_rat_point dummy_pnt1, dummy_pnt2;
mutable int dummy_int;

public:

//in future versions of LEDA these operators should be defined
// friend inline
// Window_stream& operator<<(Window_stream& os, const Point& p){
//     return os << leda_point(p.xcoordD(),p.ycoordD());
//   }
// friend inline
// Window_stream& operator<<(Window_stream& os, const X_curve& c){
//     leda_segment s(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D());
//     return os << s;
//   }


};

CGAL_END_NAMESPACE

#endif
