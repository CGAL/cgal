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
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Arr_circles_real_traits.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_CIRCLES_REAL_TRAITS_H
#define CGAL_ARR_CIRCLES_REAL_TRAITS_H

#include <CGAL/basic.h>
#include <list>

#include <CGAL/Cartesian.h>

#include <CGAL/Circle_2.h>
#include <CGAL/Point_2.h>
//#include <CGAL/squared_distance_2.h>

//#include <CGAL/IO/Window_stream.h>
//#include <CGAL/leda_real.h>

CGAL_BEGIN_NAMESPACE

template <class NT> class Arr_circles_real_traits;

template <class NT>
class Circ_Curve {

public:
  typedef Cartesian<NT>        R;
  typedef Point_2<R>           Point;
  typedef Circle_2<R>          Circle;

  Circ_Curve(const NT& x, const NT& y, const NT& r2) : 
    c(Point(x,y),r2), s(x-CGAL::sqrt(r2),y),t(x-CGAL::sqrt(r2),y)
  {}
  Circ_Curve(const NT& x, const NT& y, const NT& r2, 
	     const Point& src , const Point& trgt) :
    c(Point(x,y),r2), s(src),t(trgt)
  {
    CGAL_precondition(c.has_on_boundary(src));
    CGAL_precondition(c.has_on_boundary(trgt));
  }
  
  //a ctr with cgal_circle
  Circ_Curve(const Circle& _c) : c(_c), 
    s(_c.center().x()-CGAL::sqrt(c.squared_radius()), _c.center().y()),
    t(_c.center().x()-CGAL::sqrt(c.squared_radius()), _c.center().y()) {}
  
  Circ_Curve(const Circle& _c, const Point& src, const Point& trgt) :
    c(_c), s(src), t(trgt)
  { 
    CGAL_precondition(c.has_on_boundary(src));
    CGAL_precondition(c.has_on_boundary(trgt));
  }  

  Circ_Curve () {}
  Circ_Curve (const Circ_Curve& cv) : c(cv.c),s(cv.s),t(cv.t)
  {}
  Circ_Curve& operator=(const Circ_Curve& cv) {
    c=cv.c; s=cv.s; t=cv.t;
    return *this;
  }
  
  //Public member functions for the users
  const Circle& circle() const {return c;}
  const Point& source() const {return s;}
  const Point& target() const {return t;}
  
  bool is_x_monotone() const {
    if (s==t) 
      return false;  //closed circle
    int point_position=CGAL::compare_y(s,c.center()) * 
      CGAL::compare_y(t,c.center());
    if (point_position < 0)
      return false; //one above and one below
    if  (orientation(s,c.center(),t)!=(c.orientation()) )
      return true; //if the same orientation or on diameter (==COLLINEAR)
    return false;
  }
  
  friend class Arr_circles_real_traits<NT>;
  
private:
  Circle c;
  Point s,t;  //source, target
};

// DEBUG
// template <class NT> 
//   ::std::ostream& operator <<  
//   (::std::ostream& os,const Circ_Curve<NT>& cv) {  
//     os << "Curve:\n" ;
//     os << "s: " << cv.source() << std::endl;
//     os << "t: " << cv.target() << std::endl;
//     os << "circle: " << cv.circle() << std::endl;
//     return os;
//   }
// DEBUG

template <class _NT>
class Arr_circles_real_traits {
public:
  
  typedef _NT                    NT;

  //the difference between Curve and X_curve is semantical only,
  // NOT syntactical
  typedef Circ_Curve<NT>         Curve;
  typedef Curve                  X_curve;
 
  // using typename to please compiler (e.g., CC with IRIX64 on mips)
  typedef typename Curve::R      R;
  typedef typename Curve::Point  Point;
  typedef typename Curve::Circle Circle;

  typedef Vector_2<R>            Vector;

  // The workaround below seems obsolete and produces an error on gcc-3.
  typedef enum {
    UNDER_CURVE = -1,
    CURVE_NOT_IN_RANGE,
    ABOVE_CURVE,
    ON_CURVE
  } Curve_point_status;

// #ifndef __GNUC__
//   enum Curve_point_status {
//     UNDER_CURVE = -1,
//     ABOVE_CURVE = 1,
//     ON_CURVE = 2,
//     CURVE_NOT_IN_RANGE = 0
//   };
// #else
//   //workaround for egcs, otherwise we get an ICE
//   typedef int Curve_point_status;
//   static const int UNDER_CURVE = -1;
//   static const int ABOVE_CURVE = 1;
//   static const int ON_CURVE = 2;
//   static const int CURVE_NOT_IN_RANGE = 0;
// #endif

  Arr_circles_real_traits() {}


  Comparison_result compare_x(const Point& p0, const Point& p1) const {
    return compare_value(p0.x(),p1.x());
  }
  Comparison_result compare_y(const Point& p0, const Point& p1) const {
    return compare_value(p0.y(),p1.y());
  }

  //no vertical segments :
  bool curve_is_vertical(const X_curve&) const {return false;} 

  bool curve_is_in_x_range(const X_curve& cv, const Point& p) const {
    CGAL_precondition(is_x_monotone(cv));

    return (compare_x(p,cv.s) * compare_x(p,cv.t)) <=0 ;
  }

  Comparison_result curve_compare_at_x(const X_curve& cv1, 
				       const X_curve& cv2, 
				       const Point& p) const {
    CGAL_precondition(is_x_monotone(cv1));
    CGAL_precondition(is_x_monotone(cv2));

    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;
    
    Point p1=curve_calc_point(cv1,p);
    Point p2=curve_calc_point(cv2,p);

    return compare_y(p1,p2);
  }


  Comparison_result curve_compare_at_x_left(const X_curve& cva, 
					    const X_curve& cvb,
					    const Point& p) const {

    CGAL_precondition(is_x_monotone(cva));
    CGAL_precondition(is_x_monotone(cvb));

    //cases in which the curve is not defined - return EQUAL 
    if (!curve_is_in_x_range(cva,p) || !curve_is_in_x_range(cvb,p) )
      return EQUAL;

    //check that both curves are defined to the left of p
    if ( (compare_x(curve_source(cva),p)!=SMALLER) &&
         (compare_x(curve_target(cva),p)!=SMALLER) )
      return EQUAL;
    if ( (compare_x(curve_source(cvb),p)!=SMALLER) &&
         (compare_x(curve_target(cvb),p)!=SMALLER) )
      return EQUAL;

    Comparison_result r = curve_compare_at_x(cva, cvb, p);
    
    if ( r != EQUAL)
      return r;     // since the curve is continous 

    //otherwise
    // <cv1> and <cv2> meet at a point with the same x-coordinate as p
    // compare their derivatives
    Point q=curve_calc_point(cva,p);

    X_curve cv1(cva),cv2(cvb);
    if (compare_x(curve_source(cv1),q) == SMALLER)
      cv1 = curve_flip(cva);
    if (compare_x(curve_source(cv2),q) == SMALLER)
      cv2 = curve_flip(cvb);
 
    Vector d1=derivative_vec(cv1,q);
    Vector d2=derivative_vec(cv2,q);

    if ((compare_value(d1[0],0)==EQUAL)||
        (compare_value(d2[0],0)==EQUAL) ) { //one or both are infinite
      if (CGAL_NTS is_negative(d1[1]*d2[1])) {
        return compare_value(d1[1],d2[1]) ;
      }
      else {
        if (compare_value(d1[0],0)!=EQUAL) //d2 is vertical
          return compare_value(0,d2[1]);
        if (compare_value(d2[0],0)!=EQUAL) //d1 is vertical
          return compare_value(d1[1],0);

        //otherwise both are vertical
        //and they share a tangent at q
        //compare the norm of tangent vector (==second derivative)
        if (CGAL_NTS is_negative(compare_x(cv1.s,cv1.t) * 
				 cv1.c.orientation()) ) { 
          //curves are on lower part of circle (if d2 has greater value then
          //it is below d1 and return LARGER)
          return 
            compare_value(d2[0]*d2[0]+d2[1]*d2[1], d1[0]*d1[0]+d1[1]*d1[1]);
        }      
        else { //upper part of circle(if d1 has greater value then
          //it is above d2 and return LARGER)
          return
            compare_value(d1[0]*d1[0]+d1[1]*d1[1], d2[0]*d2[0]+d2[1]*d2[1]);
        }
      }
    }

    //in any other case both derivatives are finite and to the left of q
    //    return compare_value(derivative(cv2,q), derivative(cv1,q));
    Comparison_result ccr=compare_value(d2[1]/d2[0], d1[1]/d1[0] );

    if (ccr!=EQUAL)
      return ccr;
    else {
      //they share a common tangent
      //compare the second derivative (norm of vectors) - needs checking
      //check if we are above or below

      bool cv1_is_on_lower=(compare_x(cv1.s,cv1.t) * cv1.c.orientation() < 0);
      bool cv2_is_on_lower=(compare_x(cv2.s,cv2.t) * cv2.c.orientation() < 0);
       
      if (cv1_is_on_lower != cv2_is_on_lower) { 
	//one is from above one from below
	if (cv1_is_on_lower) return LARGER; 
	else
	  return SMALLER;
      }

      //otherwise both are on upper or both on lower
      if (cv1_is_on_lower) { 
	//curves are on lower part of circle (if |d2| has greater value then
	//it is below d1 and return LARGER)
	return 
	  compare_value(d2[0]*d2[0]+d2[1]*d2[1], d1[0]*d1[0]+d1[1]*d1[1]);
      }      
      else { //upper part of circle(if |d1| has greater value then
	//it is above d2 and return LARGER)
	return
	  compare_value(d1[0]*d1[0]+d1[1]*d1[1], d2[0]*d2[0]+d2[1]*d2[1]);
      }

    }

  }

  


  Comparison_result curve_compare_at_x_right(const X_curve& cva, 
					     const X_curve& cvb,
					     const Point& p) const {
    CGAL_precondition(is_x_monotone(cva));
    CGAL_precondition(is_x_monotone(cvb));

    //cases in which the curve is not defined - return EQUAL 
    if (!curve_is_in_x_range(cva,p) || !curve_is_in_x_range(cvb,p) )
      return EQUAL;
    
    //check that both curves are defined to the right of p
    if ( (compare_x(curve_source(cva),p)!=LARGER) &&
         (compare_x(curve_target(cva),p)!=LARGER) )
      return EQUAL;
    if ( (compare_x(curve_source(cvb),p)!=LARGER) &&
         (compare_x(curve_target(cvb),p)!=LARGER) )
      return EQUAL;

    Comparison_result r = curve_compare_at_x(cva, cvb, p);
    
    if ( r != EQUAL)
      return r;     // since the curve is continous 
    // <cv1> and <cv2> meet at a point with the same x-coordinate as p
    // compare their derivatives

    //make both curves compared - left to right
    Point q=curve_calc_point(cva,p);
    
    X_curve cv1(cva),cv2(cvb);
    if (compare_x(curve_source(cv1),q) == LARGER)
      cv1 = curve_flip(cva);
    if (compare_x(curve_source(cv2),q) == LARGER)
      cv2 = curve_flip(cvb);
    
    Vector d1=derivative_vec(cv1,q);
    Vector d2=derivative_vec(cv2,q);

    if ((compare_value(d1[0],0)==EQUAL)||
        (compare_value(d2[0],0)==EQUAL) ) { //one or both are vertical
      if (CGAL_NTS is_negative(d1[1]*d2[1]) ) {
        return compare_value(d1[1],d2[1]) ;
      }
      else {
        if (compare_value(d1[0],0)!=EQUAL) //d2 is vertical
          return compare_value(0,d2[1]);
        if (compare_value(d2[0],0)!=EQUAL ) //d1 is vertical
          return compare_value(d1[1],0);
        
        //otherwise they share a tangent at q
        //compare the norm of tangent vector (==second derivative)
        if (compare_x(cv1.s,cv1.t) * cv1.c.orientation() < 0) { 
          //curves are on lower part of circle (if |d2| has greater value then
          //it is below d1 and return LARGER)
          return 
            compare_value(d2[0]*d2[0]+d2[1]*d2[1], d1[0]*d1[0]+d1[1]*d1[1]);
        }      
        else { //upper part of circle(if |d1| has greater value then
          //it is above d2 and return LARGER)
          return
            compare_value(d1[0]*d1[0]+d1[1]*d1[1], d2[0]*d2[0]+d2[1]*d2[1]);
        }

      }
    }

    //in other case both derivatives are finite and to the right of q
    //return compare_value(derivative(cv1,q), derivative(cv2,q));
    Comparison_result ccr=compare_value(d1[1]/d1[0], d2[1]/d2[0] );
    if (ccr!=EQUAL)
      return ccr;
    else {
      //they share a common tangent
      //compare the second derivative (norm of vectors) - needs checking
      //check if we are above or below

      bool cv1_is_on_lower=(compare_x(cv1.s,cv1.t) * cv1.c.orientation() < 0);
      bool cv2_is_on_lower=(compare_x(cv2.s,cv2.t) * cv2.c.orientation() < 0);
       
      if (cv1_is_on_lower != cv2_is_on_lower) { 
	//one is from above one from below
	if (cv1_is_on_lower) return LARGER; 
	else
	  return SMALLER;
      }

      //otherwise both are on upper or on lower
      if (cv1_is_on_lower) { 
	//curves are on lower part of circle (if |d2| has greater value then
	//it is below d1 and return LARGER)
	return 
	  compare_value(d2[0]*d2[0]+d2[1]*d2[1], d1[0]*d1[0]+d1[1]*d1[1]);
      }      
      else { //upper part of circle(if |d1| has greater value then
	//it is above d2 and return LARGER)
	return
	  compare_value(d1[0]*d1[0]+d1[1]*d1[1], d2[0]*d2[0]+d2[1]*d2[1]);
      }

    }
  }

  Curve_point_status 
  curve_get_point_status(const X_curve &cv, const Point& p) const
  {
    CGAL_precondition(is_x_monotone(cv));

    if (!curve_is_in_x_range(cv, p))
      return CURVE_NOT_IN_RANGE;
    int res = compare_y(p, curve_calc_point(cv, p));
    if (res == SMALLER)
      return UNDER_CURVE;
    if (res == LARGER)
      return ABOVE_CURVE;
    if (res == EQUAL)
      return ON_CURVE;
    return ON_CURVE;
  }


  bool curve_is_between_cw(const X_curve& cv,const X_curve& first,
                           const X_curve& second, const Point& p) const
  {
    CGAL_precondition(is_x_monotone(cv));
    CGAL_precondition(is_x_monotone(first));
    CGAL_precondition(is_x_monotone(second));

    X_curve cv0 = first;
    X_curve cv1 = second;
    X_curve cvx = cv;

    if ( !is_same(cv0.s,p) ) cv0 = curve_flip(cv0);
    if ( !is_same(cv1.s,p) ) cv1 = curve_flip(cv1);
    if ( !is_same(cvx.s,p) ) cvx = curve_flip(cvx);

    bool cv0_is_left = (compare_x(cv0.t,cv0.s)==SMALLER);
    bool cv1_is_left = (compare_x(cv1.t,cv1.s)==SMALLER);
    bool cvx_is_left = (compare_x(cvx.t,cvx.s)==SMALLER);
    
    //4 cases - 
    if (cv0_is_left && cv1_is_left) {
      if (curve_compare_at_x_left(cv0,cv1,p)==LARGER) //cv0 above cv1
        return ( !((curve_compare_at_x_left(cv1,cvx,p)==SMALLER)&&
                   (curve_compare_at_x_left(cv0,cvx,p)==LARGER)) );
      else { //cv1 above cv0 
        //below we assume that if cvx is not defined to the left of p
        //the result is EQUAL (as defined in PM specs)
        return ((curve_compare_at_x_left(cv0,cvx,p)==SMALLER)&&
                (curve_compare_at_x_left(cv1,cvx,p)==LARGER));
      }
    }

    if (!cv0_is_left && !cv1_is_left) {
      if (curve_compare_at_x_right(cv0,cv1,p)==LARGER) //cv0 above cv1
        return ((curve_compare_at_x_right(cv1,cvx,p)==SMALLER)&&
                (curve_compare_at_x_right(cv0,cvx,p)==LARGER));
      else //cv1 above cv0
        return ( !((curve_compare_at_x_right(cv0,cvx,p)==SMALLER)&&
                   (curve_compare_at_x_right(cv1,cvx,p)==LARGER)) );
    }

    if (cv0_is_left && !cv1_is_left) {
      if (cvx_is_left)
        return (curve_compare_at_x_left(cv0,cvx,p)==SMALLER);
      else
        return (curve_compare_at_x_right(cv1,cvx,p)==SMALLER); 
    }

    if (!cv0_is_left && cv1_is_left) {
      if (cvx_is_left)
        return (curve_compare_at_x_left(cv1,cvx,p)==LARGER);
      else
        return (curve_compare_at_x_right(cv0,cvx,p)==LARGER); 
    }

    //shouldn't get here
    CGAL_assertion(false);
    return false;

  }

  bool curve_is_same(const X_curve& cv1, const X_curve& cv2) const {
    CGAL_precondition(is_x_monotone(cv1));
    CGAL_precondition(is_x_monotone(cv2));

    return (is_same( cv1.s,cv2.s) && is_same( cv1.t,cv2.t) &&
            ( cv1.c.orientation()==cv2.c.orientation()) &&
            is_same( cv1.c.center(), cv2.c.center()) &&
            compare_value( cv1.c.squared_radius(),
			   cv2.c.squared_radius()) == EQUAL);
  }

  Point curve_source(const X_curve& cv) const {
    return cv.s;
  }
  Point curve_target(const X_curve& cv) const {
    return cv.t;
  }

  Point point_to_left(const Point& p) const {
    return Point(p.x()-NT(1),p.y());
  }
  Point point_to_right(const Point& p) const {
    return Point(p.x()+NT(1),p.y());
  }
  

  ///////////////////////////////////////////////////////
  //         ARRANGEMENT FUNCS



  X_curve curve_flip(const X_curve& cv) const {
    X_curve xc(cv.c.center().x(),cv.c.center().y(),
	       cv.c.squared_radius(),cv.t, cv.s);
    xc.c=cv.c.opposite();
    return xc;
  }

  bool is_x_monotone(const Curve& cv) const {
    return cv.is_x_monotone();
  }

  void make_x_monotone(const Curve& cv, std::list<Curve>& l) const {
    // Require:
    CGAL_precondition( ! is_x_monotone(cv) );
    bool   switch_orientation  = false;
    
    // is cv a closed circle ?
    if (cv.s==cv.t) {
      // for arrangements of circles this is all that is needed
      Point src(cv.c.center().x()-CGAL::sqrt(cv.c.squared_radius()), 
		cv.c.center().y());
      Point trg(cv.c.center().x()+CGAL::sqrt(cv.c.squared_radius()), 
		cv.c.center().y());

      // bug fix, Shai, 12 Feb. 2001
      // x-monotone curves did not respect the original orintation
      typename Curve::Circle circ(cv.circle().center(), 
				  cv.circle().squared_radius(),
				  cv.circle().orientation());

      Curve top_arc(circ, src, trg);
      l.push_back(top_arc);

      Curve bottom_arc(circ, trg, src);
      l.push_back(bottom_arc);
    }
    else { 
      //if we get a curve that is not a closed circle - for completeness
      // bug fix, Shai, 12 Feb. 2001
      // curves that are split to 3 x-monotone sub curves were not handled
      const Point &center(cv.circle().center());
      Point  mid1, mid2;
      NT     sq_radius(cv.c.squared_radius());
      Curve  work_cv;
      bool   two_cuts            = false,
             left_cut_is_first   = false;
      
      // for simplicity work on CCW curve
      if (cv.c.orientation() == CLOCKWISE) {
	work_cv = curve_flip(cv);
	switch_orientation = true;
      } 
      else {
	work_cv = Curve(cv);
      }

      CGAL_assertion(work_cv.circle().orientation() == COUNTERCLOCKWISE);
      const Point &src(work_cv.source()),
	&trg(work_cv.target());

      // now we work on a CCW circular curve which is, by precondition
      // NOT x-monotone. 
      // There are four cases, denote the quadrants: II  I  
      // denote s - source, t - target               III IV
      // In two of them there is ONE spliting point, in the other two
      // there are TWO split points.

      // First, we check in which scenario we are
      if ( compare_y(src, center) == LARGER ) {
	left_cut_is_first = true;
	if ( compare_y(trg, center) == LARGER ) {
	  // s is in II, t is in I
	  two_cuts = true;
	}
	else {
	  // s is in II, t is in III or IV
	}
      }
      else {
	// source is lower then center
	if ( compare_y(trg, center) == SMALLER ) {
	  // s is in IV, t is in III
	  two_cuts = true;
	}
	else {
	  // s is in IV, t is in I or II
	}
      }

      // Second, we calculate the two or three split points
      if ( left_cut_is_first ){
	mid1 = Point(center.x() - CGAL::sqrt(sq_radius), center.y());
	if ( two_cuts ) {
	  mid2 = Point(center.x() + CGAL::sqrt(sq_radius), center.y());;
	}
	else {
	}
      }
      else {
	mid1 = Point(center.x() + CGAL::sqrt(sq_radius), center.y());
	if ( two_cuts ) {
	  mid2 = Point(center.x() - CGAL::sqrt(sq_radius), center.y());
	}
      }

      // Third, we build the split curves
      l.push_back(Curve(work_cv.circle(), src, mid1));
      if ( two_cuts ) {
	l.push_back(Curve(work_cv.circle(), mid1, mid2));
	l.push_back(Curve(work_cv.circle(), mid2, trg));
      }
      else {
	l.push_back(Curve(work_cv.circle(), mid1, trg));
      }

      // If we switched the orientation, we have to switch back
      if ( switch_orientation ) {
	for (typename std::list<Curve>::iterator lit = l.begin(); 
	     lit != l.end(); 
	     lit++) {
 	  *lit = curve_flip(*lit);
	}
      } 
    }

    // Ensure:
    // There are indeed 2 or 3 split points
    CGAL_postcondition(l.size() >= 2 && l.size() <= 3);
    // The orientations of the split curves are the same as of cv
    CGAL_postcondition_code(
			    if ( switch_orientation ) l.reverse();
			    Orientation cv_or = cv.circle().orientation();
			    typename std::list<Curve>::iterator lit;
			    typename std::list<Curve>::iterator next_it; );
    // Check consistency of end points
    CGAL_postcondition( l.begin()->source() == cv.source() );
    CGAL_postcondition_code( lit = l.end(); lit--; );
    CGAL_postcondition( lit->target() == cv.target() );

    CGAL_postcondition_code(//for all x-monotone parts
			    for(lit = l.begin();
				lit != l.end();
				lit++) {
			      next_it = lit; next_it++;  );

    CGAL_postcondition( lit->circle().orientation()  == cv_or );
    // Consistency of split points
    CGAL_postcondition( next_it == l.end() || 
			lit->target() == next_it->source() );
    // Split points are on circle
    CGAL_postcondition( cv.circle().has_on_boundary(lit->target()) );
    // parts are indeed x-monotone
    CGAL_postcondition( is_x_monotone(*lit) );
    CGAL_postcondition_code( } ); // end of for
  }

    

  void curve_split(const X_curve& cv, X_curve& c1, X_curve& c2, 
                   const Point& split_pt) const {
    CGAL_precondition(is_x_monotone(cv));

    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_get_point_status(cv,split_pt)==ON_CURVE);
    CGAL_precondition(compare_x(curve_source(cv),split_pt)!=EQUAL);
    CGAL_precondition(compare_x(curve_target(cv),split_pt)!=EQUAL);

    c1=cv;
    c2=cv;
    c1.t=split_pt;
    c2.s=split_pt;

  }


  //returns true iff the intersection is lexicographically strictly right 
  //of pt.
  bool do_intersect_to_right(const X_curve& c1, const X_curve& c2,
                             const Point& pt) const 
  {
    CGAL_precondition(is_x_monotone(c1));
    CGAL_precondition(is_x_monotone(c2));

    //two arcs from the same circle
    if ( is_same(c1.c.center(),c2.c.center()) &&
         compare_value(c1.c.squared_radius(),c2.c.squared_radius())==EQUAL ) {
      if ((is_same(c1.s,c2.s)&&(compare_x(c1.s,pt)==LARGER))||
          (is_same(c1.t,c2.t)&&(compare_x(c1.t,pt)==LARGER))||
          (is_same(c1.s,c2.t)&&(compare_x(c1.s,pt)==LARGER))||
          (is_same(c1.t,c2.s)&&(compare_x(c1.t,pt)==LARGER)))
        return true; //meeting endpoints

      bool c1_is_on_lower=(compare_x(c1.s,c1.t) * c1.c.orientation() < 0);
      bool c2_is_on_lower=(compare_x(c2.s,c2.t) * c2.c.orientation() < 0);
      if (c1_is_on_lower!=c2_is_on_lower) return false; //the case where the 
      //endpoints meet has already been taken care of, so 
      //if they are on different parts return false

      //check overlaps of x-monotone curves
      Point leftmost1,rightmost1;
      if (compare_x(c1.s,c1.t)==SMALLER) {
        leftmost1=c1.s; rightmost1=c1.t;
      }
      else {
        leftmost1=c1.t; rightmost1=c1.s;
      }

      Point leftmost2,rightmost2;
      if (compare_x(c2.s,c2.t)==SMALLER) {
        leftmost2=c2.s; rightmost2=c2.t;
      }
      else {
        leftmost2=c2.t; rightmost2=c2.s;
      }

      if ( (compare_x(rightmost1,pt)!=LARGER) ||
           (compare_x(rightmost2,pt)!=LARGER) )
        return false; //the overlap can't be right of pt (or end at)

      //otherwise if there is an overlap it has a point right of pt
      if ( compare_x(rightmost1,leftmost2)==SMALLER ||
           compare_x(rightmost2,leftmost1)==SMALLER ) { //no overlap
        return false;
      }
       
      
      return true; //there is an overlap and it has a point right of pt!

    } //end of case of arcs from same circle

    Point first;
    Point last;

    if (!circle_intersection(c1.c,c2.c,&first,&last)) return false;

    if (compare_x(first,pt)==LARGER) {
      if ((curve_get_point_status(c1,first) == ON_CURVE) &&
          (curve_get_point_status(c2,first) == ON_CURVE) )
        return true;
    }
    
    if (compare_x(last,pt)==LARGER) {
      if ((curve_get_point_status(c1,last) == ON_CURVE) &&
          (curve_get_point_status(c2,last) == ON_CURVE) )
        return true;
    }

    
    return false;
    
  }



  bool nearest_intersection_to_right(const X_curve& c1,
				     const X_curve& c2,
				     const Point& pt,
                                     Point& p1,
                                     Point& p2) const {

    CGAL_precondition(is_x_monotone(c1));
    CGAL_precondition(is_x_monotone(c2));

    Point rgt,lft;

    //case where the arcs are from the same circle
    if ( is_same(c1.c.center(),c2.c.center()) &&
         compare_value(c1.c.squared_radius(),c2.c.squared_radius())==EQUAL ) {
      //can intersect only at endpoints       
      Point rightmost1,leftmost1;
      if (compare_x(c1.s,c1.t)==LARGER) {
        rightmost1=c1.s;leftmost1=c1.t;
      }
      else {
        rightmost1=c1.t;leftmost1=c1.s;
      }

      Point rightmost2,leftmost2;
      if (compare_x(c2.s,c2.t)==LARGER) {
        rightmost2=c2.s;leftmost2=c2.t;
      }
      else {
        rightmost2=c2.t;leftmost2=c2.s;
      }

      bool c1_is_on_lower=(compare_x(c1.s,c1.t) * c1.c.orientation() < 0);
      bool c2_is_on_lower=(compare_x(c2.s,c2.t) * c2.c.orientation() < 0);
      if (c1_is_on_lower!=c2_is_on_lower) {
        //an intersection can occure only at end points
        if (is_same(rightmost1,rightmost2)) {
          if (compare_x(rightmost1,pt)==LARGER) {
            p1=p2=rightmost1;
            return true;
          }
        }
        if (is_same(leftmost1,leftmost2)) {
          if (compare_x(leftmost1,pt)==LARGER) {
            p1=p2=leftmost1;
            return true;
          }
        }
        return false;
      }

      //now we are dealing with two x-curves on the same side of circle
      if ( (compare_x(rightmost1,pt) != LARGER) ||
           (compare_x(rightmost2,pt) != LARGER) )
        return false; //the intersection can't be right of pt
 
      //now, if there is an intersection it has a point right of pt
      if ( compare_x(rightmost1,leftmost2)==SMALLER ||
           compare_x(rightmost2,leftmost1)==SMALLER ) { //no intersection
        return false;
      }
 
      //now we know there is an intersection, find p1,p2
      //p2 is the leftmost of the 2 rightmost points
      if (compare_x(rightmost1,rightmost2)==SMALLER) {
        p2=rightmost1;
      }
      else {
        p2=rightmost2;
      }

      //p1 is the rightmost of the 2 leftmost (if it is right of pt) 
      if (compare_x(leftmost1,leftmost2)==LARGER) {
        p1=leftmost1;
      }
      else {
        p1=leftmost2;
      }
      if (compare_x(p1,pt)==SMALLER) {
	//this assumes pt is on the curve, maybe we 
	//need to have p1=point_on_curve (pt.x())...?
        p1=pt; 
      }
       
      return true;
    } //end of case where arcs come fromsame circle



    Point first;
    Point last;

    circle_intersection(c1.c,c2.c,&first,&last);

    if (compare_x(first,last)==SMALLER) {
      rgt=first;
      lft=last;
    }
    else {
      rgt=last;
      lft=first;
    }

    if (compare_x(rgt,pt)==LARGER) {
      if ((curve_get_point_status(c1,rgt) == ON_CURVE) &&
          (curve_get_point_status(c2,rgt) == ON_CURVE) ) {
        p1=p2=rgt;
        return true;
      }
    }

    if (compare_x(lft,pt)==LARGER) {
      if ((curve_get_point_status(c1,lft) == ON_CURVE) &&
          (curve_get_point_status(c2,lft) == ON_CURVE) ) {
        p1=p2=lft;
        return true;
      }
    }

    //can be done differently (the check first)
    return false;

  }

  Point point_reflect_in_x_and_y (const Point& pt) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    Point reflected_pt( -pt.hx(), -pt.hy(), pt.hw());
    return reflected_pt;
  }
      

  X_curve curve_reflect_in_x_and_y (const X_curve& cv) const
  {
    Circle circ( point_reflect_in_x_and_y (cv.circle().center()),
		 cv.circle().squared_radius(), 
		 //reflection in two axes means no change in orientation
		 cv.circle().orientation());  
    // 		 CGAL::opposite( cv.circle().orientation()));
    X_curve reflected_cv( circ,
			  point_reflect_in_x_and_y (cv.source()),
			  point_reflect_in_x_and_y (cv.target()));
    return reflected_cv;
  }


  //currently we assume that no two circles overlap (might change in future)
  bool curves_overlap(const X_curve& c1, const X_curve& c2) const {
    CGAL_precondition(is_x_monotone(c1));
    CGAL_precondition(is_x_monotone(c2));

    //case where the arcs are from the same circle (otherwise return false)
    if ( is_same(c1.c.center(),c2.c.center()) &&
         compare_value(c1.c.squared_radius(),c2.c.squared_radius())==EQUAL ) {

      bool c1_is_on_lower=(compare_x(c1.s,c1.t) * c1.c.orientation() < 0);
      bool c2_is_on_lower=(compare_x(c2.s,c2.t) * c2.c.orientation() < 0);
      if (c1_is_on_lower!=c2_is_on_lower)
        return false;
      
      //check overlaps of x-monotone curves
      Point leftmost1,rightmost1;
      if (compare_x(c1.s,c1.t)==SMALLER) {
        leftmost1=c1.s; rightmost1=c1.t;
      }
      else {
        leftmost1=c1.t; rightmost1=c1.s;
      }

      Point leftmost2,rightmost2;
      if (compare_x(c2.s,c2.t)==SMALLER) {
        leftmost2=c2.s; rightmost2=c2.t;
      }
      else {
        leftmost2=c2.t; rightmost2=c2.s;
      }

      if ( compare_x(rightmost1,leftmost2)!=LARGER ||
           compare_x(rightmost2,leftmost1)!=LARGER ) { //no overlap
        return false;
      }
      else {
        return true;
      }
    }
    

    return false; //circles don't overlap
  }


  ////////////////////////////////////////////////////////////////////
  // PRIVATE FUNCS
private:
  Comparison_result  compare_value (const NT& a, const NT& b) const {
    return CGAL_NTS compare(a,b);
  }

  
  //calculates the point on the X_curve with the same x coordinate as p
  Point curve_calc_point(const X_curve& cv, const Point& p) const {
    //simple cases 
    if (compare_x(cv.s,p)==EQUAL)
      return cv.s;
    if (compare_x(cv.t,p)==EQUAL)
      return cv.t;

    NT px(p.x());
    NT sqr = (CGAL::sqrt(cv.c.squared_radius() - 
			 (px-cv.c.center().x())*(px-cv.c.center().x()) ));
    
    Point lst1_first(px,cv.c.center().y() + sqr);
    Point lst1_last(px,cv.c.center().y() - sqr);
    

    Point p1;
    if (compare_x(cv.s,cv.t) * cv.c.orientation() < 0) { //lower part of circle
      if (compare_y(lst1_first,lst1_last) == LARGER)
        p1=lst1_last;
      else
        p1=lst1_first;
    }      
    else { //upper part of circle
      if (compare_y(lst1_first,lst1_last) == LARGER)
        p1=lst1_first;
      else
        p1=lst1_last;
    }
    
    return p1;
  }
  

  Vector derivative_vec(const X_curve& cv, const Point& p) const {
    if (cv.c.orientation()==COUNTERCLOCKWISE) { //ccw - (-y,x)
      return Vector((cv.c.center().y()-p.y()), (p.x()-cv.c.center().x())); 
    } 
    else 
      return Vector((p.y()-cv.c.center().y()), (cv.c.center().x())-p.x());
  }      

  bool is_same(const Point &p1, const Point &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }


  bool circle_intersection(const Circle& ca, const Circle& cb,
			   Point* p1, Point* p2) const {
    //function checks if the circles ca,cb intersect,
    //if they don't - returns false
    //if they do p1,p2 will hold the intersection points

    NT l2=squared_distance(ca.center(),cb.center());
    NT l=CGAL::sqrt(l2);
    
    NT ra=CGAL::sqrt(ca.squared_radius());
    NT rb=CGAL::sqrt(cb.squared_radius());
    
    if ( (compare_value(l, ra+rb) == LARGER) ||
         (compare_value(ra, l+rb) == LARGER) ||
         (compare_value(rb, l+ra) == LARGER) )
      return false;

    //x is the distance on the segment-of-centers from ca.center()
    //y is the distance from the segment-of-centers to the intersection point
    NT x = (ca.squared_radius()-cb.squared_radius()+l2) / NT(2*l);
    NT y = CGAL::sqrt(ca.squared_radius() - x*x);

    //debug
    Vector v_ab=cb.center()-ca.center();
    //Vector_2<Cartesian<NT> > v_ab=cb.center()-ca.center();

    v_ab = v_ab/(CGAL::sqrt(v_ab.x()*v_ab.x()+v_ab.y()*v_ab.y())); //normalize

    Vector v_ab_perp(-v_ab.y(),v_ab.x());

    *p1 = ca.center() + x*v_ab +  y*v_ab_perp;
    *p2 = ca.center() + x*v_ab -  y*v_ab_perp;
    
    return true;

  }

};




///////////////////////////////////////////////////////////////
// auxilary output functions


/*
  //#define CGAL_ARR_IDDO_DEBUG
  #ifdef CGAL_ARR_IDDO_DEBUG
  //debug - specialization of output for reals, to make points more readable
  //(otherwise I get these lon..g reals in printouts)
  template<class NT> 
  inline ::std::ostream& operator<<(::std::ostream& o, 
  const Point_2<Cartesian<NT> >& p)
  {
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  o << " (" << x << "," << y << ") " ;
  return o;
  }
  #endif



  #ifndef CGAL_ARR_IDDO_DEBUG
//a simple version of the windowstream operator (sufficient for X_curve)
template <class NT>
//friend
Window_stream& operator<<(Window_stream& os,
  const typename Arr_circles_real_traits<NT>::Curve &cv)
{
//This is not good enough - it assumes s and t have different x coord, 
//but for x-monotone arcs it is sufficient (and that's all I need).
//runs faster than above
double px= CGAL::to_double((cv.source().x()+cv.target().x())/2);
double R2= CGAL::to_double(cv.circle().squared_radius());
double sqr = CGAL::sqrt(R2 - 
(CGAL::to_double(px-cv.circle().center().x())*
CGAL::to_double(px-cv.circle().center().x())));
  
double py;
//under part
if ((cv.source().x()-cv.target().x()) * cv.circle().orientation() < 0) 
  py= CGAL::to_double(cv.circle().center().y())-sqr;
else
py= CGAL::to_double(cv.circle().center().y())+sqr;
  
  
os.draw_arc(leda_point(CGAL::to_double(cv.source().x()),
			     CGAL::to_double(cv.source().y())),
			     leda_point(px,py),
			     leda_point(CGAL::to_double(cv.target().x()),
					CGAL::to_double(cv.target().y())));
  
return os;
}

#else 
//CGAL_ARR_IDDO_DEBUG defined - use the complicated version for general Curves
template <class NT>
Window_stream& operator<<(Window_stream& os,
  const Arr_circles_real_traits<NT>::Curve &cv)
{
  double px,py; //middle point coordinates
  double R2= CGAL::to_double(cv.circle().squared_radius());

  //checking for X-monotone case
  //the folowing is equivelent to "if (curve is x-monotone)"
  if (cv.is_x_monotone()) {
    px= CGAL::to_double((cv.source().x()+cv.target().x()))/2;
    double sqr = CGAL::sqrt(R2 - 
			    (CGAL::to_double(px-cv.circle().center().x())*
			     CGAL::to_double(px-cv.circle().center().x())));
    if (CGAL_NTS sign(cv.source().x()-cv.target().x()) * 
	cv.circle().orientation() < 0) //under part
      py= CGAL::to_double(cv.circle().center().y())-sqr;
    else
      py= CGAL::to_double(cv.circle().center().y())+sqr;
  }
  else { //if not x-monotone the above is not good enough
    if (cv.source()==cv.target()) { //closed circle
      return os << cv.circle() ;
    }
        
    py=CGAL::to_double(cv.circle().center().y());          
    if (CGAL::compare_y(cv.source(),cv.circle().center())*
	cv.circle().orientation() >0) {
      //either s is under center and orient is cw or
      //s is above and orient is ccw
      px=CGAL::to_double(cv.circle().center().x())-CGAL::sqrt(R2);
    }
    else
      if (CGAL::compare_y(cv.source(),cv.circle().center())*
	  cv.circle().orientation() < 0) {
	//either s is under center and orient is ccw or
	//s is above and orient is cw
	px=CGAL::to_double(cv.circle().center().x())+CGAL::sqrt(R2);
      }
      else 
	{ //s is one of the endpoints of the circle choos other endpoint
	  if (CGAL::compare_x(cv.source(),cv.circle().center())==SMALLER)
	    px=CGAL::to_double(cv.circle().center().x())+CGAL::sqrt(R2);
	  else
	    px=CGAL::to_double(cv.circle().center().x())-CGAL::sqrt(R2);
	}
  }

  os.draw_arc(leda_point(CGAL::to_double(cv.source().x()),
			 CGAL::to_double(cv.source().y())),
	      leda_point(px,py),
	      leda_point(CGAL::to_double(cv.target().x()),
			 CGAL::to_double(cv.target().y())));


  return os;
}
#endif //CGAL_ARR_IDDO_DEBUG
*/

CGAL_END_NAMESPACE



#endif






