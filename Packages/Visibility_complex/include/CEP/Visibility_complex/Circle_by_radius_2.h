#ifndef CGAL_CIRCLE_BY_RADIUS_2_H
#define CGAL_CIRCLE_BY_RADIUS_2_H

#include <CGAL/Circle_2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template <class R_>
class Circle_by_radius_2 : public R_::Circle_2
{
public:
    // -------------------------------------------------------------------------
    typedef  R_                   R;
    typedef typename R::RT        RT;
    typedef typename R::FT        FT;
    typedef typename R::Circle_2  Circle_2;
    typedef typename R::Point_2   Point_2;
    // -------------------------------------------------------------------------
    Circle_by_radius_2()
	: Circle_2() , _radius(0) { }
    Circle_by_radius_2(const Circle_2 &t)
	: Circle_2(t) , _radius(CGAL_NTS sqrt(t.squared_radius())) { }
    Circle_by_radius_2(const Point_2&center,
		       const FT &radius,
		       const Orientation &orientation)
      : Circle_2(center, radius*radius, orientation) , _radius(radius) { }
    Circle_by_radius_2(const Point_2&center,
		       const FT &radius)
      : Circle_2(center, radius*radius, COUNTERCLOCKWISE) , _radius(radius) { }
    Circle_by_radius_2(const Point_2&p, const Point_2&q, const Point_2&r)
      : Circle_2(p,q,r) , 
	_radius(CGAL_NTS sqrt(squared_radius())) { }
    Circle_by_radius_2(const Point_2& p, const Point_2& q,
		       const Orientation &orientation)
      : Circle_2(p,q,orientation) , 
	_radius(CGAL_NTS sqrt(squared_radius())) { }
    Circle_by_radius_2(const Point_2& p,
		       const Point_2& q)
      : Circle_2(p,q,COUNTERCLOCKWISE) ,
	_radius(CGAL_NTS sqrt(squared_radius())) { }
    Circle_by_radius_2(const Point_2& center,
		       const Orientation& orientation)
      : Circle_2(center,FT(0),orientation) ,
	_radius(CGAL_NTS sqrt(squared_radius())) { }
    Circle_by_radius_2(const Point_2& center)
      : Circle_2(center,FT(0),COUNTERCLOCKWISE) ,
	_radius(CGAL_NTS sqrt(squared_radius())) { }
    // -------------------------------------------------------------------------
    FT radius() const { return _radius; }
private:
    FT _radius;
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif  // CGAL_CIRCLE_BY_RADIUS_2_H
