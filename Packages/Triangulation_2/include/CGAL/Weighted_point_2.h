// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : $RCSfile$
// source        : $Source$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_WEIGHTED_POINT_2_H
#define CGAL_WEIGHTED_POINT_2_H

#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

template < class Pt, class We >
class CGAL_Weighted_point_2 : public Pt
{
public:
  typedef We Weight;
  typedef Pt Point;
  typedef typename  Point::RT RT;
private:
	Weight _weight;

public: //constructors and destructors
  CGAL_Weighted_point_2 ()
    : Point ()
  {
    _weight=Weight ( 0 ) ;
  }


  CGAL_Weighted_point_2	( const CGAL_Weighted_point_2 &p0)
  {
    Point::operator=(p0.point() );
    _weight=Weight( p0.weight() );
  }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
	template <class Weight_0 >
	CGAL_Weighted_point_2
	( const CGAL_Weighted_point_2<  Point, Weight_0 > &p0) : Point(p0.point())
	{	_weight=Weight( p0.weight() );
	}
#endif

	CGAL_Weighted_point_2 ( const Point &p )
		: Point ( p )
	{	_weight=Weight ( 0 ) ;
	}

	CGAL_Weighted_point_2 ( const RT &x, const RT &y )
		: Point ( x, y )
	{	_weight=Weight ( 0 ) ;
	}

	CGAL_Weighted_point_2 ( const Point &p, const Weight &_weight_ )
		: Point ( p )
	{	_weight=_weight_;
	}


public:
  // conversion from Weighted_point to Weight
	operator Weight() const
	{return weight();}

	Point point() const
	{	return (Point)*this;
	}

	Weight weight() const
	{	return _weight;
	}

	Weight power(const Point &p)
	{	return ((p.x()-x())*(p.x()-x())+(p.y()-y())*(p.y()-y())-weight()*weight());
	}

};


template < class Point, class Weight >
ostream &operator<<(ostream &os, const CGAL_Weighted_point_2<Point,Weight> &p)
{
	return os << p.point() << " " << p.weight() ;
}

template < class Point, class Weight >
istream  &operator>>(istream  &is,  CGAL_Weighted_point_2<Point,Weight> &p)
{
	Weight _weight_;
	Point _point_;
	is >> _point_  >> _weight_ ;
	p=CGAL_Weighted_point_2<Point,Weight>( _point_, _weight_ );
	return is;
}




#endif
