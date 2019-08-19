// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>
//                 Saurabh Singh <ssingh@cs.iitr.ac.in>

#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H

#include <QRectF>
#include <vector>

// TODO: should be included in PainterOstream.h
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>

#include "Utils.h"

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

class QPainter;

namespace CGAL {
namespace Qt {

template < typename ArrTraits >
class ArrangementPainterOstreamBase : public QGraphicsSceneMixin
{
public:
  // typedefs
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Triangle_2                   Triangle_2;
  typedef typename Kernel::Iso_rectangle_2              Iso_rectangle_2;
  typedef typename Kernel::Circle_2                     Circle_2;

public:
  /*! Constructor */
  ArrangementPainterOstreamBase( QPainter* p,
                                 QRectF clippingRectangle = QRectF( ) ) :
    painterOstream( p, clippingRectangle ),
    qp( p ),
    convert( clippingRectangle ),
    // scene( NULL ),
    clippingRect( QRectF( ) ), // null rectangle
    scale( 1.0 )
  {
    if ( p != 0 )
    {
      this->scale = p->worldTransform( ).m11( );
    }
  }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstreamBase() {}

  // methods
  template < typename T >
  ArrangementPainterOstreamBase& operator<<( const T& t )
  {
    this->painterOstream << t;
    return *this;
  }

  void setScene( QGraphicsScene* scene_ )
  {
    this->scene = scene_;

    // set the clipping rectangle
    if ( scene_ == NULL )
    {
      return;
    }

    // std::cout<<"In setScene: scene_ != NULL\n";
    this->clippingRect = this->viewportRect( );
    this->convert = Converter< Kernel >( this->clippingRect );
  }

#if 0
  void setScene( QGraphicsScene* scene_ )
  {
    this->scene = scene_;

    // set the clipping rectangle
    if ( scene_ == NULL )
    {
      return;
    }
    this->clippingRect = this->getViewportRect( );
  }
#endif

protected: // methods
#if 0
  QRectF getViewportRect( ) const
  {
    // assumes scene is not null and attached to exactly one view
    QGraphicsView* view = this->scene->views( ).first( );
    QPointF p1 = view->mapToScene( 0, 0 );
    QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
    QRectF clipRect = QRectF( p1, p2 );

    return clipRect;
  }
#endif

protected:
  // fields
  PainterOstream< Kernel > painterOstream;
  QPainter* qp;
  Converter< Kernel > convert;
  // QGraphicsScene* scene;
  QRectF clippingRect;
  double scale;

}; // class ArrangementPainterOstreamBase

template < typename ArrTraits >
class ArrangementPainterOstream :
    public ArrangementPainterOstreamBase< ArrTraits >
{
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    ArrangementPainterOstreamBase< ArrTraits >( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}
};

template < typename Kernel_ >
class ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >:
  public ArrangementPainterOstreamBase<CGAL::Arr_segment_traits_2<Kernel_> >
{
public: // typedefs
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass(p, clippingRectangle)
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename SegmentTraits >
class ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> > :
  public ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                         SegmentTraits> >
{
public: // typedefs
  typedef ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                          SegmentTraits> > Superclass;
  typedef typename Superclass::Traits                   Traits;
  typedef typename Superclass::Kernel                   Kernel;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

    /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                         NtTraits > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_conic_traits_2<RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits> >
{
public: // typedefs
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Construct_x_monotone_curve_2
    Construct_x_monotone_curve_2;
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::FT                           FT;

public: // inner classes
  // utility class to use with std::sort on an Intersect_2 result set.
  class Compare_intersection_point_result
  {
  public:
    typedef std::pair< Intersection_point_2, Multiplicity > Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()( const Result& o1, const Result& o2 )
    {
      Point_2 p1 = o1.first;
      Point_2 p2 = o2.first;
      return ( p1.x( ) < p2.x( ) );
    }
  };

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle ),
    //intersect_2( this->traits.intersect_2_object( ) ),
    // Why doesn't this work?
    construct_x_monotone_curve_2(this->
                                 traits.construct_x_monotone_curve_2_object())
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    // std::cout<< "In ArrangementPainterOstream& operator T"<<std::endl;
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected: // methods
  // Returns subcurves of curve that are actually visible in the view.
  // Assumes that clippingRect is valid.
  std::vector< X_monotone_curve_2 > visibleParts( X_monotone_curve_2 curve )
  {
    // see if we intersect the bottom edge of the viewport
    Intersect_2 intersect_2 = this->traits.intersect_2_object( );
    Point_2 bottomLeft = this->convert( this->clippingRect.bottomLeft( ) );
    Point_2 bottomRight = this->convert( this->clippingRect.bottomRight( ) );
    Point_2 topLeft = this->convert( this->clippingRect.topLeft( ) );
    Point_2 topRight = this->convert( this->clippingRect.topRight( ) );
    X_monotone_curve_2 bottom =
      this->construct_x_monotone_curve_2( bottomLeft, bottomRight );
    X_monotone_curve_2 left =
      this->construct_x_monotone_curve_2( bottomLeft, topLeft );
    X_monotone_curve_2 top =
      this->construct_x_monotone_curve_2( topLeft, topRight );
    X_monotone_curve_2 right =
      this->construct_x_monotone_curve_2( topRight, bottomRight );

    std::vector< CGAL::Object > bottomIntersections;
    std::vector< CGAL::Object > leftIntersections;
    std::vector< CGAL::Object > topIntersections;
    std::vector< CGAL::Object > rightIntersections;
    std::vector< CGAL::Object > intersections;

    intersect_2( bottom, curve, std::back_inserter( bottomIntersections ) );
    intersect_2( left, curve, std::back_inserter( leftIntersections ) );
    intersect_2( top, curve, std::back_inserter( topIntersections ) );
    intersect_2( right, curve, std::back_inserter( rightIntersections ) );
    // int total = bottomIntersections.size( )
    //   + leftIntersections.size( )
    //   + topIntersections.size( )
    //   + rightIntersections.size( );

    intersect_2( bottom, curve, std::back_inserter( intersections ) );
    intersect_2( left, curve, std::back_inserter( intersections ) );
    intersect_2( top, curve, std::back_inserter( intersections ) );
    intersect_2( right, curve, std::back_inserter( intersections ) );

    this->filterIntersectionPoints( intersections );
    //std::cout << "total intersections: " << intersections.size( )
    //          << std::endl;
    //this->printIntersectResult( intersections );

    Point_2 leftEndpt = curve.source( );
    Point_2 rightEndpt = curve.target( );

    if ( leftEndpt.x( ) > rightEndpt.x( ) )
    {
      std::swap( leftEndpt, rightEndpt );
    }

    QPointF qendpt1 = this->convert( leftEndpt );
    QPointF qendpt2 = this->convert( rightEndpt );

    std::list< Point_2 > pointList;
    for ( unsigned int i = 0; i < intersections.size( ); ++i )
    {
      CGAL::Object o = intersections[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, o ) )
      {
        Point_2 pt = pair.first;
        pointList.push_back( pt );
      }
    }

    bool includeLeftEndpoint = this->clippingRect.contains( qendpt1 );
    bool includeRightEndpoint = this->clippingRect.contains( qendpt2 );
    if ( includeLeftEndpoint )
    {
      pointList.push_front( leftEndpt );
    }

    if ( includeRightEndpoint )
    {
      pointList.push_back( rightEndpt );
    }

    Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
    std::vector< X_monotone_curve_2 > clippings;
    typename std::list< Point_2 >::iterator pointListItr = pointList.begin( );
    for ( unsigned int i = 0; i < pointList.size( ); i += 2 )
    {
      Point_2 p1 = *pointListItr++;
      Point_2 p2 = *pointListItr++;
      X_monotone_curve_2 subcurve =
        construct_x_monotone_subcurve_2( curve, p1, p2 );
      clippings.push_back( subcurve );
    }

#if 0
    // std::cout << "pointList size: " << pointList.size( ) << std::endl;
    // if ( intersections.size( ) % 2 == 0 )
    // {
    //   // either both curve endpoints are in view or both are out
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         if ( this->clippingRect.contains( qendpt2 ) )
    //         {
    //             std::cout << "both endpoints are in view" << std::endl;
    //         }
    //     }
    //     else if ( !this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "both endpoints are out of view" << std::endl;
    //     }
    // }
    // else
    // { // one curve endpoint is in view
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         std::cout << "left endpoint is in view" << std::endl;
    //     }
    //     else if ( this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "right endpoint is in view" << std::endl;
    //     }
    // }

    std::vector< X_monotone_curve_2 > res;
    res.push_back( curve );
    return res;
#endif
    return clippings;
  }

  // keep only the intersection points ie. throw out overlapping curve segments
  void filterIntersectionPoints( std::vector< CGAL::Object >& res )
  {
    std::vector< std::pair< Intersection_point_2, Multiplicity > > tmp;

    // filter out the non-intersection point results
    for ( unsigned int i = 0; i < res.size( ); ++i )
    {
      CGAL::Object obj = res[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        tmp.push_back( pair );
      }
    }
    res.clear( );

    // sort the intersection points by x-coord
    Compare_intersection_point_result compare_intersection_point_result;
    std::sort( tmp.begin( ), tmp.end( ), compare_intersection_point_result );

    // box up the sorted elements
    for ( unsigned int i = 0; i < tmp.size( ); ++i )
    {
      std::pair< Intersection_point_2, Multiplicity > pair = tmp[ i ];
      CGAL::Object o = CGAL::make_object( pair );
      res.push_back( o );
    }
  }

  void printIntersectResult( const std::vector< CGAL::Object >& res )
  {
    for ( std::vector< CGAL::Object >::const_iterator it = res.begin( );
          it != res.end( ); ++it )
    {
      CGAL::Object obj = *it;
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        Point_2 pt = pair.first;
        /* QPointF qpt = */ this->convert( pt );
        // std::cout << "(" << pt.x( ) << " " << pt.y( ) < ")" << std::endl;
      }
    }
  }

protected: // members
  Traits traits;
  //Intersect_2 intersect_2;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
};

template < typename RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                         NtTraits > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_Bezier_curve_traits_2<RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits> > {

public: // typedefs
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
//  typedef typename Traits::Construct_x_monotone_curve_2
//    Construct_x_monotone_curve_2;
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::FT                           FT;

#if 0
public: // inner classes
  // utility class to use with std::sort on an Intersect_2 result set.
  class Compare_intersection_point_result
  {
  public:
    typedef std::pair< Intersection_point_2, Multiplicity > Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()( const Result& o1, const Result& o2 )
    {
      Point_2 p1 = o1.first;
      Point_2 p2 = o2.first;
      return ( p1.x( ) < p2.x( ) );
    }
  };
#endif
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle )
    { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    // std::cout<< "In ArrangementPainterOstream& operator T"<<std::endl;
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
#if 0
protected: // methods
  // Returns subcurves of curve that are actually visible in the view.
  // Assumes that clippingRect is valid.
  std::vector< X_monotone_curve_2 > visibleParts( X_monotone_curve_2 curve )
  {
    // see if we intersect the bottom edge of the viewport
    Intersect_2 intersect_2 = this->traits.intersect_2_object( );
    Point_2 bottomLeft = this->convert( this->clippingRect.bottomLeft( ) );
    Point_2 bottomRight = this->convert( this->clippingRect.bottomRight( ) );
    Point_2 topLeft = this->convert( this->clippingRect.topLeft( ) );
    Point_2 topRight = this->convert( this->clippingRect.topRight( ) );
    X_monotone_curve_2 bottom =
      this->construct_x_monotone_curve_2( bottomLeft, bottomRight );
    X_monotone_curve_2 left =
      this->construct_x_monotone_curve_2( bottomLeft, topLeft );
    X_monotone_curve_2 top =
      this->construct_x_monotone_curve_2( topLeft, topRight );
    X_monotone_curve_2 right =
      this->construct_x_monotone_curve_2( topRight, bottomRight );

    std::vector< CGAL::Object > bottomIntersections;
    std::vector< CGAL::Object > leftIntersections;
    std::vector< CGAL::Object > topIntersections;
    std::vector< CGAL::Object > rightIntersections;
    std::vector< CGAL::Object > intersections;

    intersect_2( bottom, curve, std::back_inserter( bottomIntersections ) );
    intersect_2( left, curve, std::back_inserter( leftIntersections ) );
    intersect_2( top, curve, std::back_inserter( topIntersections ) );
    intersect_2( right, curve, std::back_inserter( rightIntersections ) );
    // int total = bottomIntersections.size( )
    //   + leftIntersections.size( )
    //   + topIntersections.size( )
    //   + rightIntersections.size( );

    intersect_2( bottom, curve, std::back_inserter( intersections ) );
    intersect_2( left, curve, std::back_inserter( intersections ) );
    intersect_2( top, curve, std::back_inserter( intersections ) );
    intersect_2( right, curve, std::back_inserter( intersections ) );

    this->filterIntersectionPoints( intersections );
    //std::cout << "total intersections: " << intersections.size( )
    //          << std::endl;
    //this->printIntersectResult( intersections );

    Point_2 leftEndpt = curve.source( );
    Point_2 rightEndpt = curve.target( );

    if ( leftEndpt.x( ) > rightEndpt.x( ) )
    {
      std::swap( leftEndpt, rightEndpt );
    }

    QPointF qendpt1 = this->convert( leftEndpt );
    QPointF qendpt2 = this->convert( rightEndpt );

    std::list< Point_2 > pointList;
    for ( unsigned int i = 0; i < intersections.size( ); ++i )
    {
      CGAL::Object o = intersections[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, o ) )
      {
        Point_2 pt = pair.first;
        pointList.push_back( pt );
      }
    }

    bool includeLeftEndpoint = this->clippingRect.contains( qendpt1 );
    bool includeRightEndpoint = this->clippingRect.contains( qendpt2 );
    if ( includeLeftEndpoint )
    {
      pointList.push_front( leftEndpt );
    }

    if ( includeRightEndpoint )
    {
      pointList.push_back( rightEndpt );
    }

    Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
    std::vector< X_monotone_curve_2 > clippings;
    typename std::list< Point_2 >::iterator pointListItr = pointList.begin( );
    for ( unsigned int i = 0; i < pointList.size( ); i += 2 )
    {
      Point_2 p1 = *pointListItr++;
      Point_2 p2 = *pointListItr++;
      X_monotone_curve_2 subcurve =
        construct_x_monotone_subcurve_2( curve, p1, p2 );
      clippings.push_back( subcurve );
    }

#if 0
    // std::cout << "pointList size: " << pointList.size( ) << std::endl;
    // if ( intersections.size( ) % 2 == 0 )
    // {
    //   // either both curve endpoints are in view or both are out
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         if ( this->clippingRect.contains( qendpt2 ) )
    //         {
    //             std::cout << "both endpoints are in view" << std::endl;
    //         }
    //     }
    //     else if ( !this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "both endpoints are out of view" << std::endl;
    //     }
    // }
    // else
    // { // one curve endpoint is in view
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         std::cout << "left endpoint is in view" << std::endl;
    //     }
    //     else if ( this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "right endpoint is in view" << std::endl;
    //     }
    // }

    std::vector< X_monotone_curve_2 > res;
    res.push_back( curve );
    return res;
#endif
    return clippings;
  }

  // keep only the intersection points ie. throw out overlapping curve segments
  void filterIntersectionPoints( std::vector< CGAL::Object >& res )
  {
    std::vector< std::pair< Intersection_point_2, Multiplicity > > tmp;

    // filter out the non-intersection point results
    for ( unsigned int i = 0; i < res.size( ); ++i )
    {
      CGAL::Object obj = res[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        tmp.push_back( pair );
      }
    }
    res.clear( );

    // sort the intersection points by x-coord
    Compare_intersection_point_result compare_intersection_point_result;
    std::sort( tmp.begin( ), tmp.end( ), compare_intersection_point_result );

    // box up the sorted elements
    for ( unsigned int i = 0; i < tmp.size( ); ++i )
    {
      std::pair< Intersection_point_2, Multiplicity > pair = tmp[ i ];
      CGAL::Object o = CGAL::make_object( pair );
      res.push_back( o );
    }
  }

  void printIntersectResult( const std::vector< CGAL::Object >& res )
  {
    for ( std::vector< CGAL::Object >::const_iterator it = res.begin( );
          it != res.end( ); ++it )
    {
      CGAL::Object obj = *it;
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        Point_2 pt = pair.first;
        /* QPointF qpt = */ this->convert( pt );
        // std::cout << "(" << pt.x( ) << " " << pt.y( ) < ")" << std::endl;
      }
    }
  }
#endif
protected: // members
  Traits traits;
  //Intersect_2 intersect_2;
  //Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
};

template < typename Kernel_ >
class ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public: // typedefs
  typedef Kernel_                                       Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel >           Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename CircularKernel >
class ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<
                                   CircularKernel > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_circular_arc_traits_2<
                                          CircularKernel > >
{
public:
  typedef CircularKernel                                Kernel;
  typedef CGAL::Arr_circular_arc_traits_2< Kernel >     Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename Coefficient_ >
class ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<
                                   Coefficient_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_algebraic_segment_traits_2<
                                          Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2< Coefficient >
                                                        Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }
  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods

  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  ArrangementPainterOstream& operator<<( const Point_2& p );

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected:
  void setupFacade( );
  
};

} // namespace Qt
} // namespace CGAL

// #if CGAL_EXPLICIT_INSTANTIATION == 0

// #include "ArrangementPainterOstream.cpp"

// #endif

#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
