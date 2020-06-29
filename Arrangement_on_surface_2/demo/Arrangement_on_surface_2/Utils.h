// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H

#include <CGAL/iterator.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>

#include "ArrangementDemoGraphicsView.h"
#include "ArrangementTypes.h"
#include "GraphicsSceneMixin.h"

class QGraphicsScene;
class QGraphicsSceneMouseEvent;

BOOST_MPL_HAS_XXX_TRAIT_DEF( Approximate_2 )

template <typename Arr_, bool b = has_Approximate_2< Arr_ >::value>
struct Supports_landmarks
{
  typedef CGAL::Boolean_tag< b > Tag;
};

template <typename Arr_>
struct Supports_landmarks< Arr_, true >
{
  typedef CGAL::Tag_true Tag;
};

/**
   Support for new ArrTraits should specify types:

   * Kernel - a not-necessarily-exact kernel to represent the arrangement
   graphically. We'll use the Point_2 type provided by this kernel for
   computing distances
   * Point_2 - the point type used in the particular arrangement
   * CoordinateType - the coordinate type used by the point type
   */
template <typename ArrTraits>
class ArrTraitsAdaptor
{ };

template <typename Kernel_>
class ArrTraitsAdaptor< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Kernel_>
class ArrTraitsAdaptor< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename SegmentTraits>
class ArrTraitsAdaptor< CGAL::Arr_polyline_traits_2< SegmentTraits > >
{
public:
  typedef CGAL::Arr_polyline_traits_2< SegmentTraits > ArrTraits;
  typedef typename SegmentTraits::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits >
class ArrTraitsAdaptor< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel,
                                                  NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef AlgKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class ArrTraitsAdaptor<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>
    ArrTraits;
  typedef RatKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Coefficient_>
class ArrTraitsAdaptor< CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > >
{
public:
  typedef Coefficient_ Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        ArrTraits;
  typedef typename ArrTraits::Point_2                   Point_2; // CKvA_2
  typedef typename ArrTraits::Algebraic_real_1          CoordinateType;
  typedef CGAL::Cartesian< typename ArrTraits::Bound >  Kernel;
  //typedef typename ArrTraits::CKvA_2                  Kernel;
};

template <typename ArrTraits >
class Arr_compute_y_at_x_2 : public QGraphicsSceneMixin
{
public:
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  // typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( );

  CoordinateType
  operator()(const X_monotone_curve_2& curve, const CoordinateType& x);

  double approx(const X_monotone_curve_2& curve, const CoordinateType& x);

protected:
  template <typename TTraits>
  CoordinateType operator()(
    const X_monotone_curve_2& curve, const CoordinateType& x, TTraits traits_,
    CGAL::Arr_oblivious_side_tag);

  template <typename TTraits>
  CoordinateType operator()(
    const X_monotone_curve_2& curve, const CoordinateType& x, TTraits traits_,
    CGAL::Arr_open_side_tag);

protected:
  Traits traits;
  Intersect_2 intersectCurves;
};

template <typename Coefficient_>
class Arr_compute_y_at_x_2< CGAL::Arr_algebraic_segment_traits_2<
                              Coefficient_ > > : public QGraphicsSceneMixin
{
public:
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
  typedef typename Traits::Algebraic_real_1             CoordinateType;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  CoordinateType operator()(
    const X_monotone_curve_2& curve, const CoordinateType& x,
    Point_2* out = nullptr);

  double approx(const X_monotone_curve_2& curve, const CoordinateType& x);

protected:
  X_monotone_curve_2 makeVerticalLine(const CoordinateType& x);
  Traits traits;
};

template <typename RatKernel, class AlgKernel, class NtTraits>
struct Arr_compute_y_at_x_2<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>> :
    public QGraphicsSceneMixin
{
  typedef CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>
    Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Rational Rational;
  typedef typename Traits::Algebraic Algebraic;
  typedef typename Traits::Point_2 Point_2;

  Algebraic operator()(
    const X_monotone_curve_2& curve, const Rational& x, Point_2* out = nullptr);

  Algebraic get_t(const X_monotone_curve_2& curve, const Rational& x);
  double approx(const X_monotone_curve_2& curve, const Rational& x);
};

template <class ArrTraits>
class Compute_squared_distance_2_base : public QGraphicsSceneMixin
{
public:
  typedef CGAL::Cartesian< double > InexactKernel;

public:
  // ctors
  Compute_squared_distance_2_base( ) { }

public: // methods

  template < class T1, class T2 >
  double operator() ( const T1& t1, const T2& t2 )
  {
    return this->squaredDistance( t1, t2 );
  }

protected: // fields
  typename Kernel::Compute_squared_distance_2 squared_distance;
  InexactKernel::Compute_squared_distance_2 squaredDistance;
};

template <typename ArrTraits >
class Compute_squared_distance_2 :
  public Compute_squared_distance_2_base< ArrTraits >
{ };

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_segment_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_segment_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_linear_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_linear_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_polyline_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_polyline_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_polyline_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Curve_2::Subcurve_const_iterator Seg_const_it;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Compute_squared_distance_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                            AlgKernel,
                                                            NtTraits > > :
  public Compute_squared_distance_2_base< CGAL::Arr_conic_traits_2< RatKernel,
                                                                    AlgKernel,
                                                                    NtTraits > >
{
public:
  typedef AlgKernel                                     Kernel;
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef Compute_squared_distance_2_base< Traits >     Superclass;
  // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
  typedef typename Traits::Point_2                      Conic_point_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public: // methods
  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Compute_squared_distance_2<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>> :
    public Compute_squared_distance_2_base<
      CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
{
public:
  typedef RatKernel Kernel;
  typedef CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>
    Traits;
  typedef Compute_squared_distance_2_base<Traits> Superclass;
  // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
  typedef typename Traits::Point_2 Conic_point_2;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // methods
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;
};

template <typename Coefficient_ >
class Compute_squared_distance_2< CGAL::Arr_algebraic_segment_traits_2<
                                    Coefficient_ > > :
  public Compute_squared_distance_2_base< CGAL::Arr_algebraic_segment_traits_2<
                                            Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        Traits;
  typedef typename Traits::Bound                        FT;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel     Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef std::pair< double, double >                   Coord_2;
  typedef std::vector< Coord_2 >                        Coord_vec_2;
  typedef typename X_monotone_curve_2::Curve_analysis_2 Curve;
  typedef typename Curve::Polynomial_traits_2           Polynomial_traits_2;
  typedef typename Curve::Polynomial_2                  Polynomial_2;
  typedef typename Polynomial_traits_2::Multivariate_content
                                                        Multivariate_content;
  typedef typename Polynomial_traits_2::Substitute      Substitute;
  typedef typename Polynomial_traits_2::
	  Construct_innermost_coefficient_const_iterator_range
		  ConstructInnerCoeffIter;

public:
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;

};

template <typename ArrTraits>
class Construct_x_monotone_subcurve_2
{
public:
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType
                                                        CoordinateType;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  Construct_x_monotone_subcurve_2( );

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const boost::optional<Point_2>& pLeft,
                                  const boost::optional<Point_2>& pRight );

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2


/*
 * This specialization for conic traits makes use of X_monotone_curve_2::trim,
 * which is not necessarily available.
 */
template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Construct_x_monotone_subcurve_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>
                                                        ArrTraits;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType
                                                        CoordinateType;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const boost::optional<Point_2>& pLeft,
                                  const boost::optional<Point_2>& pRight );

}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Construct_x_monotone_subcurve_2<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>
                                                        ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename ArrTraits::Point_2                   Point_2;

  Construct_x_monotone_subcurve_2();

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator()(
    const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
    const boost::optional<Point_2>& pRight);

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

// FIXME: return Traits::Point_2 instead of Kernel::Point_2
template <typename ArrTraits >
class SnapStrategy : public QGraphicsSceneMixin
{
public:
  //typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;

  virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event ) = 0;

protected:
  SnapStrategy( QGraphicsScene* scene_ );
}; // class SnapStrategy

template <typename ArrTraits >
SnapStrategy< ArrTraits >::SnapStrategy( QGraphicsScene* scene_ )
{
  this->scene = scene_;
}

template < typename ArrTraits>
class SnapToGridStrategy : public SnapStrategy<ArrTraits>
{
public:
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef SnapStrategy< ArrTraits >                     Superclass;

  /*! Constructors */
  SnapToGridStrategy( ) :
    Superclass( NULL ),
    gridSize( 50 )
  { }

  SnapToGridStrategy( QGraphicsScene* scene ) :
    Superclass( scene ),
    gridSize( 50 )
  { }

  /*! Destructors (virtual) */
  ~SnapToGridStrategy() {}

  Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
  {
    return this->snapPoint( event, ArrTraits( ) );
  }

  template <typename TTraits>
  Point_2 snapPoint(QGraphicsSceneMouseEvent* event, TTraits /* traits */)
  {
    QPointF clickedPoint = event->scenePos( );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    { // fallback case; we usually shouldn't end up here
      Kernel_point_2 res = this->convert( event->scenePos( ) );
      return Point_2( CGAL::to_double(res.x( )), CGAL::to_double(res.y()) );
    }

    qreal d( this->gridSize / 2.0 );
    int left = int( viewportRect.left( ) ) -
      (int( viewportRect.left( ) ) % this->gridSize);
    int right = int( viewportRect.right( ) ) +
      (this->gridSize - int( viewportRect.right( ) ) % this->gridSize);
    int x = int(clickedPoint.x( ));
    int y = int(clickedPoint.y( ));
    for ( int i = left - this->gridSize; i <= right; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.x( ) && clickedPoint.x( ) <= i + d )
      {
        x = i;
        break;
      }
    }
    int top = int( viewportRect.top( ) ) -
      (int( viewportRect.top( ) ) % this->gridSize);
    int bottom = int( viewportRect.bottom( ) ) +
      (this->gridSize - int( viewportRect.bottom( ) ) % this->gridSize);
    for ( int i = top - this->gridSize; i <= bottom; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.y( ) && clickedPoint.y( ) <= i + d )
      {
        y = i;
        break;
      }
    }
    //return this->convert( QPointF( x, y ) );
    Point_2 res( x, y );
    return res;
  }

  void setGridSize( int size )
  {
    this->gridSize = size;
  }

protected:
  int gridSize;
  CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToGridStrategy

template <typename Arr_>
class SnapToArrangementVertexStrategy:
  public SnapStrategy< typename Arr_::Geometry_traits_2 >
{
public:
  typedef Arr_                                          Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef SnapStrategy< Traits >                        Superclass;
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Kernel::Compute_squared_distance_2
    Compute_squared_distance_2;
  typedef typename Kernel::FT                           FT; typedef typename Traits::Point_2                      Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  SnapToArrangementVertexStrategy( ):
    Superclass( NULL ),
    arrangement( NULL )
  { }

  SnapToArrangementVertexStrategy( Arrangement* arr, QGraphicsScene* scene_ ):
    Superclass( scene_ ),
    arrangement( arr )
  { }

  Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
  {
    Kernel_point_2 clickedPoint = this->convert( event->scenePos( ) );
    return this->snapPoint( clickedPoint, Traits( ) );
  }

  template <typename TTraits>
  Point_2 snapPoint(const Kernel_point_2& clickedPoint, TTraits /* traits */)
  {
    Point_2 initialPoint( CGAL::to_double(clickedPoint.x()),
                          CGAL::to_double(clickedPoint.y()) );
    Point_2 closestPoint( CGAL::to_double(clickedPoint.x()),
                          CGAL::to_double(clickedPoint.y()) );
    bool first = true;
    FT minDist( 0 );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    {
      return initialPoint;
    }

    FT maxDist( ( viewportRect.right( ) - viewportRect.left( ) ) / 4.0 );
    for ( Vertex_iterator vit = this->arrangement->vertices_begin( );
          vit != this->arrangement->vertices_end( ); ++vit )
    {
      Point_2 point = vit->point( );
      Kernel_point_2 thisPoint( CGAL::to_double(point.x()),
                                CGAL::to_double(point.y()) );
      FT dist = this->compute_squared_distance_2( clickedPoint, thisPoint );
      if ( first || ( dist < minDist ) )
      {
        first = false;
        minDist = dist;
        closestPoint = point;
      }
    }
    if ( ! first && minDist < maxDist )
    {
      return closestPoint;
    }
    else
    {
      return initialPoint;
    }
  }

  void setArrangement( Arrangement* arr )
  {
    this->arrangement = arr;
  }

protected:
  Arrangement* arrangement;
  Compute_squared_distance_2 compute_squared_distance_2;
  CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToArrangementVertexStrategy

/**
   Converts between Kernel points and Arrangement points.

   The conversion is not necessarily exact.
*/
template <typename ArrTraits>
class Arr_construct_point_2
{
  typedef typename ArrTraits::Point_2                            Point_2;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType CoordinateType;
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel         Kernel;
  typedef typename Kernel::Point_2                               Kernel_point_2;

public:
  Point_2 operator()( const Kernel_point_2& pt );

  template <typename T >
  Point_2 operator()( const T& x, const T& y );

protected:
  template <typename T, typename TTraits >
  Point_2 operator()(const T& x, const T& y, TTraits /* traits */);
};

class Find_nearest_edge_base : public QGraphicsSceneMixin
{
public:
  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge_base() {}
};

template <typename Arr_, typename ArrTraits = typename Arr_::Geometry_traits_2>
class Find_nearest_edge : public Find_nearest_edge_base
{
public: // typedefs
  typedef Arr_ Arrangement;
  //typedef typename Arrangement::Geometry_traits_2 ArrTraits;
  typedef Compute_squared_distance_2< ArrTraits > Point_curve_distance;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Point_curve_distance::FT             FT;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

public:
  /*! constructor */
  Find_nearest_edge( Arrangement* arr_ ) :
    Find_nearest_edge_base( ),
    arr( arr_ )
  { }

  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge() {}

public: // member methods
  Halfedge_const_handle operator()( const Point_2& queryPt );

  virtual void setScene( QGraphicsScene* scene_ )
  {
    this->pointCurveDistance.setScene( scene_ );
    Find_nearest_edge_base::setScene( scene_ );
  }

protected: // member methods
  Face_const_handle getFace( const CGAL::Object& obj );

protected: // member fields
  Arrangement* arr;
  Point_curve_distance pointCurveDistance;
  Arr_construct_point_2< ArrTraits > toArrPoint;

}; // class Find_nearest_edge

#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
