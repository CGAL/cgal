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

#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include "GraphicsSceneMixin.h"
#include "Utils.h"

class QPainter;

namespace CGAL {
namespace Qt {

template < typename ArrTraits >
class ArrangementPainterOstreamBase : public QGraphicsSceneMixin
{
public:
  // typedefs
  typedef ArrTraits                                     Traits;
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
    clippingRect( QRectF( ) ) // null rectangle
  {
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

  void setScene( QGraphicsScene* scene_ ) override
  {
    QGraphicsSceneMixin::setScene(scene_);

    // set the clipping rectangle
    if ( scene_ )
    {
      this->clippingRect = this->viewportRect( );
      this->convert = Converter< Kernel >( this->clippingRect );
    }
  }

protected:
  // fields
  PainterOstream< Kernel > painterOstream;
  QPainter* qp;
  Converter< Kernel > convert;
  QRectF clippingRect;
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
  std::vector< X_monotone_curve_2 > visibleParts( X_monotone_curve_2 curve );

  // keep only the intersection points ie. throw out overlapping curve segments
  void filterIntersectionPoints( std::vector< CGAL::Object >& res );

  void printIntersectResult( const std::vector< CGAL::Object >& res );

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
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits >
                                                        Traits;
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
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::FT                           FT;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle )
    { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  std::vector<std::pair<double, double>>
  getPoints(const X_monotone_curve_2& curve);

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    return *this;
  }
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
  typedef ArrangementPainterOstreamBase<Traits>         Super;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef std::pair<double, double>                     Coord_2;
  typedef std::vector<Coord_2>                          Coord_vec_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
      Superclass(p, clippingRectangle)
  {
  }

  void setScene(QGraphicsScene* scene_) override
  {
    Super::setScene(scene_);
    this->setupFacade();
  }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods

  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve );

  template <typename EdgesIterator>
  void paintEdges(EdgesIterator first, EdgesIterator last)
  {
    this->qp->save();
    this->remapFacadePainter();
    for (auto it = first; it != last; ++it)
    {
      X_monotone_curve_2 curve = it->curve();
      this->paintCurve(curve);
    }
    this->qp->restore();
  }

  // Maybe move these functions to someplace else?
  std::list<Coord_vec_2> getPointsList(const X_monotone_curve_2& curve);
  QTransform getPointsListMapping();
  QTransform getPointsListMapping(
    const QTransform& worldTransform, const QGraphicsView* view);
  void setupFacade();

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected:
  void paintCurve(const X_monotone_curve_2& curve);
  void remapFacadePainter();
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
