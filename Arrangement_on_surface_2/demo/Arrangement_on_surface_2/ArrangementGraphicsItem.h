// Copyright (c) 2008, 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
//#include <CGAL/apply_to_range.h>
// TODO: should be included in PainterOstream.h
//#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include <QGraphicsScene>
#include <QPainter>
//#include <QStyleOption>

#include "ArrangementPainterOstream.h"
#include "Utils.h"
#include <iostream>

class QGraphicsScene;

namespace CGAL {
namespace Qt {

class ArrangementGraphicsItemBase :
    public GraphicsItem, public QGraphicsSceneMixin
{
public:
  ArrangementGraphicsItemBase( );

  const QPen& getVerticesPen( ) const;
  const QPen& getEdgesPen( ) const;
  void setVerticesPen( const QPen& pen );
  void setEdgesPen( const QPen& pen );
  bool visibleVertices( ) const;
  void setVisibleVertices( const bool b );
  bool visibleEdges( ) const;
  void setVisibleEdges( const bool b );
  void setBackgroundColor( QColor color );
  // void setScene( QGraphicsScene* scene_ );

protected:
  // QRectF getViewportRect( ) const;

  CGAL::Bbox_2 bb;
  bool bb_initialized;

  QPen verticesPen;
  QPen edgesPen;
  bool visible_edges;
  bool visible_vertices;

  QGraphicsScene* scene;

  QColor backgroundColor;

}; // class ArrangementGraphicsItemBase

template <typename Arr_, typename ArrTraits = typename Arr_::Geometry_traits_2>
class ArrangementGraphicsItem : public ArrangementGraphicsItemBase
{
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement::Curve_iterator          Curve_iterator;
  typedef typename Arrangement::Edge_iterator           Edge_iterator;
  typedef typename Arrangement::Halfedge                Halfedge;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Face_iterator           Face_iterator;
  typedef typename Arrangement::Hole_iterator           Holes_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;
 
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  //typedef typename Kernel::Segment_2 Segment_2;

  typedef ArrangementGraphicsItemBase                   Superclass;
  typedef typename Kernel::Segment_2                    Segment_2;

public:
  /*! Constructor */
  ArrangementGraphicsItem( Arrangement* t_ );

  /*! Destructor (virtual) */
  ~ArrangementGraphicsItem() {}
  
public:
  void modelChanged( );
  QRectF boundingRect( ) const;
  virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget);

protected:
  template < typename TTraits >
  void paint( QPainter* painter, TTraits traits );

  template < typename Kernel_ >
  void paint( QPainter* painter, CGAL::Arr_linear_traits_2< Kernel_ > traits );

  template < typename CircularKernel >
  void paint( QPainter* painter,
              CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

  // template < typename Coefficient_ >
  // void paint( QPainter* painter,
  //             CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > traits );


  void paintFaces( QPainter* painter )
  {
    typename Traits::Left_side_category category;
    this->paintFaces( painter, category );
  }

  void paintFaces( QPainter* painter, CGAL::Arr_oblivious_side_tag )
  {
    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      fi->set_visited( false );
    }

    Face_handle unboundedFace = this->arr->unbounded_face( );
    this->paintFace( unboundedFace, painter );
  }

  void paintFaces( QPainter* painter, CGAL::Arr_open_side_tag )
  {
    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      // if ( fi->is_fictitious( ) )
      // {
      //     std::cout << "setting fictitious face not visited" << std::endl;
      // }
      fi->set_visited( false );
    }

    Face_handle fictitiousFace = this->arr->fictitious_face( );
    this->paintFace( fictitiousFace, painter );
  }

#if 0
  template < typename TTraits >
  void paintFaces( QPainter* painter, TTraits traits ) { }

  template < typename Kernel_ >
  void paintFaces(QPainter* painter, CGAL::Arr_segment_traits_2<Kernel_>) {}

  template < typename Kernel_ >
  void paintFaces(QPainter* painter, CGAL::Arr_polyline_traits_2<Kernel_>) {}

  template < typename RatKernel, class AlgKernel, class NtTraits >
  void paintFaces(QPainter* painter,
                  CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits >)
  {}

  template < typename CircularKernel >
  void paintFaces(QPainter* painter,
                  CGAL::Arr_circular_arc_traits_2<CircularKernel > traits) {}

  template < typename Coefficient_ >
  void paintFaces(QPainter* painter,
                  CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits) {}
#endif

  void paintFace( Face_handle f, QPainter* painter )
  {
    if (! f->visited( ) )
    {
      int holes = 0;
      int inner_faces = 0;
      Holes_iterator hit; // holes iterator
      this->paintFace( f, painter, Traits( ) );
      f->set_visited( true );
      for ( hit = f->holes_begin(); hit != f->holes_end(); ++hit )
      {
        Ccb_halfedge_circulator cc = *hit;
        do {
          Halfedge_handle he = cc;
          Halfedge_handle he2 = he->twin();
          Face_handle inner_face = he2->face();
          if ( this->antenna( he ) )
            continue;

          // move on to next hole
          if ( ! inner_face->visited( ) )
            inner_faces++;
          this->visit_ccb_faces( inner_face, painter );
        } while ( ++cc != *hit );
        holes++;
      }// for
      // if ( f->is_unbounded( ) )
      // {
      //   std::cout << "unbounded face has " << holes << " holes" << std::endl;
      //   std::cout << "unbounded face has " << inner_faces << " inner faces"
      //             << std::endl;
      // }
      // if ( f->is_fictitious( ) )
      // {
      //   std::cout << "fictitious face has " << holes << " holes"
      //             << std::endl;
      //   std::cout << "fictitious face has " << inner_faces << " inner faces"
      //             << std::endl;
      // }
    }
  }

  void visit_ccb_faces( Face_handle & fh, QPainter* painter )
  {
    this->paintFace( fh, painter );

    Ccb_halfedge_circulator cc=fh->outer_ccb();
    do {
      Halfedge he = *cc;
      if (! he.twin()->face()->visited())
      {
        Face_handle nei = (Face_handle) he.twin()->face();
        this->visit_ccb_faces( nei , painter );
      }
      //created from the outer boundary of the face
    } while (++cc != fh->outer_ccb());
  }

  /*! antenna - return true if the halfedge and its
   *  twin point to the same face.
   */
  bool antenna(Halfedge_handle h)
  {
    Halfedge_handle twin = h->twin();
    return (twin->face() == h->face());
  }

  template < typename Traits >
  void paintFace(Face_handle /* f */, QPainter* /* painter */,
                 Traits /* traits */)
  { }

  template < typename Kernel_ >
  void paintFace( Face_handle f, QPainter* painter,
                  CGAL::Arr_segment_traits_2< Kernel_ > )
  {
    if (!f->is_unbounded())  // f is not the unbounded face
    {
      QVector< QPointF > pts; // holds the points of the polygon

      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        double x = CGAL::to_double(cc->source()->point().x());
        double y = CGAL::to_double(cc->source()->point().y());
        QPointF coord_source(x , y);
        pts.push_back(coord_source );
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      QPolygonF pgn (pts);

      // FIXME: get the bg color
      QColor color = this->backgroundColor;
      if ( f->color().isValid() )
      {
        color = f->color();
      }
      QBrush oldBrush = painter->brush( );
      painter->setBrush( color );
      painter->drawPolygon( pgn );
      painter->setBrush( oldBrush );
    }
    else
    {
      QRectF rect = this->viewportRect( );
      QColor color = this->backgroundColor;
      painter->fillRect( rect, color );
    }
  }

  template < typename Kernel_ >
  void paintFace( Face_handle f, QPainter* painter,
                  CGAL::Arr_polyline_traits_2< Kernel_ > )
  {
    if (!f->is_unbounded())  // f is not the unbounded face
    {
      QVector< QPointF > pts; // holds the points of the polygon
      typename X_monotone_curve_2::const_iterator           pt_itr;
      typename X_monotone_curve_2::const_reverse_iterator   pt_rev_itr;
      X_monotone_curve_2 cv;

      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc = f->outer_ccb();
      do {
        cv = cc->curve();
        bool curve_has_same_direction =
          ( *(cc->curve().begin()) == cc->source()->point() );
        if ( curve_has_same_direction )
        {
          for( pt_itr = cv.begin() , ++pt_itr ; pt_itr != cv.end(); ++pt_itr)
          {
            double x = CGAL::to_double((*pt_itr).x());
            double y = CGAL::to_double((*pt_itr).y());
            QPointF coord_source(x , y);
            pts.push_back(coord_source );
          }
        }
        else
        {
          for (pt_rev_itr = cv.rbegin() , ++pt_rev_itr; pt_rev_itr != cv.rend();
               ++pt_rev_itr)
          {
            double x = CGAL::to_double((*pt_rev_itr).x());
            double y = CGAL::to_double((*pt_rev_itr).y());
            QPointF coord_source(x , y);
            pts.push_back(coord_source );
          }
        }
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      QPolygonF pgn( pts );

      // fill the face according to its color (stored at any of her
      // incidents curves)
      QBrush oldBrush = painter->brush( );
      QColor def_bg_color = this->backgroundColor;
      if (! f->color().isValid())
      {
        painter->setBrush( def_bg_color );
      }
      else
      {
        painter->setBrush( f->color( ) );
      }
      painter->drawPolygon( pgn );
      painter->setBrush( oldBrush );
    }
    else
    {
      QRectF rect = this->viewportRect( );
      QColor color = this->backgroundColor;
      painter->fillRect( rect, color );
    }
  }

  template < typename RatKernel, typename AlgKernel, typename NtTraits >
  void paintFace( Face_handle f, QPainter* painter,
                  CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > )
  {
    if (! f->is_unbounded())  // f is not the unbounded face
    {
      QVector< QPointF > pts; // holds the points of the polygon
      /* running with around the outer of the face and generate from it
       * polygon
       */
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do 
      {
        if (this->antenna(cc))
          continue;

        Halfedge_handle he = cc;
        X_monotone_curve_2 c = he->curve();
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double(he->source()->point().x()),
          sy = CGAL::to_double(he->source()->point().y()),
          tx = CGAL::to_double(he->target()->point().x()),
          ty = CGAL::to_double(he->target()->point().y());

        QPointF coord_source(sx, sy);
        QPointF coord_target(tx, ty);
        QPoint coord_source_viewport = this->fromScene( coord_source );
        QPoint coord_target_viewport = this->fromScene( coord_target );

        if (c.orientation() == CGAL::COLLINEAR)
          pts.push_back(coord_source );
        else
        {
          // If the curve is monotone, than its source and its target has the
          // extreme x co-ordinates on this curve.
          bool is_source_left = (sx < tx);
          int  x_min = is_source_left ?
            coord_source_viewport.x( ) : coord_target_viewport.x( );
          int  x_max = is_source_left ?
            coord_target_viewport.x( ) : coord_source_viewport.x( );
          double curr_x, curr_y;
          int x;

          Arr_conic_point_2 px;

          pts.push_back(coord_source );

          const int DRAW_FACTOR = 5;
          if (is_source_left) 
          {
            for (x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR) 
            {
              //= COORD_SCALE)
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object()(curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object()(curr_p, c.right()) !=
                    CGAL::LARGER))
                continue;
              px = c.point_at_x (curr_p);
              curr_y = CGAL::to_double(px.y());
              QPointF curr( curr_x, curr_y );
              pts.push_back( curr );
            }// for
          }
          else 
          {
            for (x = x_max; x > x_min; x-=DRAW_FACTOR) 
            {
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object() (curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object() (curr_p, c.right()) !=
                    CGAL::LARGER))
                continue;
              px = c.point_at_x (curr_p);
              curr_y = CGAL::to_double(px.y());
              QPointF curr( curr_x, curr_y );
              pts.push_back( curr );
            }// for
          }// else
          pts.push_back(coord_target );
        }
        //created from the outer boundary of the face
      } while (++cc != f->outer_ccb());

      // make polygon from the outer ccb of the face 'f'
      QPolygonF pgn( pts );

      // fill the face according to its color (stored at any of her
      // incidents curves)
      QBrush oldBrush = painter->brush( );
      QColor def_bg_color = this->backgroundColor;
      if (! f->color().isValid())
      {
        painter->setBrush( def_bg_color );
      }
      else
      {
        painter->setBrush( f->color( ) );
      }
      painter->drawPolygon( pgn );
      painter->setBrush( oldBrush );
    }
    else
    {
      QRectF rect = this->viewportRect( );
      QColor color = this->backgroundColor;
      painter->fillRect( rect, color );
    }
  }

  template < typename CircularKernel >
  void paintFace(Face_handle f, QPainter* painter,
                 CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
  {
    if ( f->is_unbounded( ) )
    {
      QRectF rect = this->viewportRect( );
      QColor color = this->backgroundColor;
      if ( f->color().isValid() )
        color = f->color();
      painter->fillRect( rect, color );
      return;
    }
    QPainterPath path;
    bool isFirstArc = true;
    Ccb_halfedge_circulator cc=f->outer_ccb();
    do
    {
      if ( this->antenna( cc ) )
        continue;
            
      if ( isFirstArc )
      {
        isFirstArc = false;
        X_monotone_curve_2 c = cc->curve( );

        QPointF source( to_double(c.source().x()), to_double(c.source().y()) );
        QPointF target( to_double(c.target().x()), to_double(c.target().y()) );
        if ( ! this->isProperOrientation( cc ) )
        {
          std::swap( source, target );
        }

        QPointF circleCenter( to_double(c.supporting_circle().center().x()),
                              to_double(c.supporting_circle().center().y()) );
        //this->drawDiagnosticArc( circleCenter, source, target, painter );

        std::swap( source, target );
        double asource = atan2( -(source - circleCenter).y(),
                                (source - circleCenter).x() );
        double atarget = atan2( -(target - circleCenter).y(),
                                (target - circleCenter).x() );
        double aspan = atarget - asource;
        std::swap( source, target );

        path.moveTo( source );
        path.arcTo( convert(c.supporting_circle().bbox()), asource * 180/M_PI,
                    aspan * 180/M_PI );
        path.lineTo( target );
      }
      else
      {
        X_monotone_curve_2 c = cc->curve( );
        QPointF source( to_double(c.source().x()), to_double(c.source().y()) );
        QPointF target( to_double(c.target().x()), to_double(c.target().y()) );
        if ( ! this->pathTouchingSource( path, c ) )
          std::swap( source, target );

        QPointF circleCenter( to_double(c.supporting_circle().center().x()),
                              to_double(c.supporting_circle().center().y()) );
        //this->drawDiagnosticArc( circleCenter, source, target, painter );

        std::swap( source, target );
        double asource = atan2( -(source - circleCenter).y(),
                                (source - circleCenter).x() );
        double atarget = atan2( -(target - circleCenter).y(),
                                (target - circleCenter).x() );
        double aspan = atarget - asource;
        std::swap( source, target );

        path.arcTo( convert(c.supporting_circle().bbox()), asource * 180/M_PI,
                    aspan * 180/M_PI );
        path.lineTo( target );
      }
    } while (++cc != f->outer_ccb());
        
    if ( f->color().isValid() )
    {
      QPen savePen = painter->pen();
      QBrush saveBrush = painter->brush();
      QPen pen = painter->pen();
      pen.setColor( f->color() );

      painter->setPen( pen );
      painter->setBrush( f->color() );
      painter->drawPath( path );
      painter->setPen( savePen );
      painter->setBrush( saveBrush );
    }
  }

  void drawDiagnosticArc( QPointF c, QPointF s, QPointF t, QPainter* qp )
  {
    QBrush saveBrush = qp->brush( );
    double r = QLineF( c, s ).length();
    QRectF bb( c.x() - r, c.y() - r, 2*r, 2*r );

    QPainterPath path;
    // because we flipped the y-axis in our view, we need to flip our angles
    // with respect to the y-axis. we can do this by flipping s and t, then
    // negating their y-values
    std::swap( s, t );
    double as = atan2( -(s - c).y(), (s - c).x() );
    double at = atan2( -(t - c).y(), (t - c).x() );
    double aspan = at - as;

    path.moveTo( c );
    path.lineTo( s );
    path.arcTo( bb, as * 180/M_PI, aspan * 180/M_PI );
    path.closeSubpath( );
    qp->drawEllipse( bb );

    qp->setBrush( ::Qt::red );
    qp->fillPath( path, qp->brush() );
    qp->setBrush( saveBrush );

  }

  /**
     Return false if the tip of the given curve doesn't align with either of the
     endpoints of the next curve.
  */
  bool isProperOrientation( Ccb_halfedge_circulator cc )
  {
    Ccb_halfedge_circulator ccnext = cc;
    Halfedge_handle he = cc;
    X_monotone_curve_2 thisCurve = he->curve( );
    ccnext++;
    while ( this->antenna( ccnext ) ) ccnext++;
    Halfedge_handle next_he = ccnext;
    X_monotone_curve_2 nextCurve = next_he->curve( );

    QPointF thisTarget( to_double(thisCurve.target().x()),
                        to_double(thisCurve.target().y()) );
    QPointF nextSource( to_double(nextCurve.source().x()),
                        to_double(nextCurve.source().y()) );
    QPointF nextTarget( to_double(nextCurve.target().x()),
                        to_double(nextCurve.target().y()) );
    double dist1 = QLineF( thisTarget, nextSource ).length();
    double dist2 = QLineF( thisTarget, nextTarget ).length();
    bool res = ( dist1 < 1e-2 || dist2 < 1e-2 );

    return res;
  }

  bool pathTouchingSource( const QPainterPath& path, X_monotone_curve_2 c )
  {
    QPointF a = path.currentPosition( );
    QPointF b( to_double(c.source().x()), to_double(c.source().y()) );
    QPointF d( to_double(c.target().x()), to_double(c.target().y()) );
    bool res = (QLineF( a, b ).length() < 1e-2);

    return res;
  }

#if 0
  template < typename CircularKernel >
  void paintFace(Face_handle f, QPainter* painter,
                 CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
  {
    // std::cout << "face begin" << std::endl;
    if (! f->is_unbounded())  // f is not the unbounded face
    {
      QBrush oldBrush = painter->brush( );
      QColor def_bg_color = this->backgroundColor;
      if (! f->color().isValid())
      {
        painter->setBrush( def_bg_color );
      }
      else
      {
        painter->setBrush( f->color( ) );
      }

      // approach: draw a closed path with arcs
      QPainterPath path;
      //path.setFillRule( ::Qt::WindingFill );

      bool isFirstArc = true;
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do 
      {
        if (this->antenna(cc))
          continue;

        Halfedge_handle he = cc;
        X_monotone_curve_2 c = he->curve(); // circular arc

        const typename CircularKernel::Circle_2& circ = c.supporting_circle();
        const typename CircularKernel::Point_2& center = circ.center();
        QPointF source( to_double(c.source().x()), to_double(c.source().y()) );
        QPointF target( to_double(c.target().x()), to_double(c.target().y()) );

        if ( isFirstArc )
        {
          isFirstArc = false;
          // make sure we are going the right way
          Ccb_halfedge_circulator ccnext = cc;
          ccnext++;
          while ( this->antenna( ccnext ) ) ccnext++;
          Halfedge_handle next_he = ccnext;
          X_monotone_curve_2 next_c = next_he->curve( );
          QPointF next_source( to_double(next_c.source().x()),
                               to_double(next_c.source().y()) );
          QPointF next_target( to_double(next_c.target().x()),
                               to_double(next_c.target().y()) );
          // std::cout << "next curve's points: " << std::endl
          //           << "  s: " << next_source.x( )
          //           << " " << next_source.y( ) << std::endl
          //           << "  t: " << next_target.x( )
          //           << " " << next_target.y( ) << std::endl;
          double dist1 = QLineF( target, next_source ).length();
          double dist2 = QLineF( target, next_target ).length();
          if ( dist1 > 1e-2 && dist2 > 1e-2 )
          {
            std::cout << "swapping first curve points" << std::endl;
            std::swap( source, target );
          }

          path.moveTo( source );
        }
        else
        {
          //bool source_is_left = (source.x() < target.x());
          double dist = QLineF( path.currentPosition(), source ).length( );
          // if ( dist > 1e-2 )
          // {
          //     std::cout << "swapping source and target " << dist
          //               << std::endl;
          //     std::swap( source, target );
          // }
        }
        // std::cout << "currentPosition: " << path.currentPosition().x()
        //           << ", " << path.currentPosition().y() << std::endl;
        std::stringstream ss;
        ss << "  source: " << to_double(source.x( )) << ", "
           << to_double(source.y()) << std::endl
           << "  target: " << to_double(target.x( )) << ", "
           << to_double(target.y()) << std::endl;

#if 0
        source = QPointF( to_double(c.source().x()), to_double(c.source().y()));
        target = QPointF( to_double(c.target().x()), to_double(c.target().y()));
        double asource = std::atan2( -to_double(source.y() - center.y()),
                                     to_double(source.x() - center.x())); 
        double atarget = std::atan2( -to_double(target.y() - center.y()),
                                     to_double(target.x() - center.x()));
#endif
        double asource = std::atan2( to_double(source.y() - center.y()),
                                     to_double(source.x() - center.x()) );
        double atarget = std::atan2( to_double(target.y() - center.y()),
                                     to_double(target.x() - center.x()) );

        ss << "  asource: " << (asource * 180/M_PI) << std::endl
           << "  atarget: " << (atarget * 180/M_PI) << std::endl;

        //std::swap(asource, atarget);
        double aspan = atarget - asource;
        ss << "  aspan: " << (aspan * 180/M_PI) << std::endl;

        if(aspan < 0.)
          aspan += 2 * CGAL_PI;

        const double coeff = 180*16/CGAL_PI;
        double circR = sqrt(to_double(circ.squared_radius()));

#if 0
        painter->drawEllipse( to_double( center.x() - circR ),
                              to_double( center.y() - circR ),
                              2*circR, 2*circR );
#endif
        QPointF curPos = path.currentPosition( );
        //path.lineTo( to_double(target.x( )), to_double(target.y()) );

#if 0
        ss << "drawing line from " << curPos.x( ) << ", " << curPos.y( )
           << " to "
           << to_double(target.x( )) << ", " << to_double(target.y());
#endif
        // std::cout << ss.str( ) << std::endl;
        path.arcTo( convert(circ.bbox()), asource * 180/CGAL_PI,
                    aspan *180/CGAL_PI );
#if 0
        qp->drawArc(convert(circ.bbox()), 
                    (int)(asource * coeff), 
                    (int)(aspan * coeff));
#endif
      } while (++cc != f->outer_ccb());

      // fill the face according to its color (stored at any of her
      // incidents curves)
#if 0
      painter->drawPolygon( pgn );
#endif
      path.closeSubpath( );
#if 0
      QPainterPathStroker stroker;
      QPainterPath fillablePath = stroker.createStroke( path );
#endif
      painter->fillPath( path, painter->brush() );
#if 0
      QPolygonF fillIt = path.toFillPolygon( );
      painter->drawPolygon( fillIt );
#endif
      //painter->drawPath( path );
      //painter->fillPath( path, painter->brush() );
      painter->setBrush( oldBrush );
    }
    else
    {
      QRectF rect = this->viewportRect( );
      QColor color = this->backgroundColor;
      painter->fillRect( rect, color );
    }
  }
#endif

  //void cacheCurveBoundingRects( );

  void updateBoundingBox();

  template < typename TTraits>
  void updateBoundingBox(TTraits traits );

  template < typename Kernel_>
  void updateBoundingBox(CGAL::Arr_linear_traits_2<Kernel_> traits);

  // template < typename Coefficient_>
  // void updateBoundingBox(CGAL::Arr_algebraic_segment_traits_2<Coefficient_>
  //                        traits);

  Arrangement* arr;
  ArrangementPainterOstream< Traits > painterostream;
  CGAL::Qt::Converter< Kernel > convert;
  std::map< Curve_iterator, CGAL::Bbox_2 > curveBboxMap;
}; // class ArrangementGraphicsItem

  template < typename Arr_, class ArrTraits >
  ArrangementGraphicsItem< Arr_, ArrTraits >::
  ArrangementGraphicsItem( Arrangement* arr_ ):
    arr( arr_ ),
    painterostream( 0 )
  {
    if ( this->arr->number_of_vertices( ) == 0 ) {
      this->hide( );
    }
    this->updateBoundingBox( );
    this->setZValue( 3 );
  }

  template < typename Arr_, typename ArrTraits >
  QRectF 
  ArrangementGraphicsItem< Arr_, ArrTraits >::
  boundingRect( ) const
  {
    QRectF rect = this->convert( this->bb );
    return rect;
  }

  template < typename Arr_, typename ArrTraits >
  void 
  ArrangementGraphicsItem< Arr_, ArrTraits >::
  paint(QPainter* painter,
        const QStyleOptionGraphicsItem* /* option */,
        QWidget*  /*widget*/)
  {
    this->paint( painter, ArrTraits( ) );
  }

  template < typename Arr_, typename ArrTraits >
  template < typename TTraits >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  paint(QPainter* painter, TTraits /* traits */)
  {
    this->paintFaces( painter );

    painter->setPen( this->verticesPen );
    this->painterostream =
      ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
    this->painterostream.setScene( this->scene );

    for ( Vertex_iterator it = this->arr->vertices_begin( );
          it != this->arr->vertices_end( ); ++it )
    {
      Point_2 p = it->point( );
      Kernel_point_2 pt( p.x( ), p.y( ) );
      this->painterostream << pt;
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( );
          it != this->arr->edges_end( ); ++it )
    {
      X_monotone_curve_2 curve = it->curve( );
      this->painterostream << curve;
    }
  }

  template < typename Arr_, typename ArrTraits >
  template < typename Kernel_ >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  paint(QPainter* painter, CGAL::Arr_linear_traits_2< Kernel_ > /* traits */)
  {
    this->updateBoundingBox( );
    this->paintFaces( painter );
    painter->setPen( this->verticesPen );

    this->painterostream =
      ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
    this->painterostream.setScene( this->scene );

    for ( Vertex_iterator it = this->arr->vertices_begin( );
          it != this->arr->vertices_end( ); ++it )
    {
      Point_2 pt = it->point( );
      this->painterostream << pt;
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( );
          it != this->arr->edges_end( ); ++it )
    {
      X_monotone_curve_2 curve = it->curve( );
      this->painterostream << curve;
    }
  }

  template < typename Arr_, typename ArrTraits >
  template < typename CircularKernel >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  paint(QPainter* painter,
        CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
  {
    typedef Kernel_point_2 Non_arc_point_2;
    typedef typename Traits::Point_2 Arc_point_2;

    this->paintFaces( painter );

    painter->setPen( this->verticesPen );
    this->painterostream =
      ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
    this->painterostream.setScene( this->scene );

    for ( Vertex_iterator it = this->arr->vertices_begin( );
          it != this->arr->vertices_end( ); ++it )
    {
      Arc_point_2 pt = it->point( );
      Non_arc_point_2 pt2(CGAL::to_double(pt.x( )), CGAL::to_double(pt.y()) );
      this->painterostream << pt2;
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( );
          it != this->arr->edges_end( ); ++it )
    {
      X_monotone_curve_2 curve = it->curve( );
      this->painterostream << curve;
    }
  }

#if 0
  template < typename Arr_, typename ArrTraits >
  template < typename Coefficient_ >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  paint(QPainter* painter,
        CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
  {
    painter->setPen( this->verticesPen );
    QRectF clipRect = this->boundingRect( );
    if ( std::isinf(clipRect.left( )) || 
         std::isinf(clipRect.right( )) || 
         std::isinf(clipRect.top( )) || 
         std::isinf(clipRect.bottom( )) )
    {
      clipRect = this->viewportRect( );
    }

    this->painterostream =
      ArrangementPainterOstream< Traits >( painter, clipRect );
    this->painterostream.setScene( this->scene );

    for ( Vertex_iterator it = this->arr->vertices_begin( );
          it != this->arr->vertices_end( ); ++it )
    {
      Point_2 p = it->point( );
      //std::pair< double, double > approx = p.to_double( );
      //Kernel_point_2 pt( approx.first, approx.second );
      //this->painterostream << pt;
      this->painterostream << p;
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( );
          it != this->arr->edges_end( ); ++it )
    {
      X_monotone_curve_2 curve = it->curve( );
      this->painterostream << curve;
    }
  }
#endif

  // We let the bounding box only grow, so that when vertices get removed
  // the maximal bbox gets refreshed in the GraphicsView
  template < typename Arr_, typename ArrTraits >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::updateBoundingBox( )
  {
    this->updateBoundingBox( ArrTraits( ) );
  }

  template < typename Arr_, typename ArrTraits >
  template < typename TTraits >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  updateBoundingBox(TTraits /* traits */)
  {
    this->prepareGeometryChange( );
    if ( this->arr->number_of_vertices( ) == 0 )
    {
      this->bb = Bbox_2( 0, 0, 0, 0 );
      this->bb_initialized = false;
      return;
    }
    else
    {
      this->bb = this->arr->vertices_begin( )->point( ).bbox( );
      this->bb_initialized = true;
    }

    for ( Curve_iterator it = this->arr->curves_begin( );
          it != this->arr->curves_end( );
          ++it )
    {
      if ( this->curveBboxMap.count( it ) == 0 )
      {
        this->curveBboxMap[ it ] = it->bbox( );
      }
      this->bb = this->bb + this->curveBboxMap[ it ];
    }
  }

  template < typename Arr_, typename ArrTraits >
  template < typename Kernel_ >
  void
  ArrangementGraphicsItem< Arr_, ArrTraits >::
  updateBoundingBox(CGAL::Arr_linear_traits_2< Kernel_ > /* traits */)
  {
    this->prepareGeometryChange( );
    QRectF clipRect = this->viewportRect( );
    this->convert = Converter<Kernel>( clipRect );

    if ( ! clipRect.isValid( ) /*|| this->arr->number_of_vertices( ) == 0*/ )
    {
      this->bb = Bbox_2( 0, 0, 0, 0 );
      this->bb_initialized = false;
      return;
    }
    else
    {
      this->bb = this->convert( clipRect ).bbox( );
      this->bb_initialized = true;
    }

    for ( Curve_iterator it = this->arr->curves_begin( );
          it != this->arr->curves_end( );
          ++it )
    {
      if ( it->is_segment( ) )
      {
        this->bb = this->bb + it->segment( ).bbox( );
      }
      else if ( it->is_ray( ) )
      {
        QLineF qclippedRay = this->convert( it->ray( ) );
        Segment_2 clippedRay = this->convert( qclippedRay );
        this->bb = this->bb + clippedRay.bbox( );
      }
      else // ( it->is_line( ) )
      {
        QLineF qclippedLine = this->convert( it->line( ) );
        Segment_2 clippedLine = this->convert( qclippedLine );
        this->bb = this->bb + clippedLine.bbox( );
      }
    }
  }

#if 0
  template < typename Arr_, typename ArrTraits >
  template < typename Coefficient_ >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::
  updateBoundingBox(CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits)
  {
    this->prepareGeometryChange( );
    if ( this->arr->number_of_vertices( ) == 0 )
    {
      this->bb = Bbox_2( 0, 0, 0, 0 );
      this->bb_initialized = false;
      return;
    }
    else
    {
      //std::pair< double, double > approx =
      //  this->arr->vertices_begin( )->point( ).to_double( );
      //this->bb = CGAL::Bbox_2( approx.first, approx.second, 
      //                         approx.first, approx.second );
      this->bb = CGAL::Bbox_2( 0, 0, 0, 0 );
      this->bb_initialized = true;
    }
    typename Traits::Make_x_monotone_2 make_x_monotone_2 =
      traits.make_x_monotone_2_object( );
    for ( Curve_iterator it = this->arr->curves_begin( );
          it != this->arr->curves_end( );
          ++it )
    {
      std::vector< CGAL::Object > cvs;
      make_x_monotone_2( *it, std::back_inserter( cvs ) );
      for ( unsigned int i = 0 ; i < cvs.size( ); ++i )
      {
        X_monotone_curve_2 cv;
        CGAL::assign( cv, cvs[ i ] );
        this->bb = this->bb + cv.bbox( );
      }
    }
  }
#endif

  template < typename Arr_, typename ArrTraits >
  void ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
  {
    if ( this->arr->is_empty( ) )
    {
      this->hide( );
    }
    else
    {
      this->show( );
    }
    this->updateBoundingBox( );
    this->update( );
  }

  /**
     Specialized methods:
     updateBoundingBox
  */
#if 0
  template < typename Arr_, typename Kernel_ >
  class ArrangementGraphicsItem< Arr_, CGAL::Arr_linear_traits_2< Kernel_ > > :
    public ArrangementGraphicsItemBase
  {
    typedef Arr_ Arrangement;
    typedef ArrangementGraphicsItemBase Superclass;
    typedef typename Arrangement::Geometry_traits_2     Traits;
    typedef typename Arrangement::Vertex_iterator       Vertex_iterator;
    typedef typename Arrangement::Curve_iterator        Curve_iterator;
    typedef typename Arrangement::Edge_iterator         Edge_iterator;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::X_monotone_curve_2         X_monotone_curve_2;
    typedef typename Kernel::Point_2                    Point_2;
    typedef typename Kernel::Segment_2                  Segment_2;

  public:
    ArrangementGraphicsItem( Arrangement* arr_ ):
      arr( arr_ ),
      painterostream( 0 )
    {
      if ( this->arr->number_of_vertices( ) == 0 ) {
        this->hide( );
      }
      this->updateBoundingBox( );
      this->setZValue( 3 );
    }

  public: // methods
#if 0
    void modelChanged( )
    {
      if ( this->arr->is_empty( ) )
      {
        this->hide( );
      }
      else
      {
        this->show( );
      }
      this->updateBoundingBox( );
      this->update( );
    }

    // @override QGraphicsItem::boundingRect
    QRectF boundingRect( ) const
    {
      QRectF rect = this->convert( this->bb );
      return rect;
    }
#endif

    // @override QGraphicsItem::paint
#if 0
    virtual void paint( QPainter* painter,
                        const QStyleOptionGraphicsItem* option,
                        QWidget* widget )
    {
      this->updateBoundingBox( );
      painter->setPen( this->verticesPen );
      this->painterostream =
        ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
      this->painterostream.setScene( this->scene );

      for ( Vertex_iterator it = this->arr->vertices_begin( );
            it != this->arr->vertices_end( ); ++it )
      {
        Point_2 pt = it->point( );
        this->painterostream << pt;
      }
      painter->setPen( this->edgesPen );
      for ( Edge_iterator it = this->arr->edges_begin( );
            it != this->arr->edges_end( ); ++it )
      {
        X_monotone_curve_2 curve = it->curve( );
        this->painterostream << curve;
      }
    }
#endif

  protected: // methods
#if 0
    void updateBoundingBox( )
    {
      this->prepareGeometryChange( );
      QRectF clipRect = this->viewportRect( );
      this->convert = Converter<Kernel>( clipRect );

      if ( ! clipRect.isValid( ) /*|| this->arr->number_of_vertices( ) == 0*/ )
      {
        this->bb = Bbox_2( 0, 0, 0, 0 );
        this->bb_initialized = false;
        return;
      }
      else
      {
        this->bb = this->convert( clipRect ).bbox( );
        this->bb_initialized = true;
      }

      for ( Curve_iterator it = this->arr->curves_begin( );
            it != this->arr->curves_end( );
            ++it )
      {
        if ( it->is_segment( ) )
        {
          this->bb = this->bb + it->segment( ).bbox( );
        }
        else if ( it->is_ray( ) )
        {
          QLineF qclippedRay = this->convert( it->ray( ) );
          Segment_2 clippedRay = this->convert( qclippedRay );
          this->bb = this->bb + clippedRay.bbox( );
        }
        else // ( it->is_line( ) )
        {
          QLineF qclippedLine = this->convert( it->line( ) );
          Segment_2 clippedLine = this->convert( qclippedLine );
          this->bb = this->bb + clippedLine.bbox( );
        }
      }
    }
#endif

  protected: // fields
    Arrangement* arr;
    ArrangementPainterOstream< Traits > painterostream;
    CGAL::Qt::Converter< Kernel > convert;
  }; // class ArrangementGraphicsItem
#endif

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
