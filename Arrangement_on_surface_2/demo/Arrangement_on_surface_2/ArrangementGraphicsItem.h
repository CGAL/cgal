// Copyright (c) 2008, 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

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
#include <QApplication>
#include <QKeyEvent>
#include <QPainter>
//#include <QStyleOption>

#include "ArrangementPainterOstream.h"
#include "Utils.h"

#include <iostream>
#include <limits>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "ui_ArrangementDemoWindow.h"
#include "ArrangementDemoGraphicsView.h"

class QGraphicsScene;
class QApplication;

extern Ui::ArrangementDemoWindow* getCurrentDemoWindowUi();
extern ArrangementDemoGraphicsView* getCurrentView();

namespace CGAL {
namespace Qt {

class ArrangementGraphicsItemBase :
    public GraphicsItem, public QGraphicsSceneMixin
{
public:
  ArrangementGraphicsItemBase( ):
    bb( 0, 0, 0, 0 ),
    verticesPen( QPen( ::Qt::blue, 3. ) ),
    edgesPen( QPen( ::Qt::blue, 1. ) ),
    visible_edges( true ),
    visible_vertices( true ),
    scene( NULL ),
    backgroundColor( ::Qt::white ),
    scalingFactor(1.0)
  {
    this->verticesPen.setCosmetic( true );
    this->verticesPen.setCapStyle( ::Qt::SquareCap );
    this->edgesPen.setCosmetic( true );
  }

  const QPen& getVerticesPen( ) const
  {
    return this->verticesPen;
  }

  const QPen& getEdgesPen( ) const
  {
    return this->edgesPen;
  }

  void setVerticesPen( const QPen& pen )
  {
    this->verticesPen = pen;
  }

  void setEdgesPen( const QPen& pen )
  {
    this->edgesPen = pen;
  }

  bool visibleVertices( ) const
  {
    return this->visible_vertices;
  }

  void setVisibleVertices( const bool b )
  {
    this->visible_vertices = b;
    this->update( );
  }

  bool visibleEdges( ) const
  {
    return this->visible_edges;
  }

  void setVisibleEdges( const bool b )
  {
    this->visible_edges = b;
    this->update( );
  }

  void setBackgroundColor( QColor color )
  {
    this->backgroundColor = color;
  }

  void setScene( QGraphicsScene* scene_ )
  {
    this->QGraphicsSceneMixin::setScene( scene_ );
    this->scene = scene_;
  }

protected:
  QRectF viewportRect( ) const
  {
    QRectF res;
    if ( this->scene == NULL )
    {
      return res;
    }

    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
      return res;
    }
    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );
    QPointF p1 = viewport->mapToScene( 0, 0 );
    QPointF p2 = viewport->mapToScene(viewport->width(), viewport->height());

    double xmin = (std::min)(p1.x(), p2.x());
    double xmax = (std::max)(p1.x(), p2.x());
    double ymin = (std::min)(p1.y(), p2.y());
    double ymax = (std::max)(p1.y(), p2.y());

    res = QRectF( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );

    return res;
  }

  CGAL::Bbox_2 bb;

  QPen verticesPen;
  QPen edgesPen;
  bool visible_edges;
  bool visible_vertices;

  QGraphicsScene* scene;

  QColor backgroundColor;
  double scalingFactor;

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
  typedef typename Arrangement::Unbounded_face_iterator Unbounded_face_iterator;
  typedef typename Arrangement::Hole_iterator           Holes_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;

  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  //typedef typename Kernel::Segment_2 Segment_2;

  typedef ArrangementGraphicsItemBase                   Superclass;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Traits::Curve_2                      Curve_2;

public:
  /*! Constructor */
  ArrangementGraphicsItem( Arrangement* t_ );

  /*! Destructor (virtual) */
  ~ArrangementGraphicsItem() {}

public:
  void modelChanged( );
  QRectF boundingRect( ) const override;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget) override;

protected:

  template < typename TTraits >
  void paint( QPainter* painter, TTraits traits );

  template < typename CircularKernel >
  void paint( QPainter* painter,
              CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

  template < typename Coefficient_ >
  void paint( QPainter* painter,
              CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > traits );

  void updateBoundingBox();

  template < typename TTraits>
  void updateBoundingBox(TTraits traits );

  template < typename Kernel_>
  void updateBoundingBox(CGAL::Arr_linear_traits_2<Kernel_> traits);

  template < typename RatKernel, class AlgKernel, class NtTraits >
  void updateBoundingBox(CGAL::Arr_Bezier_curve_traits_2<
                         RatKernel,
                         AlgKernel,
                         NtTraits > traits);
  template < typename Coefficient_>
  void updateBoundingBox(CGAL::Arr_algebraic_segment_traits_2<Coefficient_>
                         traits);

  Arrangement* arr;
  ArrangementPainterOstream< Traits > painterostream;
  CGAL::Qt::Converter< Kernel > convert;
  std::map< Curve_iterator, CGAL::Bbox_2 > curveBboxMap;

  void paintFaces( QPainter* painter )
  {
    painter->setPen(::Qt::transparent);
    typename Traits::Left_side_category category;
    this->paintFaces( painter, category );
    painter->setPen(edgesPen);
  }

  void paintFaces( QPainter* painter, CGAL::Arr_oblivious_side_tag )
  {
    // Prepare all faces for painting
    int Face_iterator_cnt = 0;
    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      fi->set_visited( false );
      Face_iterator_cnt++;
    }

    Face_iterator_cnt = 0;

    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      Face_iterator_cnt++;
      Face_handle f_handle = fi;
      this->paintFace( f_handle, painter );
    }
  }

  void paintFaces( QPainter* painter, CGAL::Arr_open_side_tag )
  {
    // Prepare all faces for painting
    int Face_iterator_cnt = 0;

    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      Face_iterator_cnt++;
      fi->set_visited( false );
    }

    Face_iterator_cnt = 0;

    for( Face_iterator fi = this->arr->faces_begin( );
         fi != this->arr->faces_end( ); ++fi )
    {
      Face_iterator_cnt++;
      Face_handle f_handle = fi;
      this->paintFace( f_handle, painter );
    }
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

  void paintFace( Face_handle f, QPainter* painter );

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
  {
  }

  template < typename Kernel_ >
  void paintFace( Face_handle f, QPainter* painter,
                  CGAL::Arr_segment_traits_2< Kernel_ > );

  template < typename Kernel_ >
  void paintFace( Face_handle f, QPainter* painter,
                  CGAL::Arr_polyline_traits_2< Kernel_ > );


#ifdef CGAL_USE_CORE
  template < typename RatKernel, typename AlgKernel, typename NtTraits >
  void
  paintFace( Face_handle f, QPainter* painter,
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
        {
          continue;
        }

        Halfedge_handle he = cc;
        X_monotone_curve_2 c = he->curve();
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double(he->source()->point().x()),
          sy = CGAL::to_double(he->source()->point().y()),
          tx = CGAL::to_double(he->target()->point().x()),
          ty = CGAL::to_double(he->target()->point().y());

        QPointF coord_source(sx, sy);
        QPointF coord_target(tx, ty);

        // Transform the point coordinates from general coordinate system to
        // Qt scene coordinate system
        QPoint coord_source_viewport = this->fromScene( coord_source );
        QPoint coord_target_viewport = this->fromScene( coord_target );

        if (c.orientation() == CGAL::COLLINEAR)
        {
          pts.push_back(coord_source );
        }
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

          // Draw the curve as pieces of small segments
          const int DRAW_FACTOR = 5;
          if (is_source_left)
          {
            for ( x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR )
            {
              //= COORD_SCALE)
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);

              // If curr_x > x_max or curr_x < x_min
              if (!(ker.compare_x_2_object()(curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object()(curr_p, c.right()) !=
                    CGAL::LARGER))
              {
                continue;
              }

              px = c.point_at_x (curr_p);
              curr_y = CGAL::to_double(px.y());
              QPointF curr( curr_x, curr_y );
              pts.push_back( curr );
            }// for
          }
          else
          {
            for ( x = x_max; x > x_min; x-=DRAW_FACTOR )
            {
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Arr_conic_point_2 curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object() (curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object() (curr_p, c.right()) !=
                    CGAL::LARGER))
              {
                continue;
              }

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
      painter->drawPolygon(pgn);
    }
    else
    {
      QRectF rect = this->viewportRect( );

      QColor color = this->backgroundColor;
      if ( f->color().isValid() )
      {
        color = f->color();
      }
      QBrush oldBrush = painter->brush( );
      painter->setBrush( color );
      painter->drawRect(rect);
      painter->setBrush( oldBrush );

    }
  }

#endif

#ifdef CGAL_USE_CORE
  template <typename RatKernel, typename AlgKernel, typename NtTraits >
  void
  paintFace( Face_handle f, QPainter* painter,
             CGAL::Arr_Bezier_curve_traits_2 <RatKernel, AlgKernel, NtTraits > ) {
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
        {
          continue;
        }

        Halfedge_handle he = cc;
        X_monotone_curve_2 c = he->curve();
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double(he->source()->point().x()),
          sy = CGAL::to_double(he->source()->point().y()),
          tx = CGAL::to_double(he->target()->point().x()),
          ty = CGAL::to_double(he->target()->point().y());

        QPointF coord_source(sx, sy);
        QPointF coord_target(tx, ty);

        // Transform the point coordinates from general coordinate system to
        // Qt scene coordinate system
        QPoint coord_source_viewport = this->fromScene( coord_source );
        QPoint coord_target_viewport = this->fromScene( coord_target );

        if (c.orientation() == CGAL::COLLINEAR)
        {
          pts.push_back(coord_source );
        }
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

          Bezier_point px;

          pts.push_back(coord_source );

          // Draw the curve as pieces of small segments
          const int DRAW_FACTOR = 5;
          if (is_source_left)
          {
            for ( x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR )
            {
              //= COORD_SCALE)
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Bezier_point curr_p(curr_x, 0);

              // If curr_x > x_max or curr_x < x_min
              if (!(ker.compare_x_2_object()(curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object()(curr_p, c.right()) !=
                    CGAL::LARGER))
              {
                continue;
              }

              px = c.point_at_x (curr_p);
              curr_y = CGAL::to_double(px.y());
              QPointF curr( curr_x, curr_y );
              pts.push_back( curr );
            }// for
          }
          else
          {
            for ( x = x_max; x > x_min; x-=DRAW_FACTOR )
            {
              curr_x = this->toScene( x );
              Alg_kernel   ker;
              Bezier_point curr_p(curr_x, 0);
              if (!(ker.compare_x_2_object() (curr_p, c.left()) !=
                    CGAL::SMALLER &&
                    ker.compare_x_2_object() (curr_p, c.right()) !=
                    CGAL::LARGER))
              {
                continue;
              }

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
      if ( f->color().isValid() )
      {
        color = f->color();
      }
      QBrush oldBrush = painter->brush( );
      painter->setBrush( color );
      painter->drawRect(rect);
      painter->setBrush( oldBrush );

    }
  }

  template < typename CircularKernel >
  void paintFace(Face_handle f, QPainter* painter,
                 CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */);

  template < typename Coefficient_ >
  void paintFace(Face_handle f, QPainter* painter,
                 CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
 {
    typedef Coefficient_ Coefficient;
    typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
    typedef typename Traits::X_monotone_curve_2         X_monotone_curve_2;

    Traits traits;

    if (f->is_unbounded()) return;

    // Only bounded face is drawn
    QVector< QPointF > pts; // holds the points of the polygon
    typedef typename Traits::CKvA_2                     CKvA_2;
    typedef std::pair< double, double >                 Coord_2;
    boost::optional < Coord_2 > p1, p2;
    typedef std::vector< Coord_2 >                      Coord_vec_2;
    std::list<Coord_vec_2> points;

    typedef Curve_renderer_facade<CKvA_2>               Facade;
    QGraphicsView* view = this->scene->views( ).first( );
    int height = view->height();

    QRectF viewport = this->viewportRect( );
    CGAL::Bbox_2 bbox = this->convert( viewport ).bbox( );
    Facade::setup(bbox, view->width(), view->height());
    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc=f->outer_ccb();
    do {
      double src_x = CGAL::to_double(cc->source()->point().x());
      // double src_y = CGAL::to_double(cc->source()->point().y());
      double tgt_x = CGAL::to_double(cc->target()->point().x());
      // double tgt_y = CGAL::to_double(cc->target()->point().y());

      X_monotone_curve_2 curve = cc->curve();
      Facade::instance().draw( curve, points, &p1, &p2 );

      if(points.empty())
      {
        continue;
      }

      QVector< QPointF > face_curve_points;
      typename std::list<Coord_vec_2>::const_iterator lit = points.begin();
      while(lit != points.end())
      {
        const Coord_vec_2& vec = *lit;
        // QPointF first = view->mapToScene( QPoint(vec[0].first, height-vec[0].second) );

        // QPointF last = view->mapToScene( QPoint(vec[vec.size()-1].first, height-vec[vec.size()-1].second) );

        typename Coord_vec_2::const_iterator vit = vec.begin();

        double sceneRectWidth = this->scene->width();
        double sceneRectHeight = this->scene->height();

        bool isPushBack = true;
        while ( vit != vec.end() )
        {
          QPoint coord( vit->first + sceneRectWidth/2, height - vit->second -sceneRectHeight/2);
          QPointF qpt = view->mapToScene( coord );
          if ( src_x < tgt_x )
          {
            face_curve_points.push_back( qpt );
            isPushBack = true;
          }
          else if ( src_x > tgt_x )
          {
            face_curve_points.push_front( qpt );
            isPushBack = false;
          }
          else // src_x == tgt_x
          {
            if (isPushBack)
            {
              face_curve_points.push_back( qpt );
            }
            else
            {
              face_curve_points.push_front( qpt );
            }
          }

          vit++;
          // std::cout << qpt.x() << "\t" << qpt.y() << std::endl;
        }
        lit++;
      }
      pts += face_curve_points;
      points.clear();
      face_curve_points.clear();
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
    QPen oldPen = painter->pen();
    painter->setBrush( color );

    painter->drawPolygon( pgn );

    painter->setBrush( oldBrush );
    painter->setPen( oldPen );

 }

  template < typename Kernel_ >
  void paintFace(Face_handle f, QPainter* painter,
                 CGAL::Arr_linear_traits_2< Kernel_ > /* traits */)
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
      // QRectF rect = this->viewportRect( );
      // QColor color = this->backgroundColor;
      // painter->fillRect( rect, color );
#if 0
      QRectF rect = this->viewportRect( );
      std::cout<<rect.left()<<'\t';
      std::cout<<rect.right()<<'\t';
      std::cout<<rect.top()<<'\t';
      std::cout<<rect.bottom()<<'\n';

      QColor color = this->backgroundColor;
      if ( f->color().isValid() )
      {
        color = f->color();
      }
      QBrush oldBrush = painter->brush( );
      painter->setBrush( color );
      painter->drawRect(rect);
      painter->setBrush( oldBrush );
#endif

    }
  }
#endif
#if 0
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
#endif

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

}; // class ArrangementGraphicsItem


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
