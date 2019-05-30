// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>


#include "ArrangementGraphicsItem.h"

namespace CGAL {
namespace Qt {

#if 0

QRectF ArrangementGraphicsItemBase::getViewportRect( ) const
{
  QRectF clipRect;
  if ( this->scene == NULL || this->scene->views( ).size( ) == 0 )
  {
    return clipRect;
  }

  QGraphicsView* view = this->scene->views( ).first( );
  QPointF p1 = view->mapToScene( 0, 0 );
  QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
  clipRect = QRectF( p1, p2 );

  return clipRect;
}

#endif

template < typename Arr_, class ArrTraits >
ArrangementGraphicsItem< Arr_, ArrTraits >::
ArrangementGraphicsItem( Arrangement* arr_ ):
  arr( arr_ ),
  painterostream( 0 )
{
  if ( this->arr->number_of_vertices( ) == 0 )
  {
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
  std::cout<<"In paint ArrTraits"<<std::endl;

  this->paintFaces( painter );

  painter->setPen( this->verticesPen );

  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
  this->painterostream.setScene( this->scene );

  QRectF rect = this->boundingRect( );
  std::cout<<"Curve boundingRect rect\n";
  std::cout<<"left, right, bottom, top:\n";
  std::cout<<rect.left()<<", "<<rect.right()<<", "<<rect.bottom()<<", "<<rect.top()<<std::endl;

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

    Bbox_2 bbox = curve.bbox();
    std::cout<<"Curve bounding box\n";
    std::cout<<"xmin, xmax, ymin, ymax:\n";
    std::cout<<bbox.xmin()<<", "<<bbox.xmax()<<", "<<bbox.ymin()<<", "<<bbox.ymax()<<std::endl;
    this->painterostream << curve;
  }
}

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
  std::cout<<"In paint Arr_circular_arc_traits_2"<<std::endl;
  this->paintFaces( painter );

  typedef Kernel_point_2 Non_arc_point_2;
  typedef typename Traits::Point_2 Arc_point_2;

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


template < typename Arr_, typename ArrTraits >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
{
  std::cout<<"In paint Arr_algebraic_segment_traits_2\n";
  painter->setPen( this->verticesPen );
  QRectF clipRect = this->boundingRect( );

  std::cout<<"left, right, bottom, top:\n";
  std::cout<<clipRect.left()<<", "<<clipRect.right()<<", "<<clipRect.bottom()<<", "<<clipRect.top()<<std::endl;

  if ( std::isinf(clipRect.left( )) ||
       std::isinf(clipRect.right( )) ||
       std::isinf(clipRect.top( )) ||
       std::isinf(clipRect.bottom( )) )
  {
    std::cout<<"In If with infinite bound\n";
    clipRect = this->viewportRect( );
  }

  std::cout<<"left, right, bottom, top:\n";
  std::cout<<clipRect.left()<<", "<<clipRect.right()<<", "<<clipRect.bottom()<<", "<<clipRect.top()<<std::endl;

  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, clipRect );
  this->painterostream.setScene( this->scene );

  std::cout<<"After initializing painterostream\n";

  for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( ); ++it )
  {
    Point_2 p = it->point( );
    //std::pair< double, double > approx = p.to_double( );
    //Kernel_point_2 pt( approx.first, approx.second );
    //this->painterostream << pt;
    this->painterostream << p;
  }

  std::cout<<"After done with vertices\n";

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->painterostream << curve;
  }
}

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
  std::cout<<"In updateBoundingBox Arr_linear_traits_2\n";
  this->prepareGeometryChange( );
  QRectF clipRect = this->viewportRect( );
  this->convert = Converter<Kernel>( clipRect );

  std::cout<<"left, right, bottom, top:\n";
  std::cout<<clipRect.left()<<", "<<clipRect.right()<<", "<<clipRect.bottom()<<", "<<clipRect.top()<<std::endl;

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

template < typename Arr_, typename ArrTraits >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
updateBoundingBox(CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits)
{
  std::cout<<"In updateBoundingBox Arr_algebraic_segment_traits_2\n";
  this->prepareGeometryChange( );
  if ( this->arr->number_of_vertices( ) == 0 )
  {
    this->bb = Bbox_2( 0, 0, 0, 0 );
    this->bb_initialized = false;
    std::cout<<"Leaving updateBoundingBox no vertex\n";
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
#if 0
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
#endif

  int curve_cnt = 0;

  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->bb = this->bb + curve.bbox( );
    std::cout<<"In updateBoundingBox for"<<std::endl;
    std::cout<<curve.bbox( ).xmin()<<"\t";
    std::cout<<curve.bbox( ).xmax()<<"\t";
    std::cout<<curve.bbox( ).ymin()<<"\t";
    std::cout<<curve.bbox( ).ymax()<<"\t"<<std::endl;

    curve_cnt++;

  }

  std::cout<<"curve_cnt\t"<<curve_cnt<<std::endl;
  std::cout<<"Leaving updateBoundingBox at the end\n";
}

template < typename Arr_, typename ArrTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
{
  std::cout<<"In ArrangementGraphicsItem modelChanged"<<std::endl;
  if ( this->arr->is_empty( ) )
  {
    this->hide( );
  }
  else
  {
    this->show( );
  }

  std::cout<<"In ArrangementGraphicsItem modelChanged after if"<<std::endl;
  this->updateBoundingBox( );
  this->update( );
}

template < typename Arr_, typename ArrTraits >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter )
{
  if ( f->visited( ) )
  {
    return;
  }

  std::cout<<"In paintFace f->visited( ) == false"<<std::endl;

  int holes = 0;
  int inner_faces = 0;
  Holes_iterator hit; // holes iterator
  this->paintFace( f, painter, Traits( ) );
  f->set_visited( true );

  for ( hit = f->holes_begin(); hit != f->holes_end(); ++hit )
  {
    // Traverse in clockwise order
    Ccb_halfedge_circulator cc = *hit;
    do {
      Halfedge_handle he = cc;
      Halfedge_handle he2 = he->twin();
      Face_handle inner_face = he2->face();
      if ( this->antenna( he ) )
      {
        continue;
      }

      // move on to next hole
      if ( ! inner_face->visited( ) )
      {
        inner_faces++;
      }

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

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
bool
ArrangementGraphicsItem< Arr_, ArrTraits >::
isProperOrientation( Ccb_halfedge_circulator cc )
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

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
bool 
ArrangementGraphicsItem< Arr_, ArrTraits >::
pathTouchingSource( const QPainterPath& path, X_monotone_curve_2 c )
{
  QPointF a = path.currentPosition( );
  QPointF b( to_double(c.source().x()), to_double(c.source().y()) );
  // QPointF d( to_double(c.target().x()), to_double(c.target().y()) );
  bool res = (QLineF( a, b ).length() < 1e-2);

  return res;
}

template < typename Arr_, typename ArrTraits >
template < typename Kernel_ >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter,
                CGAL::Arr_segment_traits_2< Kernel_ > )
{
  std::cout<<"In paintFace Arr_segment_traits_2"<<std::endl;

  if (!f->is_unbounded())  // f is not the unbounded face
  {
    std::cout<<"In paintFace Arr_segment_traits_2 bounded"<<std::endl;
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
    std::cout<<"In paintFace Arr_segment_traits_2 unbounded"<<std::endl;
    QRectF rect = this->viewportRect( );
    QColor color = this->backgroundColor;
    painter->fillRect( rect, color );
  }
}

template < typename Arr_, typename ArrTraits >
template < typename Kernel_ >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter,
                CGAL::Arr_polyline_traits_2< Kernel_ > )
{
  std::cout<<"In paintFace Arr_polyline_traits_2"<<std::endl;

  if (!f->is_unbounded())  // f is not the unbounded face
  {
    std::cout<<"In paintFace Arr_polyline_traits_2 bounded"<<std::endl;
    typedef typename CGAL::Arr_polyline_traits_2<Kernel_> Arr_poly_traits;
    typedef typename Arr_poly_traits::Compare_endpoints_xy_2 Comp_end_pts_2;
    typedef typename Arr_poly_traits::Construct_min_vertex_2 Poly_const_min_v;
    typedef typename Arr_poly_traits::Construct_max_vertex_2 Poly_const_max_v;

    // Obtain a polyline traits class and construct the needed functors
    Arr_poly_traits poly_tr;
    Comp_end_pts_2 comp_end_pts = poly_tr.compare_endpoints_xy_2_object();
    Poly_const_min_v poly_const_min_v=poly_tr.construct_min_vertex_2_object();
    Poly_const_max_v poly_const_max_v=poly_tr.construct_max_vertex_2_object();

    // Construct needed functors from the segment traits
    typedef typename Arr_poly_traits::Subcurve_traits_2      Subcurve_traits;
    typedef typename Subcurve_traits::Construct_min_vertex_2 Seg_const_min_v;
    typedef typename Subcurve_traits::Construct_max_vertex_2 Seg_const_max_v;
    Seg_const_min_v construct_min_v = poly_tr.subcurve_traits_2()->
      construct_min_vertex_2_object();
    Seg_const_max_v construct_max_v = poly_tr.subcurve_traits_2()->
      construct_max_vertex_2_object();

    // Iterator of the segments of an x-monotone polyline
    typename X_monotone_curve_2::Subcurve_const_iterator seg_it;

    QVector< QPointF > pts; // holds the points of the polygon
    X_monotone_curve_2 cv;

    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do {
      // The drawing is actually identical to segment
      double x = CGAL::to_double(cc->source()->point().x());
      double y = CGAL::to_double(cc->source()->point().y());
      QPointF coord_source(x , y);
      pts.push_back(coord_source );

      // The code below is problematic
      // cv = cc->curve();

      // // Determine the direction of cv (left-to-right or right-to-left)
      // Comparison_result dir = comp_end_pts(cv);

      // for (seg_it = cv.subcurves_begin();
      //      seg_it != cv.subcurves_end() ; ++seg_it)
      //   {
      //     if (dir == SMALLER)
      //       {
      //         // cv is directed from left-to-right
      //         // Adding the left-min vertex of the current segment
      //         double x = CGAL::to_double((construct_min_v(*seg_it)).x());
      //         double y = CGAL::to_double((construct_min_v(*seg_it)).y());
      //         QPointF coord_source(x , y);
      //         pts.push_back(coord_source );
      //       }
      //     else
      //       {
      //         // cv is directed from right-to-left
      //         // Adding the right-max vertex of the current segment
      //         double x = CGAL::to_double((construct_max_v(*seg_it)).x());
      //         double y = CGAL::to_double((construct_max_v(*seg_it)).y());
      //         QPointF coord_source(x , y);
      //         pts.push_back(coord_source );
      //       }
      //   }

      // if (dir == SMALLER)
      //   {
      //     // Add the right-most point of cv
      //     double x = CGAL::to_double((poly_const_max_v(cv)).x());
      //     double y = CGAL::to_double((poly_const_max_v(cv)).y());
      //     QPointF coord_source(x , y);
      //     pts.push_back(coord_source );
      //   }
      // else
      //   {
      //     // Add the left-most point of cv
      //     double x = CGAL::to_double((poly_const_min_v(cv)).x());
      //     double y = CGAL::to_double((poly_const_min_v(cv)).y());
      //     QPointF coord_source(x , y);
      //     pts.push_back(coord_source );
      //   }
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
    std::cout<<"In paintFace Arr_polyline_traits_2 unbounded"<<std::endl;

    // Draw a infinite bounding box
    QRectF rect = this->viewportRect( );
    QColor color = this->backgroundColor;
    painter->fillRect( rect, color );
  }
}

template < typename Arr_, typename ArrTraits >
template < typename RatKernel, typename AlgKernel, typename NtTraits >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter,
                CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > )
{
  std::cout<<"In paintFace Arr_conic_traits_2"<<std::endl;

  if (! f->is_unbounded())  // f is not the unbounded face
  {
    std::cout<<"In paintFace Arr_conic_traits_2 bounded"<<std::endl;
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
    painter->drawPolygon( pgn );
    painter->setBrush( oldBrush );
  }
  else
  {
    std::cout<<"In paintFace Arr_conic_traits_2 unbounded"<<std::endl;
    QRectF rect = this->viewportRect( );
    QColor color = this->backgroundColor;
    painter->fillRect( rect, color );
  }
}

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
void 
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace(Face_handle f, QPainter* painter,
               CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
{
  std::cout<<"In paintFace Arr_circular_arc_traits_2"<<std::endl;

  if ( f->is_unbounded( ) )
  {

    std::cout<<"In paintFace Arr_circular_arc_traits_2 unbounded"<<std::endl;
    QRectF rect = this->viewportRect( );
    QColor color = this->backgroundColor;
    if ( f->color().isValid() )
    {
      color = f->color();
    }
    painter->fillRect( rect, color );
    return;
  }

  std::cout<<"In paintFace Arr_circular_arc_traits_2 bounded"<<std::endl;

  QPainterPath path;
  bool isFirstArc = true;
  Ccb_halfedge_circulator cc=f->outer_ccb();
  do
  {
    if ( this->antenna( cc ) )
    {
      continue;
    }

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

} // namespace Qt
} // namespace CGAL

