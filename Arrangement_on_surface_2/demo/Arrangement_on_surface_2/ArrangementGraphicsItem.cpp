#include "ArrangementGraphicsItem.h"
#include "ArrangementTypes.h"

namespace CGAL {
namespace Qt {

template < typename Arr_, class ArrTraits >
ArrangementGraphicsItem< Arr_, ArrTraits >::
ArrangementGraphicsItem( Arrangement* arr_ ):
  arr( arr_ ),
  painterostream( 0 )
{
  this->updatePointsItem();
  this->updateBoundingBox( );
  this->setZValue( 3 );
}

template < typename Arr_, typename ArrTraits >
QRectF
ArrangementGraphicsItem< Arr_, ArrTraits >::
boundingRect( ) const
{
  double xmin = -std::numeric_limits<double>::max() / 4;
  double ymin = -std::numeric_limits<double>::max() / 4;
  double xmax = std::numeric_limits<double>::max() / 4;
  double ymax = std::numeric_limits<double>::max() / 4;
  if (this->bb.xmin() > xmin) xmin = this->bb.xmin();
  if (this->bb.ymin() > ymin) ymin = this->bb.ymin();
  if (this->bb.xmax() < xmax) xmax = this->bb.xmax();
  if (this->bb.ymax() < ymax) ymax = this->bb.ymax();
  if (xmin > xmax || ymin > ymax)
  {
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
  }
  QRectF rect = this->convert(Bbox_2{xmin, ymin, xmax, ymax});
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

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );

    // Bbox_2 bbox = curve.bbox();
    this->painterostream << curve;
  }
}

template < typename Arr_, typename ArrTraits >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
{
  QRectF clipRect = this->boundingRect( );

  // paint the faces for the purpose of brushing
  this->paintFaces( painter );

  // paint the curve itself
  painter->setPen( this->verticesPen );
  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, clipRect );
  this->painterostream.setScene( this->scene );

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->painterostream << curve;
  }
}

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

  this->bb = {};
  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
  {
    // can throws "CGAL::internal::Zero_resultant_exception"
    try {
      this->bb += it->curve().bbox();
    } catch(...) {}
  }
}

template <typename Arr_, typename ArrTraits>
void ArrangementGraphicsItem<Arr_, ArrTraits>::updatePointsItem()
{
  this->pointsGraphicsItem.clear();
  for (auto it = this->arr->vertices_begin(); it != this->arr->vertices_end();
       ++it)
  {
    Point_2 p = it->point();
    this->pointsGraphicsItem.insert(p);
  }
}

template < typename Arr_, typename ArrTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
{
  this->updatePointsItem();
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
template < typename Kernel_ >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter,
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
    if ( f->color().isValid() )
    {
      color = f->color();
    }
    QBrush oldBrush = painter->brush( );
//    QPen pen = painter->pen();
//    pen.setCosmetic(true);
//    painter->setPen(pen);
    painter->setBrush( color );
    painter->drawRect(rect);
    painter->fillRect(rect, color);
    painter->setBrush( oldBrush );

  }
}

template < typename Arr_, typename ArrTraits >
template < typename Kernel_ >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace( Face_handle f, QPainter* painter,
                CGAL::Arr_polyline_traits_2< Kernel_ > )
{
  if (!f->is_unbounded())  // f is not the unbounded face
  {
    typedef typename CGAL::Arr_polyline_traits_2<Kernel_> Arr_poly_traits;
    typedef typename Arr_poly_traits::Compare_endpoints_xy_2 Comp_end_pts_2;
    typedef typename Arr_poly_traits::Construct_min_vertex_2 Poly_const_min_v;
    typedef typename Arr_poly_traits::Construct_max_vertex_2 Poly_const_max_v;

    // Obtain a polyline traits class and construct the needed functors
    Arr_poly_traits poly_tr;
    // Comp_end_pts_2 comp_end_pts = poly_tr.compare_endpoints_xy_2_object();
    // Poly_const_min_v poly_const_min_v=poly_tr.construct_min_vertex_2_object();
    // Poly_const_max_v poly_const_max_v=poly_tr.construct_max_vertex_2_object();

    // Construct needed functors from the segment traits
    typedef typename Arr_poly_traits::Subcurve_traits_2      Subcurve_traits;
    typedef typename Subcurve_traits::Construct_min_vertex_2 Seg_const_min_v;
    typedef typename Subcurve_traits::Construct_max_vertex_2 Seg_const_max_v;
    // Seg_const_min_v construct_min_v = poly_tr.subcurve_traits_2()->
    //   construct_min_vertex_2_object();
    // Seg_const_max_v construct_max_v = poly_tr.subcurve_traits_2()->
    //   construct_max_vertex_2_object();

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
    QPen pen = painter->pen();
    pen.setCosmetic(true);
    painter->setPen(pen);

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
    QPen pen = painter->pen();
    pen.setCosmetic(true);
    painter->setPen(pen);
    painter->setBrush( color );
    painter->drawRect(rect);
    painter->setBrush( oldBrush );

  }
}

template class ArrangementGraphicsItem<Seg_arr>;
template class ArrangementGraphicsItem<Pol_arr>;
template class ArrangementGraphicsItem<Conic_arr>;
template class ArrangementGraphicsItem<Lin_arr>;
template class ArrangementGraphicsItem<Alg_seg_arr>;

} // namespace QT
} // namespace CGAL
