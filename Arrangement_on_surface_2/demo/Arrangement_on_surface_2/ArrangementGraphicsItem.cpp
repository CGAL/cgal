#include "ArrangementGraphicsItem.h"
#include "ArrangementPainterOstream.h"
#include "ArrangementTypes.h"
#include <CGAL/Qt/Converter.h>
#include <QPainter>
#include <limits>

namespace CGAL {
namespace Qt {

ArrangementGraphicsItemBase::ArrangementGraphicsItemBase() :
    bb(0, 0, 0, 0), verticesPen(QPen(::Qt::blue, 3.)),
    edgesPen(QPen(::Qt::blue, 1.)), backgroundColor(::Qt::white)
{
  this->verticesPen.setCosmetic(true);
  this->verticesPen.setCapStyle(::Qt::SquareCap);
  this->edgesPen.setCosmetic(true);
  this->pointsGraphicsItem.setParentItem(this);
}

const QPen& ArrangementGraphicsItemBase::getVerticesPen() const
{
  return this->verticesPen;
}

const QPen& ArrangementGraphicsItemBase::getEdgesPen() const
{
  return this->edgesPen;
}

void ArrangementGraphicsItemBase::setVerticesPen(const QPen& pen)
{
  this->verticesPen = pen;
}

void ArrangementGraphicsItemBase::setEdgesPen(const QPen& pen)
{
  this->edgesPen = pen;
}

void ArrangementGraphicsItemBase::setBackgroundColor(QColor color)
{
  this->backgroundColor = color;
}

template < typename Arr_>
ArrangementGraphicsItem< Arr_>::
ArrangementGraphicsItem( Arrangement* arr_ ):
  arr( arr_ )
{
  this->updatePointsItem();
  this->updateBoundingBox( );
  this->setZValue( 3 );
}

template < typename Arr_ >
QRectF
ArrangementGraphicsItem< Arr_ >::
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
  return {QPointF{xmin, ymin}, QPointF{xmax, ymax}};
}

template < typename Arr_ >
void
ArrangementGraphicsItem< Arr_ >::
paint(QPainter* painter,
      const QStyleOptionGraphicsItem* /* option */,
      QWidget*  /*widget*/)
{
  this->paint( painter, Traits( ) );
}

template <typename Arr_>
template <typename TTraits>
void ArrangementGraphicsItem<Arr_>::paint(
  QPainter* painter, TTraits /* traits */)
{
  this->paintFaces(painter);

  painter->setPen(this->verticesPen);

  auto painterOstream =
    ArrangementPainterOstream<Traits>(painter, this->boundingRect());
  painterOstream.setScene(this->getScene());

  painter->setPen(this->edgesPen);
  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
  {
    X_monotone_curve_2 curve = it->curve();
    painterOstream << curve;
  }
}

template < typename Arr_ >
void ArrangementGraphicsItem< Arr_ >::updateBoundingBox( )
{
  this->updateBoundingBox( Traits( ) );
}

template < typename Arr_ >
template < typename TTraits >
void ArrangementGraphicsItem< Arr_ >::
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

template <typename Arr_>
template <typename RatK, typename AlgK, typename Nt, typename BoundingTratits>
void ArrangementGraphicsItem<Arr_>::updateBoundingBox(
  CGAL::Arr_Bezier_curve_traits_2<RatK, AlgK, Nt, BoundingTratits> /* traits */)
{
  this->prepareGeometryChange( );

  this->bb = {};
  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
  {
    // FIXME: There is no method to find bounding box of bezier x monotone curve
    // Look at Bezier_x_monotone_2.h Subcurve struct
    try {
      this->bb += it->curve().supporting_curve().bbox();
    } catch(...) {}
  }
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::updatePointsItem()
{
  this->pointsGraphicsItem.clear();
  for (auto it = this->arr->vertices_begin(); it != this->arr->vertices_end();
       ++it)
  {
    Point_2 p = it->point();
    this->pointsGraphicsItem.insert(p);
  }
}

// TODO: This is ugly. Clean it.
template <>
void ArrangementGraphicsItem<Bezier_arr>::updatePointsItem()
{
  this->pointsGraphicsItem.clear();
  for (auto it = this->arr->vertices_begin(); it != this->arr->vertices_end();
       ++it)
  {
    // Bezier_point can either be rational or exact
    // Calling Bezier_point::x() assumes that it is exact, which might not be
    // the case
    std::pair<double, double> p = it->point().approximate();
    QPointF qp {p.first, p.second};
    this->pointsGraphicsItem.insert(qp);
  }
}

template < typename Arr_ >
void ArrangementGraphicsItem< Arr_ >::modelChanged( )
{
  this->updatePointsItem();
  this->updateBoundingBox( );
  this->update( );
}

template < typename Arr_ >
void
ArrangementGraphicsItem< Arr_ >::
paintFace( Face_handle f, QPainter* painter )
{
  if (f->visited()) return;

  Holes_iterator hit; // holes iterator
  this->paintFace(f, painter, Traits());
  f->set_visited(true);

  for (hit = f->holes_begin(); hit != f->holes_end(); ++hit)
  {
    // Traverse in clockwise order
    Ccb_halfedge_circulator cc = *hit;
    do
    {
      Halfedge_handle he = cc;
      Halfedge_handle he2 = he->twin();
      Face_handle inner_face = he2->face();
      if (this->antenna(he)) continue;

      this->visit_ccb_faces(inner_face, painter);
    } while (++cc != *hit);
  }
}

template < typename Arr_ >
template < typename Kernel_ >
void
ArrangementGraphicsItem< Arr_ >::
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
    if ( f->color().isValid() ) color = f->color();
    QBrush oldBrush = painter->brush();
    painter->setBrush(color);
    painter->drawPolygon( pgn );
    painter->setBrush(oldBrush);
  }
  else
  {
    QRectF rect = this->viewportRect( );

    QColor color = this->backgroundColor;
    if ( f->color().isValid() )
    {
      color = f->color();
    }
    painter->fillRect(rect, color);
  }
}

template < typename Arr_ >
template < typename Kernel_ >
void
ArrangementGraphicsItem< Arr_ >::
paintFace( Face_handle f, QPainter* painter,
                CGAL::Arr_polyline_traits_2< Kernel_ > )
{
  if (!f->is_unbounded())  // f is not the unbounded face
  {
    // typedef typename CGAL::Arr_polyline_traits_2<Kernel_> Arr_poly_traits;
    // typedef typename Arr_poly_traits::Compare_endpoints_xy_2 Comp_end_pts_2;
    // typedef typename Arr_poly_traits::Construct_min_vertex_2 Poly_const_min_v;
    // typedef typename Arr_poly_traits::Construct_max_vertex_2 Poly_const_max_v;

    // Obtain a polyline traits class and construct the needed functors
    // Arr_poly_traits poly_tr;
    // Comp_end_pts_2 comp_end_pts = poly_tr.compare_endpoints_xy_2_object();
    // Poly_const_min_v poly_const_min_v=poly_tr.construct_min_vertex_2_object();
    // Poly_const_max_v poly_const_max_v=poly_tr.construct_max_vertex_2_object();

    // Construct needed functors from the segment traits
    // typedef typename Arr_poly_traits::Subcurve_traits_2      Subcurve_traits;
    // typedef typename Subcurve_traits::Construct_min_vertex_2 Seg_const_min_v;
    // typedef typename Subcurve_traits::Construct_max_vertex_2 Seg_const_max_v;
    // Seg_const_min_v construct_min_v = poly_tr.subcurve_traits_2()->
    //   construct_min_vertex_2_object();
    // Seg_const_max_v construct_max_v = poly_tr.subcurve_traits_2()->
    //   construct_max_vertex_2_object();

    // Iterator of the segments of an x-monotone polyline
    // typename X_monotone_curve_2::Subcurve_const_iterator seg_it;

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
    QBrush oldBrush = painter->brush();
    if (!f->color().isValid())
      painter->setBrush(this->backgroundColor);
    else
      painter->setBrush(f->color());

    painter->drawPolygon( pgn );
    painter->setBrush( oldBrush );
  }
  else
  {
    QRectF rect = this->viewportRect( );
    QColor color = this->backgroundColor;
    if (f->color().isValid()) color = f->color();
    painter->fillRect(rect, color);
  }
}

template < typename Arr_ >
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
{
  // FIXME: paint unbounded faces as well
  if (f->is_unbounded()) return;

  ArrangementPainterOstream<Traits> painterOstream{
    painter, this->boundingRect()};
  painterOstream.setScene(this->getScene());

  painter->save();
  painter->setTransform(painterOstream.getPointsListMapping());

  QVector<QPointF> pts;
  /* running with around the outer of the face and generate from it
   * polygon
   */
  Ccb_halfedge_circulator cc = f->outer_ccb();
  do
  {
    if (this->antenna(cc)) continue;

    auto points = painterOstream.getPointsList(cc->curve());
    if (points.empty()) continue;

    double src_x = CGAL::to_double(cc->source()->point().x());
    double tgt_x = CGAL::to_double(cc->target()->point().x());

    if (src_x > tgt_x)
    {
      for (auto& vec : points)
          for(auto& vit : vec)
            pts.push_back({vit.first, vit.second});
    }
    else
    {
      for (auto& vec : points)
        for (auto vit = vec.rbegin(); vit != vec.rend(); ++vit)
          pts.push_back({vit->first, vit->second});
    }
  } while (++cc != f->outer_ccb());

  // make polygon from the outer ccb of the face 'f'
  QPolygonF pgn(pts);
  QColor color = this->backgroundColor;
  if (f->color().isValid()) { color = f->color(); }
  painter->setBrush(color);
  painter->drawPolygon(pgn);

  painter->restore();
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFaces(QPainter* painter)
{
  QPen pen = painter->pen();
  painter->setPen(QColorConstants::Transparent);

  // Prepare all faces for painting
  for (auto fi = this->arr->faces_begin(); fi != this->arr->faces_end(); ++fi)
    fi->set_visited(false);

  for (auto fi = this->arr->faces_begin(); fi != this->arr->faces_end(); ++fi)
    this->paintFace(fi, painter);

  painter->setPen(pen);
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::visit_ccb_faces(
  Face_handle& fh, QPainter* painter)
{
  this->paintFace(fh, painter);

  Ccb_halfedge_circulator cc = fh->outer_ccb();
  do
  {
    Halfedge he = *cc;
    if (!he.twin()->face()->visited())
    {
      Face_handle nei = he.twin()->face();
      this->visit_ccb_faces(nei, painter);
    }
    // created from the outer boundary of the face
  } while (++cc != fh->outer_ccb());
}

template <typename Arr_>
bool ArrangementGraphicsItem<Arr_>::antenna(Halfedge_handle h)
{
  Halfedge_handle twin = h->twin();
  return (twin->face() == h->face());
}

template <typename Arr_>
template <typename ArrTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(Face_handle, QPainter*, ArrTraits)
{
}

template <typename Arr_>
template <typename RatKernel, typename AlgKernel, typename NtTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>)
{

  if (!f->is_unbounded()) // f is not the unbounded face
  {
    QVector<QPointF> pts; // holds the points of the polygon
    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      if (this->antenna(cc)) { continue; }

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
      QPoint coord_source_viewport = this->fromScene(coord_source);
      QPoint coord_target_viewport = this->fromScene(coord_target);

      if (c.orientation() == CGAL::COLLINEAR) { pts.push_back(coord_source); }
      else
      {
        // If the curve is monotone, than its source and its target has the
        // extreme x co-ordinates on this curve.
        bool is_source_left = (sx < tx);
        int x_min = is_source_left ? coord_source_viewport.x()
                                   : coord_target_viewport.x();
        int x_max = is_source_left ? coord_target_viewport.x()
                                   : coord_source_viewport.x();
        double curr_x, curr_y;
        int x;

        Arr_conic_point_2 px;

        pts.push_back(coord_source);

        // Draw the curve as pieces of small segments
        const int DRAW_FACTOR = 5;
        if (is_source_left)
        {
          for (x = x_min + DRAW_FACTOR; x < x_max; x += DRAW_FACTOR)
          {
            //= COORD_SCALE)
            curr_x = this->toScene(QPoint{x, 0}).x();
            Alg_kernel ker;
            Arr_conic_point_2 curr_p(curr_x, 0);

            // If curr_x > x_max or curr_x < x_min
            if (!(ker.compare_x_2_object()(curr_p, c.left()) != CGAL::SMALLER &&
                  ker.compare_x_2_object()(curr_p, c.right()) != CGAL::LARGER))
            { continue; }

            px = c.point_at_x(curr_p);
            curr_y = CGAL::to_double(px.y());
            QPointF curr(curr_x, curr_y);
            pts.push_back(curr);
          } // for
        }
        else
        {
          for (x = x_max; x > x_min; x -= DRAW_FACTOR)
          {
            curr_x = this->toScene(QPoint{x, 0}).x();
            Alg_kernel ker;
            Arr_conic_point_2 curr_p(curr_x, 0);
            if (!(ker.compare_x_2_object()(curr_p, c.left()) != CGAL::SMALLER &&
                  ker.compare_x_2_object()(curr_p, c.right()) != CGAL::LARGER))
            { continue; }

            px = c.point_at_x(curr_p);
            curr_y = CGAL::to_double(px.y());
            QPointF curr(curr_x, curr_y);
            pts.push_back(curr);
          } // for
        }   // else
        pts.push_back(coord_target);
      }
      // created from the outer boundary of the face
    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face 'f'
    QPolygonF pgn(pts);
    // fill the face according to its color (stored at any of her
    // incidents curves)
    QBrush oldBrush = painter->brush();
    QColor def_bg_color = this->backgroundColor;
    if (!f->color().isValid()) { painter->setBrush(def_bg_color); }
    else
    {
      painter->setBrush(f->color());
    }
    painter->drawPolygon(pgn);
  }
  else
  {
    QRectF rect = this->viewportRect();

    QColor color = this->backgroundColor;
    if (f->color().isValid()) { color = f->color(); }
    QBrush oldBrush = painter->brush();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits> traits)
{
  if (!f->is_unbounded())
  {
    ArrangementPainterOstream<Traits> painterOstream{
      painter, this->boundingRect()};
    painterOstream.setScene(this->getScene());

    QVector<QPointF> pts;
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      if (this->antenna(cc)) continue;

      const X_monotone_curve_2& curve = cc->curve();
      auto points = painterOstream.getPoints(curve);

      if (cc->source()->point().is_same(curve.source()))
        for (auto& vit : points)
          pts.push_back({vit.first, vit.second});
      else
        for (auto it = points.rbegin(); it != points.rend(); ++it)
          pts.push_back({it->first, it->second});

    } while(++cc != f->outer_ccb());

    QPolygonF pgn(pts);
    QColor color = this->backgroundColor;
    if (f->color().isValid()) { color = f->color(); }

    QBrush oldBrush = painter->brush();
    painter->setBrush(color);
    painter->drawPolygon(pgn);
    painter->setBrush(oldBrush);
  }
  else
  {
    QRectF rect = this->viewportRect();

    QColor color = this->backgroundColor;
    if (f->color().isValid()) color = f->color();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <typename Kernel_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_linear_traits_2<Kernel_> /* traits */)
{

  if (!f->is_unbounded()) // f is not the unbounded face
  {
    QVector<QPointF> pts; // holds the points of the polygon

    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      double x = CGAL::to_double(cc->source()->point().x());
      double y = CGAL::to_double(cc->source()->point().y());
      QPointF coord_source(x, y);
      pts.push_back(coord_source);
      // created from the outer boundary of the face
    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face 'f'
    QPolygonF pgn(pts);

    // FIXME: get the bg color
    QColor color = this->backgroundColor;
    if (f->color().isValid()) { color = f->color(); }
    QBrush oldBrush = painter->brush();
    painter->setBrush(color);
    painter->drawPolygon(pgn);
    painter->setBrush(oldBrush);
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
      if ( f->color().isValid() ) color = f->color();
      painter->fillRect(rect, color);
#endif
  }
}

template class ArrangementGraphicsItem<Seg_arr>;
template class ArrangementGraphicsItem<Pol_arr>;
template class ArrangementGraphicsItem<Conic_arr>;
template class ArrangementGraphicsItem<Lin_arr>;
template class ArrangementGraphicsItem<Alg_seg_arr>;
template class ArrangementGraphicsItem<Bezier_arr>;

} // namespace QT
} // namespace CGAL
