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

template < typename Arr_ >
template < typename TTraits >
void ArrangementGraphicsItem< Arr_ >::
paint(QPainter* painter, TTraits /* traits */)
{
  this->paintFaces( painter );

  painter->setPen( this->verticesPen );

  auto painterostream =
    ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
  painterostream.setScene( this->getScene() );

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );

    // Bbox_2 bbox = curve.bbox();
    painterostream << curve;
  }
}

template < typename Arr_ >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_ >::
paint(QPainter* painter,
      CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
{
  QRectF clipRect = this->boundingRect( );

  // paint the faces for the purpose of brushing
  this->paintFaces( painter );

  // paint the curve itself
  painter->setPen( this->verticesPen );
  auto painterostream =
    ArrangementPainterOstream< Traits >( painter, clipRect );
  painterostream.setScene( this->getScene() );

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    painterostream << curve;
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

template < typename Arr_ >
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
{
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient> Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;

  if (f->is_unbounded()) return;

  // Only bounded face is drawn
  QVector<QPointF> pts; // holds the points of the polygon
  typedef typename Traits::CKvA_2 CKvA_2;
  typedef std::pair<double, double> Coord_2;
  typedef std::vector<Coord_2> Coord_vec_2;
  std::list<Coord_vec_2> points;

  // TODO: Same as setupFacade. Refactor it
  typedef Curve_renderer_facade<CKvA_2> Facade;
  QGraphicsView* view = this->getView();
  QRectF viewport = this->viewportRect();
  CGAL::Qt::Converter<Kernel> convert;
  CGAL::Bbox_2 bbox = convert(viewport).bbox();
  Facade::setup(bbox, view->width(), view->height());

  // TODO: Same as remapFacadePainter. Refactor it
  painter->save();
  auto worldTransform = painter->transform();
  QPointF dxdy = worldTransform.map(viewport.topLeft());
  QPointF p1 = worldTransform.map(viewport.topRight());
  QPointF p2 = worldTransform.map(viewport.bottomLeft());
  float dx = dxdy.x();
  float dy = dxdy.y();
  float m11 = (p1.x() - dx) / view->width();
  float m21 = (p2.x() - dx) / view->height();
  float m22 = (p2.y() - dy) / view->height();
  float m12 = (p1.y() - dy) / view->width();
  painter->setTransform(QTransform{m11, m12, m21, m22, dx, dy});

  /* running with around the outer of the face and generate from it
   * polygon
   */
  Ccb_halfedge_circulator cc = f->outer_ccb();
  do
  {
    double src_x = CGAL::to_double(cc->source()->point().x());
    double tgt_x = CGAL::to_double(cc->target()->point().x());

    auto curve = cc->curve();
    Facade::instance().draw(curve, points);

    if (points.empty()) { continue; }

    QVector<QPointF> face_curve_points;
    for (auto& vec : points)
    {
      bool isPushBack = true;
      for (auto& vit : vec)
      {
        QPoint qpt(vit.first, vit.second);
        if (src_x < tgt_x)
        {
          face_curve_points.push_back(qpt);
          isPushBack = true;
        }
        else if (src_x > tgt_x)
        {
          face_curve_points.push_front(qpt);
          isPushBack = false;
        }
        else // src_x == tgt_x
        {
          if (isPushBack)
            face_curve_points.push_back(qpt);
          else
            face_curve_points.push_front(qpt);
        }
      }
    }
    pts += face_curve_points;
    points.clear();
    face_curve_points.clear();
    // created from the outer boundary of the face
  } while (++cc != f->outer_ccb());

  // make polygon from the outer ccb of the face 'f'
  QPolygonF pgn(pts);

  // FIXME: get the bg color
  QColor color = this->backgroundColor;
  if (f->color().isValid()) { color = f->color(); }

  painter->setBrush(color);
  painter->drawPolygon(pgn);

  painter->restore();
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFaces(QPainter* painter)
{
  painter->setPen(::Qt::transparent);
  typename Traits::Left_side_category category;
  this->paintFaces(painter, category);
  painter->setPen(edgesPen);
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFaces(
  QPainter* painter, CGAL::Arr_oblivious_side_tag)
{
  // Prepare all faces for painting
  int Face_iterator_cnt = 0;
  for (Face_iterator fi = this->arr->faces_begin();
       fi != this->arr->faces_end(); ++fi)
  {
    fi->set_visited(false);
    Face_iterator_cnt++;
  }

  Face_iterator_cnt = 0;

  for (Face_iterator fi = this->arr->faces_begin();
       fi != this->arr->faces_end(); ++fi)
  {
    Face_iterator_cnt++;
    Face_handle f_handle = fi;
    this->paintFace(f_handle, painter);
  }
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFaces(
  QPainter* painter, CGAL::Arr_open_side_tag)
{
  // Prepare all faces for painting
  int Face_iterator_cnt = 0;

  for (Face_iterator fi = this->arr->faces_begin();
       fi != this->arr->faces_end(); ++fi)
  {
    Face_iterator_cnt++;
    fi->set_visited(false);
  }

  Face_iterator_cnt = 0;

  for (Face_iterator fi = this->arr->faces_begin();
       fi != this->arr->faces_end(); ++fi)
  {
    Face_iterator_cnt++;
    Face_handle f_handle = fi;
    this->paintFace(f_handle, painter);
  }
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
      Face_handle nei = (Face_handle)he.twin()->face();
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
    painter->setBrush(color);
    painter->drawRect(rect);
    painter->setBrush(oldBrush);
  }
}

template <typename Arr_>
template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>)
{
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

template class ArrangementGraphicsItem<Seg_arr>;
template class ArrangementGraphicsItem<Pol_arr>;
template class ArrangementGraphicsItem<Conic_arr>;
template class ArrangementGraphicsItem<Lin_arr>;
template class ArrangementGraphicsItem<Alg_seg_arr>;
template class ArrangementGraphicsItem<Bezier_arr>;

} // namespace QT
} // namespace CGAL
