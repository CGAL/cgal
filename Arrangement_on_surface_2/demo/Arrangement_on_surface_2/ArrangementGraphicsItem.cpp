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
  this->updateBoundingBox( );
  this->setZValue( 3 );
}

template < typename Arr_, typename ArrTraits >
QRectF
ArrangementGraphicsItem< Arr_, ArrTraits >::
boundingRect( ) const
{
  QRectF rect = this->convert( this->bb );
  if (!rect.isValid() || std::isinf(rect.width()) || std::isinf(rect.height()))
    return this->viewportRect();
  else
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

  // QRectF rect = this->boundingRect( );

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

    // Bbox_2 bbox = curve.bbox();
    this->painterostream << curve;
  }
}

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
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
  QRectF clipRect = this->boundingRect( );

  // paint the faces for the purpose of brushing
  this->paintFaces( painter );

  // paint the curve itself
  painter->setPen( this->verticesPen );
  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, clipRect );
  this->painterostream.setScene( this->scene );


  for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( ); ++it )
  {
    Point_2 p = it->point( );
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
    return;
  }
  else
  {
    this->bb = this->arr->vertices_begin( )->point( ).bbox( );
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

template <typename Arr_, typename ArrTraits >
template < typename RatKernel, class AlgKernel, class NtTraits >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
updateBoundingBox(CGAL::Arr_Bezier_curve_traits_2<
                  RatKernel,
                  AlgKernel,
                  NtTraits > /* traits */) {
  this->prepareGeometryChange( );
  QRectF clipRect = this->viewportRect( );
  this->convert = Converter<Kernel>( clipRect );

  if ( ! clipRect.isValid( ) /*|| this->arr->number_of_vertices( ) == 0*/ )
  {
    this->bb = Bbox_2( 0, 0, 0, 0 );
    return;
  }
  else
  {
    this->bb = this->convert( clipRect ).bbox( );
  }

//  int curve_cnt = 0;
//  for ( Edge_iterator it = this->arr->edges_begin( );
//        it != this->arr->edges_end( ); ++it )
//  {
//    X_monotone_curve_2 curve = it->curve( );
//    this->bb = this->bb + curve.bbox( );
//    curve_cnt++;
//  }
#if 0
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
#endif
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
    return;
  }
  else
  {
    this->bb = this->convert( clipRect ).bbox( );
  }

  for ( Curve_iterator it = this->arr->curves_begin( );
        it != this->arr->curves_end( ); ++it )
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

template <typename Traits>
static const std::vector<typename Traits::X_monotone_curve_2> &get_xy_curves() {
  static std::vector<typename Traits::X_monotone_curve_2> xy_curves;
  if (xy_curves.empty()) {
    Traits traits{};
    typedef typename Traits::Polynomial_2 Polynomial_2;
    typedef typename Traits::Curve_2 Curve_2;
    typename Traits::Construct_curve_2 construct_curve =
        traits.construct_curve_2_object();
    typename Traits::Make_x_monotone_2 make_x_monotone =
        traits.make_x_monotone_2_object();

    Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
    Curve_2 x_cv = construct_curve(x);
    Curve_2 y_cv = construct_curve(y);

    std::vector<CGAL::Object> arcs;
    make_x_monotone(x_cv, std::back_inserter(arcs));
    make_x_monotone(y_cv, std::back_inserter(arcs));
    std::transform(arcs.begin(), arcs.end(), std::back_inserter(xy_curves),
                   [](auto &&arc_obj) {
                     typename Traits::X_monotone_curve_2 arc;
                     CGAL::assign(arc, arc_obj);
                     return arc;
                   });
  }
  return xy_curves;
}

static bool is_finite(const CGAL::Bbox_2 &box) {
  return !std::isinf(box.xmin()) && !std::isinf(box.xmax()) &&
         !std::isinf(box.ymin()) && !std::isinf(box.ymax());
}

static CGAL::Bbox_2 make_finite(const CGAL::Bbox_2 &box) {
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  if (!std::isinf(box.xmin()))
    xmin = box.xmin();
  if (!std::isinf(box.ymin()))
    ymin = box.ymin();
  if (!std::isinf(box.xmax()))
    xmax = box.xmax();
  if (!std::isinf(box.ymax()))
    ymax = box.ymax();
  return CGAL::Bbox_2{xmin, ymin, xmax, ymax};
}

template <typename Arr_, typename ArrTraits>
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_, ArrTraits>::updateBoundingBox(
    CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits) {

  this->prepareGeometryChange();

  // we start with an empty box
  this->bb = {};

  // include vertices (intersections & finite x-monototne subcurve boundaries)
  for (auto it = arr->vertices_begin(); it != arr->vertices_end(); it++) {
    auto xy = it->point().to_double();
    this->bb += make_finite({xy.first, xy.second, xy.first, xy.second});
  }

  // include finite bounds of edges
  for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it) {
    // can throws "CGAL::internal::Zero_resultant_exception"
    try {
      this->bb += make_finite(it->curve().bbox());
    } catch(...) {}
  }

  // we didn't find any "interesting" point to include in the bounding box
  // we find intersections with x and y axis and add those instead
  if (!is_finite(this->bb)) {
    std::vector<CGAL::Object> intersections;
    for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it) {
      for (auto &arc : get_xy_curves<ArrTraits>()) {
        if (arc.is_vertical() != it->curve().is_vertical()) {
          it->curve().intersections(arc, std::back_inserter(intersections));
        }
      }
    }
    for (auto it = intersections.begin(); it != intersections.end(); it++) {
      std::pair<typename Traits::Point_2, unsigned int> point_multiplicity;
      CGAL::assign(point_multiplicity, *it);
      auto &point = point_multiplicity.first;
      if (point.location() == CGAL::ARR_INTERIOR) {
        auto xy = point.to_double();
        this->bb += make_finite({xy.first, xy.second, xy.first, xy.second});
      } else {
        // the curve is probably the x or y axis
        this->bb += {0, 0, 0, 0};
      }
    }
  }

  // this should happen only if the arrangement is empty
  if (!is_finite(this->bb))
    this->bb += {0, 0, 0, 0};

  // add margin to bounding box
  float x_margin;
  float y_margin;
  if (this->bb.xmin() == this->bb.xmax() ||
      this->bb.ymin() == this->bb.ymax()) {
    static constexpr float const_margin = 10;
    x_margin = const_margin;
    y_margin = const_margin;
  } else {
    static constexpr float prop_margin = 0.10;
    x_margin = (this->bb.xmax() - this->bb.xmin()) * prop_margin;
    y_margin = (this->bb.ymax() - this->bb.ymin()) * prop_margin;
  }
  this->bb = Bbox_2{this->bb.xmin() - x_margin, this->bb.ymin() - y_margin,
                    this->bb.xmax() + x_margin, this->bb.ymax() + y_margin};
}

template < typename Arr_, typename ArrTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
{
  this->updateBoundingBox( );
  QList< QGraphicsView* > views = this->scene->views( );
  if ( views.size( ) != 0 )
  {
    QGraphicsView* viewport = views.first( );
    viewport->setSceneRect(this->boundingRect());
  }
  this->update( );
}

template <>
void ArrangementGraphicsItem<Alg_seg_arr>::modelChanged()
{
  auto oldBoundingRect = this->boundingRect();
  this->updateBoundingBox( );
  QList< QGraphicsView* > views = this->scene->views( );
  if ( views.size( ) != 0 )
  {
    QGraphicsView* viewport = views.first( );
    auto newBoundingRect = this->boundingRect();
    if (oldBoundingRect != newBoundingRect) {
      // calling fitInView while moving the mouse causes the program to crash
      // this happens with the curve erasing tool
      viewport->setSceneRect(newBoundingRect);
      viewport->fitInView(this, ::Qt::KeepAspectRatio);
    }
  }
  this->update();
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

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
paintFace(Face_handle f, QPainter* painter,
               CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
{

  if ( f->is_unbounded( ) )
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

    return;
  }

  QVector< QPointF > pts;
  QPainterPath path;
  Ccb_halfedge_circulator cc=f->outer_ccb();
  int curve_cnt = 0;

  typedef CGAL::Arr_circular_arc_traits_2<CircularKernel> Traits;
  typedef CircularKernel                                Kernel;
  typedef typename Kernel::Root_of_2                    Root_of_2;

  Arr_compute_y_at_x_2< Traits > compute_y_at_x_2;

  do
  {
    if (this->antenna(cc))
    {
      continue;
    }
    curve_cnt++;

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

    if (false)
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

      pts.push_back(coord_source );

      // Draw the curve as pieces of small segments
      const int DRAW_FACTOR = 5;
      if (is_source_left)
      {
        for ( x = x_min + DRAW_FACTOR; x < x_max; x+=DRAW_FACTOR )
        {
          //= COORD_SCALE)
          curr_x = this->toScene( x );

          // If curr_x > x_max or curr_x < x_min
          if (curr_x < CGAL::to_double(c.left().x()) || curr_x > CGAL::to_double(c.right().x()))
          {
            continue;
          }

          curr_y = compute_y_at_x_2.approx(c, Root_of_2(curr_x));
          pts.push_back( QPointF( curr_x, curr_y ) );
        }// for
      }
      else
      {
        for ( x = x_max; x > x_min; x-=DRAW_FACTOR )
        {
          curr_x = this->toScene( x );
          if (curr_x < CGAL::to_double(c.left().x())
            || curr_x > CGAL::to_double(c.right().x()))
          {
            continue;
          }

          curr_y = compute_y_at_x_2.approx(c, Root_of_2(curr_x));
          pts.push_back( QPointF( curr_x, curr_y ) );
        }// for
      }// else
      pts.push_back(coord_target );
    }
    //created from the outer boundary of the face
  } while (++cc != f->outer_ccb());

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

template class ArrangementGraphicsItem<Seg_arr>;
template class ArrangementGraphicsItem<Pol_arr>;
template class ArrangementGraphicsItem<Conic_arr>;
template class ArrangementGraphicsItem<Lin_arr>;
template class ArrangementGraphicsItem<Arc_arr>;
template class ArrangementGraphicsItem<Alg_seg_arr>;

} // namespace QT
} // namespace CGAL
