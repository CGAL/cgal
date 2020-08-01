#include "ArrangementPainterOstream.h"

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

namespace CGAL {
namespace Qt {

// Instantiation of Arr_segment_traits_2
template < typename Kernel_ >
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >&
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> > ::
operator<<( const X_monotone_curve_2& curve )
{
  const Point_2& p1 = curve.source( );
  const Point_2& p2 = curve.target( );
  Segment_2 seg( p1, p2 );

  // skip segments outside our view
  QRectF seg_bb = this->convert( seg.bbox( ) );
  if ( this->clippingRect.isValid( ) &&
       ! this->clippingRect.intersects( seg_bb )
       && (!seg.is_horizontal() && !seg.is_vertical()))
  {
    return *this;
  }

  this->painterOstream << seg;
  return *this;
}

// Instantiation of Arr_polyline_traits_2

template < typename SegmentTraits >
ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> >&
ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> >::
operator<<( const X_monotone_curve_2& curve )
{
  int cnt = 0;
  for (typename X_monotone_curve_2::Subcurve_const_iterator it =
         curve.subcurves_begin();
       it != curve.subcurves_end(); ++it)
  {
    cnt++;
    this->painterOstream << *it;
  }

  // TODO: implement polyline painting
#if 0
  const Point_2& p1 = curve.source( );
  const Point_2& p2 = curve.target( );
  Segment_2 seg( p1, p2 );
  this->painterOstream << seg;
#endif
  return *this;
}

// Instantiation of Arr_conic_traits_2
template <typename RatKernel, class AlgKernel, class NtTraits>
auto ArrangementPainterOstream<CGAL::Arr_conic_traits_2<
	RatKernel, AlgKernel, NtTraits>>::visibleParts(X_monotone_curve_2 curve)
	-> std::vector<X_monotone_curve_2>
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

  intersect_2( bottom, curve, std::back_inserter( intersections ) );
  intersect_2( left, curve, std::back_inserter( intersections ) );
  intersect_2( top, curve, std::back_inserter( intersections ) );
  intersect_2( right, curve, std::back_inserter( intersections ) );

  this->filterIntersectionPoints( intersections );

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
    typename Traits::Point_2 p1 = *pointListItr++;
    typename Traits::Point_2 p2 = *pointListItr++;
    X_monotone_curve_2 subcurve =
      construct_x_monotone_subcurve_2( curve, p1, p2 );
    clippings.push_back( subcurve );
  }

  return clippings;
}

template <typename RatKernel, class AlgKernel, class NtTraits>
void ArrangementPainterOstream<
	CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
	filterIntersectionPoints(std::vector<CGAL::Object>& res)
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

template <typename RatKernel, class AlgKernel, class NtTraits>
void ArrangementPainterOstream<
	CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
	printIntersectResult(const std::vector<CGAL::Object>& res)
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
template < typename RatKernel, class AlgKernel, class NtTraits >
ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >&
ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >::
operator<<( const X_monotone_curve_2& curve )
{
  CGAL::Bbox_2 bb = curve.bbox( );
  QRectF qbb = this->convert( bb );

  // quick cull
  if (this->clippingRect.isValid() && ! this->clippingRect.intersects(qbb)) {
    return *this;
  }

  // get number of segments
  QGraphicsView* view = this->scene->views().first();
  int xmin = view->mapFromScene(bb.xmin(), bb.ymin()).x();
  int xmax = view->mapFromScene(bb.xmax(), bb.ymin()).x();
  // can be negitive due to rotation trasnformation
  size_t n = static_cast<size_t>(std::abs(xmax - xmin));
  if (n == 0) { return *this; }

  auto paintCurve = [&](auto&& curve_) {
    std::vector<std::pair<double, double>> app_pts;
    app_pts.reserve(n + 1);
    curve_.polyline_approximation(n, std::back_inserter(app_pts));

    auto p_curr = app_pts.begin();
    auto end_pts = app_pts.end();
    auto p_next = p_curr + 1;
    int count = 0;
    do
    {
      QPointF p1(p_curr->first, p_curr->second);
      QPointF p2(p_next->first, p_next->second);
      this->qp->drawLine(p1, p2);
      p_curr++;
      p_next++;
      ++count;
    } while (p_next != end_pts);
  };

  if ( this->clippingRect.isValid( ) )
  {
    std::vector< X_monotone_curve_2 > visibleParts;
    if ( this->clippingRect.contains( qbb ) )
      visibleParts.push_back( curve );
    else
      visibleParts = this->visibleParts( curve );

    for (auto& visiblePart : visibleParts)
      paintCurve(visiblePart);
  }
  else
  { // draw the whole curve
    paintCurve(curve);
  }

  return *this;
}

// Instantiation of Arr_Bezier_traits_2
template <typename RatKernel, class AlgKernel, class NtTraits>
std::vector<std::pair<double, double>>
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits>>::getPoints(const X_monotone_curve_2& curve)
{
  std::pair<double, double> param_range = curve.parameter_range();
  auto&& supporting_curve = curve.supporting_curve();

  std::vector<std::pair<double, double>> sampled_points;
  // TODO: get adaptive number of samples
  size_t number_of_samples = 100 * (param_range.second - param_range.first);
  sampled_points.reserve(number_of_samples);

  supporting_curve.sample(
    param_range.first, param_range.second, number_of_samples,
    std::back_inserter(sampled_points));
  return sampled_points;
}

template <typename RatKernel, class AlgKernel, class NtTraits>
ArrangementPainterOstream<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>&
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits>>::operator<<(const X_monotone_curve_2& curve)
{
  auto sampled_points = this->getPoints(curve);
  if (sampled_points.empty()) return *this;

  QPainterPath painterPath;
  painterPath.moveTo(sampled_points[0].first, sampled_points[0].second);

  for (auto& p : sampled_points)
    painterPath.lineTo(p.first, p.second);

  this->qp->drawPath(painterPath);
  return *this;
}

// Instantiation of Arr_linear_traits_2

template < typename Kernel_ >
ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >&
ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >::
operator<<( const X_monotone_curve_2& curve )
{
  if ( curve.is_segment( ) )
  {
    Segment_2 seg = curve.segment( );

    // skip segments outside our view
    QRectF seg_bb = this->convert( seg.bbox( ) );
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( seg_bb )
         & (!seg.is_horizontal() && !seg.is_vertical()))
    {
      return *this;
    }

    this->painterOstream << seg;
  }
  else if ( curve.is_ray( ) )
  {
    Ray_2 ray = curve.ray( );
    QLineF qseg = this->convert( ray );
    if ( qseg.isNull( ) )
    { // it's out of view
      return *this;
    }
    Segment_2 seg = this->convert( qseg );
    this-> painterOstream << seg;
  }
  else // curve.is_line( )
  {
    Line_2 line = curve.line( );
    QLineF qseg = this->convert( line );
    if ( qseg.isNull( ) )
    { // it's out of view
      return *this;
    }
    Segment_2 seg = this->convert( qseg );
    this-> painterOstream << seg;
  }
  return *this;
}

// Instantiation of Arr_algebraic_segment_traits_2
template <typename Traits>
static bool lies_on_border(const ArrangementPainterOstream<Traits> *apo,
                           const QPointF &point)
{
  QGraphicsView* view = apo->getScene()->views().first();
  qreal width = view->width();
  qreal height = view->height();
  const float tol = 2;
  return std::abs(point.x() - width) < tol || point.x() < tol ||
         std::abs(point.y() - height) < tol || point.y() < tol;
}

template <typename Coefficient_>
void ArrangementPainterOstream<
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::remapFacadePainter()
{
  this->qp->setTransform(this->getPointsListMapping());
}

template <typename Coefficient_>
QTransform ArrangementPainterOstream<
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::getPointsListMapping()
{
  auto worldTransform = this->qp->transform();
  return this->getPointsListMapping(worldTransform, this->getView());
}

template <typename Coefficient_>
QTransform
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  getPointsListMapping(
    const QTransform& worldTransform, const QGraphicsView* view)
{
  QRectF viewport = ArrangementDemoGraphicsView::viewportRect(view);

  // (0, 0) ==> map(topLeft)
  QPointF dxdy = worldTransform.map(viewport.topLeft());
  // (view.width(), 0) ==> map(bottomRight)
  QPointF p1 = worldTransform.map(viewport.topRight());
  // (0, view.height()) ==> map(0, 0)
  QPointF p2 = worldTransform.map(viewport.bottomLeft());

  // x' = m11*x + m21*y + dx
  // y' = m22*y + m12*x + dy
  float dx = dxdy.x();
  float dy = dxdy.y();
  float m11 = (p1.x() - dx) / view->width();
  float m21 = (p2.x() - dx) / view->height();
  float m22 = (p2.y() - dy) / view->height();
  float m12 = (p1.y() - dy) / view->width();

  return QTransform{m11, m12, m21, m22, dx, dy};
}

template <typename Coefficient_>
auto ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<
  Coefficient_>>::getPointsList(const X_monotone_curve_2& curve)
  -> std::list<Coord_vec_2>
{
  typedef Curve_renderer_facade<CKvA_2> Facade;
  typedef std::pair<double, double> Coord_2;
  typedef std::vector<Coord_2> Coord_vec_2;

  std::list<Coord_vec_2> points;
  Facade::instance().draw(curve, points);
  return points;
}

template < typename Coefficient_ >
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >&
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >::
operator<<( const X_monotone_curve_2& curve )
{
  this->qp->save();
  this->remapFacadePainter();
  this->paintCurve(curve);
  this->qp->restore();
  return *this;
}

template <typename Coefficient_>
void ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<
  Coefficient_>>::paintCurve(const X_monotone_curve_2& curve)
{
  std::list<Coord_vec_2> points = this->getPointsList(curve);
  for (auto& vec : points)
  {
    auto vit = vec.begin();
    QPainterPath path;
    QPointF qpt(vit->first, vit->second);
    path.moveTo(qpt);

    for (auto& vit : vec)
    {
      QPointF qpt_new = QPointF(vit.first, vit.second);
      if (lies_on_border(this, qpt) && lies_on_border(this, qpt_new))
        path.moveTo(qpt_new);
      else
        path.lineTo(qpt_new);
      qpt = qpt_new;
    }
    this->qp->drawPath(path);
  }
}

template <typename Coefficient_>
void ArrangementPainterOstream<
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  setupFacade()
{
  QGraphicsView* view = this->getView();
  typedef Curve_renderer_facade<CKvA_2> Facade;
  QRectF viewport = ArrangementDemoGraphicsView::viewportRect(view);
  CGAL::Bbox_2 bbox = this->convert(viewport).bbox();
  Facade::setup(bbox, view->width(), view->height());
}

template class ArrangementPainterOstream<Seg_traits>;
template class ArrangementPainterOstream<Pol_traits>;
template class ArrangementPainterOstream<Conic_traits>;
template class ArrangementPainterOstream<Lin_traits>;
template class ArrangementPainterOstream<Alg_seg_traits>;
template class ArrangementPainterOstream<Bezier_traits>;

} // namespace Qt
} // namespace CGAL
