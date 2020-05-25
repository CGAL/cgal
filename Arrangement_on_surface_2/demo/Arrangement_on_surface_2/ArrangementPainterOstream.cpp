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

template < typename Kernel_ >
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >&
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> > ::
operator<<( const Point_2& p )
{
  QPointF qpt = this->convert( p );
  // clip the point if possible
  if ( this->clippingRect.isValid( ) &&
       ! this->clippingRect.contains( qpt ) )
  {
    return *this;
  }

  QPen savePen = this->qp->pen( );
  this->qp->setBrush( QBrush( savePen.color( ) ) );
  double radius = savePen.width( ) / 2.0;
  radius /= this->scale;

  this->qp->drawEllipse( qpt, radius, radius );

  this->qp->setBrush( QBrush( ) );
  this->qp->setPen( savePen );
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

template < typename SegmentTraits >
ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> >&
ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> >::
operator<<( const Point_2& p )
{
  QPointF qpt = this->convert( p );
  QPen savePen = this->qp->pen( );
  this->qp->setBrush( QBrush( savePen.color( ) ) );
  double radius = savePen.width( ) / 2.0;
  radius /= this->scale;

  // Draw a circle as a blue dot
  this->qp->drawEllipse( qpt, radius, radius );

  this->qp->setBrush( QBrush( ) );
  this->qp->setPen( savePen );

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
  // std::cout << "In ArrangementPainterOstream& operator curve" << std::endl;

  CGAL::Bbox_2 bb = curve.bbox( );
  QRectF qbb = this->convert( bb );
  // quick cull

  // bool is_horizontal = (qbb.top() == qbb.bottom());
  // bool is_vertical = (qbb.left() == qbb.right());

  if (this->clippingRect.isValid() && ! this->clippingRect.intersects(qbb)) {
    //std::cout << "quick culled curve" << std::endl;
    return *this;
  }

#if 0
  std::cout << "bottom: ("
            << this->clippingRect.bottomLeft( ).x( )
            << " "
            << this->clippingRect.bottomLeft( ).y( )
            << " "
            << this->clippingRect.bottomRight( ).x( )
            << " "
            << this->clippingRect.bottomRight( ).y( )
            << ")"
            << std::endl;
#endif

  if ( this->clippingRect.isValid( ) )
  {
    std::vector< X_monotone_curve_2 > visibleParts;
    if ( this->clippingRect.contains( qbb ) )
    {
      visibleParts.push_back( curve );
    }
    else
    {
      visibleParts = this->visibleParts( curve );
    }

    for ( unsigned int i = 0; i < visibleParts.size( ); ++i )
    {
      X_monotone_curve_2 subcurve = visibleParts[ i ];
      int n;
      if ( this->scene == NULL )
      {
        n = 100; // TODO: get an adaptive approximation
      }
      else
      {
        QGraphicsView* view = this->scene->views( ).first( );
        int xmin, xmax;
        xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
        xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
        n = xmax - xmin;
      }
      if ( n == 0 )
      {
        return *this;
      }

      std::pair<double, double>* app_pts =
        new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        subcurve.polyline_approximation(n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );
#if 0
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
#endif
        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      }
      while ( p_next != end_pts );
    }
  }
  else
  { // draw the whole curve
    int n;
    if ( this->scene == NULL )
    {
      n = 100; // TODO: get an adaptive approximation
    }
    else
    {
      QGraphicsView* view = this->scene->views( ).first( );
      int xmin, xmax;
      xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
      xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
      n = xmax - xmin;
    }
    if ( n == 0 )
    {
      return *this;
    }

    std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
    std::pair<double, double>* end_pts =
      curve.polyline_approximation(n, app_pts);

    std::pair<double, double>* p_curr = app_pts;
    std::pair<double, double>* p_next = p_curr + 1;
    int count = 0;
    do
    {
      QPointF p1( p_curr->first, p_curr->second );
      QPointF p2( p_next->first, p_next->second );
#if 0
      Segment_2 seg( p1, p2 );
      this->painterOstream << seg;
#endif
      this->qp->drawLine( p1, p2 );
      p_curr++;
      p_next++;
      ++count;
    }
    while ( p_next != end_pts );
    //std::cout << count << " approximation points" << std::endl;
  }
  return *this;
}

  template < typename RatKernel, class AlgKernel, class NtTraits >
  ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                     NtTraits > >&
  ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                     NtTraits > >::
operator<<( const Point_2& p )
{
  // std::cout<< "In ArrangementPainterOstream& operator Point_2"<<std::endl;

  QPointF qpt = this->convert( p );
  QPen savePen = this->qp->pen( );
  this->qp->setBrush( QBrush( savePen.color( ) ) );
  double radius = savePen.width( ) / 2.0;
  radius /= this->scale;

  this->qp->drawEllipse( qpt, radius, radius );

  this->qp->setBrush( QBrush( ) );
  this->qp->setPen( savePen );
  return *this;
}

// Instantiation of Arr_Bezier_traits_2

template < typename RatKernel, class AlgKernel, class NtTraits >
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >&
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >::
operator<<( const X_monotone_curve_2& curve )
{
  // std::cout << "In ArrangementPainterOstream& operator curve" << std::endl;

  CGAL::Bbox_2 bb = curve.bbox( );
  QRectF qbb = this->convert( bb );
  // quick cull

  // bool is_horizontal = (qbb.top() == qbb.bottom());
  // bool is_vertical = (qbb.left() == qbb.right());

  if (this->clippingRect.isValid() && ! this->clippingRect.intersects(qbb)) {
    //std::cout << "quick culled curve" << std::endl;
    return *this;
  }

  if ( this->clippingRect.isValid( ) )
  {
    std::vector< X_monotone_curve_2 > visibleParts;
    if ( this->clippingRect.contains( qbb ) )
    {
      visibleParts.push_back( curve );
    }
    else
    {
      visibleParts = this->visibleParts( curve );
    }

    for ( unsigned int i = 0; i < visibleParts.size( ); ++i )
    {
      X_monotone_curve_2 subcurve = visibleParts[ i ];
      int n;
      if ( this->scene == NULL )
      {
        n = 100; // TODO: get an adaptive approximation
      }
      else
      {
        QGraphicsView* view = this->scene->views( ).first( );
        int xmin, xmax;
        xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
        xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
        n = xmax - xmin;
      }
      if ( n == 0 )
      {
        return *this;
      }

      std::pair<double, double>* app_pts =
        new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        subcurve.polyline_approximation(n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );
#if 0
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
#endif
        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      }
      while ( p_next != end_pts );
    }
  }
  else
  { // draw the whole curve
    int n;
    if ( this->scene == NULL )
    {
      n = 100; // TODO: get an adaptive approximation
    }
    else
    {
      QGraphicsView* view = this->scene->views( ).first( );
      int xmin, xmax;
      xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
      xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
      n = xmax - xmin;
    }
    if ( n == 0 )
    {
      return *this;
    }

    std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
    std::pair<double, double>* end_pts =
      curve.polyline_approximation(n, app_pts);

    std::pair<double, double>* p_curr = app_pts;
    std::pair<double, double>* p_next = p_curr + 1;
    int count = 0;
    do
    {
      QPointF p1( p_curr->first, p_curr->second );
      QPointF p2( p_next->first, p_next->second );
#if 0
      Segment_2 seg( p1, p2 );
      this->painterOstream << seg;
#endif
      this->qp->drawLine( p1, p2 );
      p_curr++;
      p_next++;
      ++count;
    }
    while ( p_next != end_pts );
    //std::cout << count << " approximation points" << std::endl;
  }
  return *this;
}

template < typename RatKernel, class AlgKernel, class NtTraits >
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >&
ArrangementPainterOstream<CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel,
                                                   NtTraits > >::
operator<<( const Point_2& p )
  {
    // std::cout<< "In ArrangementPainterOstream& operator Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
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

template < typename Kernel_ >
ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >&
ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >::
operator<<( const Point_2& p )
{
  QPointF qpt = this->convert( p );
  // clip the point if possible
  if ( this->clippingRect.isValid( ) &&
       ! this->clippingRect.contains( qpt ) )
  {
    return *this;
  }

  QPen savePen = this->qp->pen( );
  this->qp->setBrush( QBrush( savePen.color( ) ) );
  double radius = savePen.width( ) / 2.0;
  radius /= this->scale;

  this->qp->drawEllipse( qpt, radius, radius );

  this->qp->setBrush( QBrush( ) );
  this->qp->setPen( savePen );
  return *this;
}

// Instantiation of Arr_circular_arc_traits_2

template < typename CircularKernel >
ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<CircularKernel > >&
ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<CircularKernel > >::
operator<<( const X_monotone_curve_2& curve )
{
  this->painterOstream << curve;
  return *this;
}

template < typename CircularKernel >
ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<CircularKernel > >&
ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<CircularKernel > >::
operator<<( const Point_2& p )
{
  QPointF qpt = this->convert( p );
  // clip the point if possible
  if ( this->clippingRect.isValid( ) &&
       ! this->clippingRect.contains( qpt ) )
  {
    return *this;
  }

  QPen savePen = this->qp->pen( );
  this->qp->setBrush( QBrush( savePen.color( ) ) );
  double radius = savePen.width( ) / 2.0;
  radius /= this->scale;

  this->qp->drawEllipse( qpt, radius, radius );

  this->qp->setBrush( QBrush( ) );
  this->qp->setPen( savePen );
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
  const float tol = 3;
  return std::abs(point.x() - width) < tol || point.x() < tol ||
         std::abs(point.y() - height) < tol || point.y() < tol;
}

template < typename Coefficient_ >
void
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >::
remapFacadePainter()
{
  QGraphicsView* view = this->scene->views().first();
  this->qp->resetTransform();
  this->qp->translate(0, view->height());
  this->qp->scale(1.0, -1.0);
}

template < typename Coefficient_ >
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >&
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >::
operator<<( const X_monotone_curve_2& curve )
{
  typedef Curve_renderer_facade<CKvA_2> Facade;
  typedef std::pair<double, double> Coord_2;
  typedef std::vector<Coord_2> Coord_vec_2;
  this->setupFacade( );

  boost::optional<Coord_2> p1, p2;
  std::list<Coord_vec_2> points;
  Facade::instance().draw(curve, points, &p1, &p2);

  this->qp->save();
  this->remapFacadePainter();
  for (auto &vec : points) {
    auto vit = vec.begin();
    QPainterPath path;
    QPointF qpt(vit->first, vit->second);
    path.moveTo(qpt);

    while (++vit != vec.end()) {
      QPointF qpt_new = QPointF(vit->first, vit->second);
      if (lies_on_border(this, qpt) && lies_on_border(this, qpt_new) && false)
	    path.moveTo(qpt_new);
	  else
	    path.lineTo(qpt_new);
	  qpt = qpt_new;
    }
    this->qp->drawPath(path);
  }
  this->qp->restore();
  return *this;
}

template < typename Coefficient_ >
ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<Coefficient_ > >&
ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<Coefficient_ > >::
operator<<( const Point_2& p )
{
  typedef Curve_renderer_facade<CKvA_2> Facade;
  std::pair< double, double > coord;
  this->setupFacade( );

  this->qp->save();
  this->remapFacadePainter();
  // workaround for https://github.com/CGAL/cgal/issues/4745
  coord = p.to_double();
  if (Facade::instance().draw(p, coord))
  {
    QPointF qpt(coord.first, coord.second);
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    this->qp->drawEllipse( qpt, radius, radius );
    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
  }
  this->qp->restore();
  return *this;
}

template < typename Coefficient_ >
void
ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >::
setupFacade( )
{
  typedef Curve_renderer_facade<CKvA_2> Facade;
  QGraphicsView* view = this->scene->views( ).first( );
  QRectF viewport = this->viewportRect( );
  CGAL::Bbox_2 bbox = this->convert( viewport ).bbox( );
  Facade::setup(bbox, view->width(), view->height());
}


template class ArrangementPainterOstream<Seg_traits>;
template class ArrangementPainterOstream<Pol_traits>;
template class ArrangementPainterOstream<Conic_traits>;
template class ArrangementPainterOstream<Lin_traits>;
template class ArrangementPainterOstream<Arc_traits>;
template class ArrangementPainterOstream<Alg_seg_traits>;

} // namespace Qt
} // namespace CGAL
