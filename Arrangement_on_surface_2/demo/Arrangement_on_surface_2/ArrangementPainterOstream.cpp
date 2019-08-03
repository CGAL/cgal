
// #ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_CPP
// #define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_CPP

#include "ArrangementPainterOstream.h"

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

template < typename Coefficient_ >
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >&
ArrangementPainterOstream<CGAL::Arr_algebraic_segment_traits_2<Coefficient_> >::
operator<<( const X_monotone_curve_2& curve )
{
  //std::cout << "paint curve stub (alg traits)" << std::endl;

  typedef Curve_renderer_facade<CKvA_2> Facade;
  typedef std::pair< double, double > Coord_2;
  typedef std::vector< Coord_2 > Coord_vec_2;
  this->setupFacade( );
  boost::optional < Coord_2 > p1, p2;
  std::list<Coord_vec_2> points;
  Facade::instance().draw( curve, points, &p1, &p2 );

  if (points.empty()) return *this;

  // QPainter *ppnt = this->qp;
  QGraphicsView* view = this->scene->views( ).first( );
  int height = view->height();
  // int width = view->width();

  typename std::list<Coord_vec_2>::const_iterator lit = points.begin();
  //ppnt->moveTo((*p1).first, height - (*p1).second);
  while(lit != points.end()) {
    const Coord_vec_2& vec = *lit;

    typename Coord_vec_2::const_iterator vit = vec.begin();
    //std::cerr << "(" << vit->first << "; " << vit->second << ")\n";
    //         if(lit == points.begin() &&*/ vit != vec.end()) {
    //             ppnt->lineTo(vit->first, height - vit->second);
    //             vit++;
    //         }
#if 0
    if(vit != vec.end())
      ppnt->moveTo(vit->first, height - vit->second);
    while(vit != vec.end()) {
      ppnt->lineTo(vit->first, height - vit->second);
      vit++;
      //std::cerr << "(" << vit->e0 << "; " << vit->e1 << "\n";
    }
    lit++;
#endif
    QPainterPath path;
    double sceneRectWidth = this->scene->width();
    double sceneRectHeight = this->scene->height();

    QPoint coord( vit->first + sceneRectWidth/2, height - vit->second -sceneRectHeight/2 );
    QPointF qpt = view->mapToScene( coord );

    QPoint last_coord( vec[vec.size()-1].first, height - vec[vec.size()-1].second );
    // QPointF last = view->mapToScene( last_coord );

    if ( vit != vec.end() )
    {
      path.moveTo( qpt );
    }
    while ( vit != vec.end() )
    {
      path.lineTo( qpt );
      vit++;
      coord = QPoint( vit->first + sceneRectWidth/2, height - vit->second -sceneRectHeight/2 );
      qpt = view->mapToScene( coord );
      // std::cout << qpt.x() << "\t" << qpt.y() << std::endl;
    }
    this->qp->drawPath( path );
    lit++;
  }
  //ppnt->lineTo((*p2).first, height - (*p2).second);
#if 0
  QPen old_pen = ppnt->pen();
  ppnt->setPen(QPen(Qt::NoPen)); // avoid drawing outlines
  // draw with the current brush attributes
  //std::cerr << "endpts1: (" << (*p1).first << "; " << (*p1).second << "\n";
  //std::cerr << "endpts2: (" << (*p2).first << "; " << (*p2).second << "\n";
  unsigned sz = CGAL_REND_PT_RADIUS;
  ppnt->drawEllipse((*p1).first - sz, height-(*p1).second - sz, sz*2, sz*2);
  ppnt->drawEllipse((*p2).first - sz, height-(*p2).second - sz, sz*2, sz*2);
  ppnt->setPen(old_pen);
#endif
  return *this;
}

template < typename Coefficient_ >
ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<Coefficient_ > >&
ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<Coefficient_ > >::
operator<<( const Point_2& p )
{
  typedef Curve_renderer_facade<CKvA_2> Facade;
  std::pair< double, double > coord;
  //std::cout << "draw point stub" << std::endl;
  this->setupFacade( );

  if (! Facade::instance().draw(p, coord))
  {
    return *this;
  }
  else
  {
    //std::cout << coord.first << " " << coord.second << std::endl;
    QGraphicsView* view = this->scene->views( ).first( );
    int height = view->height();
    // int width = view->width();

    double sceneRectWidth = this->scene->width();
    double sceneRectHeight = this->scene->height();

    QPoint coords( coord.first + sceneRectWidth/2,
                   height - coord.second -sceneRectHeight/2);
    QPointF qpt = view->mapToScene( coords );
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;
    this->qp->drawEllipse( qpt, radius, radius );
    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
  }
#if 0
  QPainter *ppnt = &ws.get_painter();
  QPen old_pen = ppnt->pen();
  ppnt->setPen(QPen(Qt::NoPen));
  unsigned sz = CGAL_REND_PT_RADIUS;
  ppnt->drawEllipse(coord.first - sz, ws.height() - coord.second - sz,
                      sz*2, sz*2);
  ppnt->setPen(old_pen);
#endif
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

  // std::cout<<"In setupFacade\n";
  // typedef Curve_renderer_facade<CKvA_2> Facade;

  // QGraphicsView* view = this->scene->views( ).first( );
  // QRect viewport = view->viewport()->rect();

  // CGAL::Bbox_2 bbox = this->convert( viewport ).bbox( );
  // Facade::setup(bbox, viewport.width(), viewport.height());

  // std::cout<<"Viewport\n";
  // std::cout<<viewport.left()<<"\t"<<viewport.right()<<std::endl;
  // std::cout<<viewport.top()<<"\t"<<viewport.bottom()<<std::endl;

  // // std::cout<<view->width()<<"\t"<<view->height()<<std::endl;
  // std::cout<<"BBox\n";
  // std::cout<<bbox.xmin()<<"\t"<<bbox.xmax()<<std::endl;
  // std::cout<<bbox.ymin()<<"\t"<<bbox.ymax()<<std::endl;

  // std::cout<<"Leaving setupFacade\n";
}

} // namespace Qt
} // namespace CGAL

// #endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_CPP
