#ifndef POINT_LOCATION_CALLBACK_H
#define POINT_LOCATION_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/CurveGraphicsItem.h>
#include <CGAL/Arrangement_with_history_2.h>
#include "Utils.h"

/**
Supports visualization of point location on arrangements.

The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class TArr >
class PointLocationCallback : public CGAL::Qt::Callback
{
public:
    typedef typename TArr::Halfedge_handle Halfedge_handle;
    typedef typename TArr::Halfedge_const_handle Halfedge_const_handle;
    typedef typename TArr::Halfedge_iterator Halfedge_iterator;
    typedef typename TArr::Face_handle Face_handle;
    typedef typename TArr::Face_const_handle Face_const_handle;
    typedef typename TArr::Vertex_const_handle Vertex_const_handle;
    typedef typename TArr::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
    typedef typename TArr::Geometry_traits_2 Traits;
    typedef typename TArr::Curve_handle Curve_handle;
    typedef typename TArr::Originating_curve_iterator Originating_curve_iterator;
    typedef typename TArr::Induced_edge_iterator Induced_edge_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::Point_2 Point;
    typedef typename Kernel::Segment_2 Segment;

    PointLocationCallback( TArr* arr_, QObject* parent_ );
    void setScene( QGraphicsScene* scene_ );
    QGraphicsScene* getScene( ) const;
    void reset( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    void highlightPointLocation( QGraphicsSceneMouseEvent *event );
    Face_const_handle getFace( const CGAL::Object& o );

    Compute_squared_distance_2< Traits > squaredDistance;
    CGAL::Qt::Converter< Kernel > convert;
    QGraphicsScene* scene;
    TArr* arr;
}; // class PointLocationCallback


template < class TArr >
PointLocationCallback< TArr >::
PointLocationCallback( TArr* arr_, QObject* parent_ ):
    CGAL::Qt::Callback( parent_ ),
    arr( arr_ ),
    scene( NULL )
{ }

template < class TArr >
void 
PointLocationCallback< TArr >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

template < class TArr >
QGraphicsScene* 
PointLocationCallback< TArr >::
getScene( ) const
{
    return this->scene;
}

template < class TArr >
void
PointLocationCallback< TArr >::
reset( )
{
    emit modelChanged( );
}

template < class TArr >
void 
PointLocationCallback< TArr >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    this->highlightPointLocation( event );
}

template < class TArr >
void 
PointLocationCallback< TArr >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{ }

template < class TArr >
void 
PointLocationCallback< TArr >::
highlightPointLocation( QGraphicsSceneMouseEvent* event )
{
    std::cout << "highlightPointLocation stub" << std::endl;
#if 0
    if (mode == MODE_POINT_LOCATION)
    {
      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(3);


      Point_2 temp_p (pl_point.x(), pl_point.y());
      CGAL::Object obj = locate(temp_p);

      Face_const_handle f = get_face(obj);
		
		/* more prudent color selection that selects the drawing color
		according to my_prefrance. replaced setColor(Qt::yellow)*/
		QColor my_preferance[4]= {Qt::yellow,Qt::green,Qt::red,Qt::blue};		
		setCorrectColor(f->color(),my_preferance, 4);		
		
      if (!f->is_unbounded()) // its an inside face
      {
        Ccb_halfedge_const_circulator cc = f->outer_ccb();
        do
        {
          m_tab_traits.draw_xcurve(this , cc->curve() );
        }
        while (++cc != f->outer_ccb());
      }


      //color the holes of the located face
      Holes_const_iterator hit, eit = f->holes_end();
      for (hit = f->holes_begin(); hit != eit; ++hit)
      {
        Ccb_halfedge_const_circulator cc = *hit;
        do
        {
          m_tab_traits.draw_xcurve(this , cc->curve() );
          cc++;
        }
        while (cc != *hit);
      }

      //color isolated vertices
      Isolated_vertex_const_iterator ivit = f->isolated_vertices_begin();
      for (; ivit != f->isolated_vertices_end(); ++ivit)
      {
        static_cast<CGAL::Qt_widget&>(*this) << ivit->point();
      }

      static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    }
#endif

    emit modelChanged( );
}

template < class TArr >
typename PointLocationCallback< TArr >::Face_const_handle
PointLocationCallback< TArr >::
getFace( const CGAL::Object& obj )
{
    Face_const_handle f;
    if ( CGAL::assign( f, obj ) )
        return f;

    Halfedge_const_handle he;
    if (CGAL::assign( he, obj ))
        return (he->face( ));

    Vertex_const_handle v;
    CGAL_assertion(CGAL::assign( v, obj ));
    CGAL::assign( v, obj );
    if ( v->is_isolated( ) )
        return v->face( );
    Halfedge_around_vertex_const_circulator eit = v->incident_halfedges( );
    return  (eit->face( ));
}

#endif // POINT_LOCATION_CALLBACK_H
