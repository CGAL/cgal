#ifndef VERTICAL_RAY_SHOOT_CALLBACK_H
#define VERTICAL_RAY_SHOOT_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include "CurveGraphicsItem.h"
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <vector>

#include "Utils.h"

class VerticalRayShootCallbackBase : public CGAL::Qt::Callback
{
public:
    void setShootingUp( bool isShootingUp );

protected:
    VerticalRayShootCallbackBase( QObject* parent_ );
    using Callback::scene;
    bool shootingUp;
}; // class VerticalRayShootCallbackBase

/**
Supports visualization of vertical ray shooting on arrangements.

The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class Arr_ >
class VerticalRayShootCallback : public VerticalRayShootCallbackBase
{
public:
    typedef VerticalRayShootCallbackBase Superclass;
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Halfedge_handle Halfedge_handle;
    typedef typename Arrangement::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
    typedef typename Arrangement::Face_handle Face_handle;
    typedef typename Arrangement::Face_const_handle Face_const_handle;
    typedef typename Arrangement::Vertex_const_handle Vertex_const_handle;
    typedef typename Arrangement::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Arrangement::Curve_handle Curve_handle;
    typedef typename Arrangement::Originating_curve_iterator Originating_curve_iterator;
    typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
    typedef typename Arrangement::Ccb_halfedge_const_circulator Ccb_halfedge_const_circulator;
    typedef typename Arrangement::Hole_const_iterator Hole_const_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
    typedef typename Traits::Intersect_2 Intersect_2;
    typedef typename Traits::Multiplicity Multiplicity;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef std::pair< typename Traits::Point_2, Multiplicity > IntersectionResult;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::FT FT;
    typedef typename CGAL::Arr_trapezoid_ric_point_location< Arrangement > TrapezoidPointLocationStrategy;
    typedef typename CGAL::Arr_simple_point_location< Arrangement > SimplePointLocationStrategy;
    typedef typename CGAL::Arr_walk_along_line_point_location< Arrangement > WalkAlongLinePointLocationStrategy;
    typedef typename CGAL::Arr_landmarks_point_location< Arrangement > LandmarksPointLocationStrategy;

    VerticalRayShootCallback( Arrangement* arr_, QObject* parent_ );
    void reset( );
    void setScene( QGraphicsScene* scene_ );

    void slotModelChanged( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    void highlightPointLocation( QGraphicsSceneMouseEvent *event );
    Face_const_handle getFace( const CGAL::Object& o );
    CGAL::Object rayShootUp( const Point_2& point );
    CGAL::Object rayShootDown( const Point_2& point );

    using Superclass::scene;
    using Superclass::shootingUp;
    Traits traits;
    Arrangement* arr;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
    Intersect_2 intersectCurves;
    CGAL::Qt::Converter< Kernel > convert;
    CGAL::Object pointLocationStrategy;
    CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
    QGraphicsLineItem* activeRay;
}; // class VerticalRayShootCallback

template < class Arr_ >
VerticalRayShootCallback< Arr_ >::
VerticalRayShootCallback( Arrangement* arr_, QObject* parent_ ):
    VerticalRayShootCallbackBase( parent_ ),
    arr( arr_ ),
    construct_x_monotone_curve_2( this->traits.construct_x_monotone_curve_2_object( ) ),
    intersectCurves( this->traits.intersect_2_object( ) ),
    pointLocationStrategy( CGAL::make_object( new WalkAlongLinePointLocationStrategy( *arr_ ) ) ),
    highlightedCurves( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
    activeRay( new QGraphicsLineItem )
{
    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this->highlightedCurves, SLOT( modelChanged( ) ) );
    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this, SLOT( slotModelChanged( ) ) );
}

template < class Arr_ >
void
VerticalRayShootCallback< Arr_ >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
    if ( this->scene )
    {
        this->scene->addItem( this->highlightedCurves );
        this->scene->addItem( this->activeRay );
    }
}


template < class Arr_ >
void
VerticalRayShootCallback< Arr_ >::
slotModelChanged( )
{
    this->activeRay->update( );
}

template < class Arr_ >
void
VerticalRayShootCallback< Arr_ >::
reset( )
{
    this->activeRay->setLine( 0, 0, 0, 0 );
    this->highlightedCurves->clear( );
    emit modelChanged( );
}

template < class Arr_ >
void 
VerticalRayShootCallback< Arr_ >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    this->highlightPointLocation( event );
}

template < class Arr_ >
void 
VerticalRayShootCallback< Arr_ >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{ }

template < class Arr_ >
void 
VerticalRayShootCallback< Arr_ >::
highlightPointLocation( QGraphicsSceneMouseEvent* event )
{
    this->highlightedCurves->clear( );
    Point_2 p1 = this->convert( event->scenePos( ) );
    CGAL::Object pointLocationResult;
    if ( this->shootingUp )
    {
        pointLocationResult = this->rayShootUp( p1 );
    }
    else
    {
        pointLocationResult = this->rayShootDown( p1 );
    }
    if ( pointLocationResult.is_empty( ) )
    {
        return;
    }
    
    QRectF viewportRect = this->viewportRect( );
    FT y2;
    if ( this->shootingUp )
    { // +y in Qt is towards the bottom
        y2 = FT( viewportRect.bottom( ) );
    }
    else
    {
        y2 = FT( viewportRect.top( ) );
    }
    Face_const_handle unboundedFace;
    Halfedge_const_handle halfedge;
    Vertex_const_handle vertex;
    if ( CGAL::assign( unboundedFace, pointLocationResult ) )
    {
        Point_2 p2( FT( p1.x( ) ), y2 );
        Segment_2 lineSegment( p1, p2 );
        QLineF qLineSegment = this->convert( lineSegment );
        this->activeRay->setLine( qLineSegment );
    }
    else if ( CGAL::assign( halfedge, pointLocationResult ) )
    {
        this->highlightedCurves->insert( halfedge->curve( ) );
        Point_2 p1c1( p1.x( ), p1.y( ) );
        Point_2 p2c1( p1.x( ), y2 );
        //std::cout << "p1.x( ): " << p1.x( ) << std::endl;
        //std::cout << "y2: " << y2 << std::endl;
        const X_monotone_curve_2 c1 =
            this->construct_x_monotone_curve_2( p1c1, p2c1 );
        const X_monotone_curve_2 c2 = halfedge->curve( );

        CGAL::Object res;
        CGAL::Oneset_iterator< CGAL::Object > oi( res );
        //std::vector< CGAL::Object > ress;

        this->intersectCurves( c1, c2, oi );
        //this->intersectCurves( c1, c2, std::back_inserter( ress ) );
        //std::cout << "num intersections: " << ress.size( ) << std::endl;
        //typedef typename 
        //std::pair< Point_2, Multiplicity > pair;
        IntersectionResult pair;
        if ( CGAL::assign( pair, res ) )
        {
            //std::cout << "bound to IntersectionResult" << std::endl;
            Point_2 p2 = pair.first;
            Segment_2 lineSegment( p1, p2 );
            QLineF qLineSegment = this->convert( lineSegment );
            this->activeRay->setLine( qLineSegment );
        }
        else
        {
            Point_2 p2 = pair.first;
            Segment_2 lineSegment( p1, p2c1 );
            QLineF qLineSegment = this->convert( lineSegment );
            this->activeRay->setLine( qLineSegment );
        }
    }
    else if ( CGAL::assign( vertex, pointLocationResult ) )
    {
        std::cout << "Hit a vertex. Now that's rare." << std::endl;
    }

    emit modelChanged( );
#if 0
    {
        Face_const_handle ubf;
        if (CGAL::assign(ubf, obj))
        {
            CGAL_assertion(ubf->is_unbounded());
            //relevant_face_color = unbounded_face_color() as initialized          
            up = Coord_point(pl_draw.x() , y_max());
            static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw, up);
        }
        // we shoot something
        else
        {
            Halfedge_const_handle he;
            if (CGAL::assign(he, obj))
            {
                Point_2 p1c1(pl_point.x() , y_max() * m_tab_traits.COORD_SCALE);
                Point_2 p2c1(pl_point.x() , pl_point.y());
                const X_monotone_curve_2 c1 =
                    m_tab_traits.curve_make_x_monotone(p1c1 , p2c1);
                const X_monotone_curve_2 c2 = he->curve();

                CGAL::Object             res;
                CGAL::Oneset_iterator<CGAL::Object> oi(res);

                m_traits.intersect_2_object()(c1, c2, oi);
                std::pair<Point_2,Multiplicity> p1;
                if (CGAL::assign(p1, res))
                {
                    Coord_type y1 =
                        CGAL::to_double(p1.first.y())/ m_tab_traits.COORD_SCALE;
                    up = Coord_point(pl_draw.x(), y1);
                }
                else
                {
                    up = pl_draw;
                }
                relevant_face_color = he->face()->color();				
                /*choose color to mark the edge that differs from the current 
                  edge_color, the background, and the relevant face color*/				
                setCorrectColor(relevant_face_color);				 
                m_tab_traits.draw_xcurve(this , he->curve() );
            }
            else
            {
                Vertex_const_handle v;
                CGAL_assertion(CGAL::assign(v, obj));
                CGAL::assign(v, obj);
                up = Coord_point(CGAL::to_double(v->point().x()) /
                        m_tab_traits.COORD_SCALE,
                        CGAL::to_double(v->point().y()) /
                        m_tab_traits.COORD_SCALE);

                //locate face that arrow will be drawn in, and retrieve its color 
                CGAL::Object obj1 = locate(temp_p);
                Face_const_handle f1 = get_face(obj1);
                relevant_face_color=f1->color();

                /*choose color to mark the vertice so that it differs from the 
                  edge_color, the background, and the relevant_face_color*/				
                setCorrectColor(relevant_face_color);				     
                static_cast<CGAL::Qt_widget&>(*this) << up;
            }
        }

        //select arrow color that differs from the color of the face it is in
        setCorrectColor(relevant_face_color);        

        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(2);
        static_cast<CGAL::Qt_widget&>(*this) << Coord_segment(pl_draw,up);

        // draw an arrow that points to 'up' point
        int x = this->x_pixel(CGAL::to_double(up.x()));
        int y = this->y_pixel(CGAL::to_double(up.y()));

        this->get_painter().drawLine(x-7 , y+7 , x , y);
        this->get_painter().drawLine(x+7 , y+7 , x , y);
        static_cast<CGAL::Qt_widget&>(*this) << CGAL::LineWidth(m_line_width);
    }
#endif
}

template < class Arr_ >
typename VerticalRayShootCallback< Arr_ >::Face_const_handle
VerticalRayShootCallback< Arr_ >::
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

template < class Arr_ >
CGAL::Object
VerticalRayShootCallback< Arr_ >::
rayShootUp( const Point_2& point )
{
    CGAL::Object pointLocationResult;
    WalkAlongLinePointLocationStrategy* walkStrategy;
    TrapezoidPointLocationStrategy* trapezoidStrategy;
    SimplePointLocationStrategy* simpleStrategy;
    LandmarksPointLocationStrategy* landmarksStrategy;
    if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = walkStrategy->ray_shoot_up( point );
    }
    else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = trapezoidStrategy->ray_shoot_up( point );
    }
    else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = simpleStrategy->ray_shoot_up( point );
    }
    else if ( CGAL::assign( landmarksStrategy, this->pointLocationStrategy ) )
    {
        // pointLocationResult = landmarksStrategy->locate( point );
        std::cerr << "Warning: landmarks point location strategy doesn't support ray shooting" << std::endl;
        return CGAL::Object( );
    }
    return pointLocationResult;
}

template < class Arr_ >
CGAL::Object
VerticalRayShootCallback< Arr_ >::
rayShootDown( const Point_2& point )
{
    CGAL::Object pointLocationResult;
    WalkAlongLinePointLocationStrategy* walkStrategy;
    TrapezoidPointLocationStrategy* trapezoidStrategy;
    SimplePointLocationStrategy* simpleStrategy;
    LandmarksPointLocationStrategy* landmarksStrategy;
    if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = walkStrategy->ray_shoot_down( point );
    }
    else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = trapezoidStrategy->ray_shoot_down( point );
    }
    else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
    {
        pointLocationResult = simpleStrategy->ray_shoot_down( point );
    }
    else if ( CGAL::assign( landmarksStrategy, this->pointLocationStrategy ) )
    {
        // pointLocationResult = landmarksStrategy->locate( point );
        std::cerr << "Warning: landmarks point location strategy doesn't support ray shooting" << std::endl;
        return CGAL::Object( );
    }
    return pointLocationResult;
}


#endif // VERTICAL_RAY_SHOOT_CALLBACK_H
