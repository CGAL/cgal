#ifndef ARRANGEMENT_CURVE_INPUT_CALLBACK_H
#define ARRANGEMENT_CURVE_INPUT_CALLBACK_H
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QEvent>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
#include "GraphicsViewCurveInput.h"
#include "Utils.h"

template < class Arr_, class ArrTraits = typename Arr_::Geometry_traits_2 >
class ArrangementCurveInputCallback:
    public CGAL::Qt::GraphicsViewCurveInput< typename Arr_::Geometry_traits_2 >
{
public:
    typedef Arr_ Arrangement;
    typedef ArrTraits Traits;
    typedef CGAL::Qt::GraphicsViewCurveInput< Traits > Superclass;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::FT FT;

    ArrangementCurveInputCallback( Arrangement* arrangement_, QObject* parent ):
        Superclass( parent ),
        arrangement( arrangement_ )
    {
        this->snapToVertexStrategy.setArrangement( arrangement_ );

        QObject::connect( this, SIGNAL( generate( CGAL::Object ) ),
            this, SLOT( processInput( CGAL::Object ) ) );
    }

    void processInput( CGAL::Object o )
    {
        Curve_2 curve;
        if ( CGAL::assign( curve, o ) )
        {
            CGAL::insert( *( this->arrangement ), curve );
        }
    
        emit CGAL::Qt::GraphicsViewInput::modelChanged( );
    }

    void setScene( QGraphicsScene* scene )
    {
        this->Superclass::setScene( scene );
        this->snapToVertexStrategy.setScene( scene );
        this->snapToGridStrategy.setScene( scene );
    }

protected:
    Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        if ( this->snapToGridEnabled )
        {
            return this->snapToGridStrategy.snapPoint( event );
        }
        else if ( this->snappingEnabled )
        {
            return this->snapToVertexStrategy.snapPoint( event );
        }
        else
        {
            return this->convert( event->scenePos( ) );
        }
    }

    Arrangement* arrangement;
    SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
    SnapToGridStrategy< Kernel > snapToGridStrategy;
}; // class ArrangementCurveInputCallback

#if 0
template < class Arr_, class Kernel_ >
class ArrangementCurveInputCallback< Arr_, CGAL::Arr_segment_traits_2< Kernel_ > >:
    public CGAL::Qt::GraphicsViewCurveInput< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef CGAL::Qt::GraphicsViewCurveInput< Traits > Superclass;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Traits::Curve_2 Curve_2;
    typedef Kernel_ Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::FT FT;

    ArrangementCurveInputCallback( Arrangement* arrangement_, QObject* parent ):
        Superclass( parent ),
        arrangement( arrangement_ )
    {
        this->snapToVertexStrategy.setArrangement( arrangement_ );

        QObject::connect( this, SIGNAL( generate( CGAL::Object ) ),
            this, SLOT( processInput( CGAL::Object ) ) );
    }

    void processInput( CGAL::Object o )
    {
        Curve_2 curve;
        if ( CGAL::assign( curve, o ) )
        {
            CGAL::insert( *( this->arrangement ), curve );
        }
    
        emit CGAL::Qt::GraphicsViewInput::modelChanged( );
    }

    void setScene( QGraphicsScene* scene )
    {
        this->Superclass::setScene( scene );
        this->snapToVertexStrategy.setScene( scene );
        this->snapToGridStrategy.setScene( scene );
    }

protected:
    Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        if ( this->snapToGridEnabled )
        {
            return this->snapToGridStrategy.snapPoint( event );
        }
        else if ( this->snappingEnabled )
        {
            return this->snapToVertexStrategy.snapPoint( event );
        }
        else
        {
            return this->convert( event->scenePos( ) );
        }
    }

    Arrangement* arrangement;
    SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
    SnapToGridStrategy< Kernel > snapToGridStrategy;
}; // class ArrangementCurveInputCallback< Arr_, CGAL::Arr_segment_traits_2< Kernel_ > >
#endif

#endif // ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
