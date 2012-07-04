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
    typedef typename ArrTraitsAdaptor< Traits >::Point_2 Point_2;
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

    void setArrangement( Arrangement* newArr )
    {
        this->arrangement = newArr;
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

#endif // ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
