#ifndef ENVELOPE_CALLBACK_H
#define ENVELOPE_CALLBACK_H
#include "Callback.h"
#include <CGAL/Qt/CurveGraphicsItem.h>
#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <list>
#include "Utils.h"

class EnvelopeCallbackBase : public CGAL::Qt::Callback
{
public:

public slots:
    virtual void showLowerEnvelope( bool b ) = 0;
    virtual void showUpperEnvelope( bool b ) = 0;

protected:
    EnvelopeCallbackBase( QObject* parent );
}; // class EnvelopeCallbackBase

template < class Arr_ >
class EnvelopeCallback : public EnvelopeCallbackBase
{
public:
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Edge_iterator Edge_iterator;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef CGAL::Envelope_diagram_1< Traits > Diagram_1;

    EnvelopeCallback( Arrangement* arr_, QObject* parent );
    void showLowerEnvelope( bool show );
    void showUpperEnvelope( bool show );
    void slotModelChanged( );
    void setScene( QGraphicsScene* scene_ );

protected:
    void updateEnvelope( bool lower );
    Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;

    Arrangement* arr;
    CGAL::Qt::CurveGraphicsItem< Traits >* lowerEnvelope;
    CGAL::Qt::CurveGraphicsItem< Traits >* upperEnvelope;
    using CGAL::Qt::Callback::scene;
}; // class EnvelopeCallback

template < class Arr_ >
EnvelopeCallback< Arr_ >::
EnvelopeCallback( Arrangement* arr_, QObject* parent ):
    EnvelopeCallbackBase( parent ),
    arr( arr_ ),
    lowerEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
    upperEnvelope( new CGAL::Qt::CurveGraphicsItem< Traits >( ) )
{
    this->lowerEnvelope->hide( );
    this->upperEnvelope->hide( );
}

template < class Arr_ >
void
EnvelopeCallback< Arr_ >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

template < class Arr_ >
void
EnvelopeCallback< Arr_ >::
slotModelChanged( )
{
    this->updateEnvelope( true );
    this->updateEnvelope( false );
}

template < class Arr_ >
void
EnvelopeCallback< Arr_ >::
updateEnvelope( bool lower )
{
    CGAL::Qt::CurveGraphicsItem< Traits >* envelopeToUpdate;
    if ( lower )
    {
        envelopeToUpdate = this->lowerEnvelope;
    }
    else
    {
        envelopeToUpdate = this->upperEnvelope;
    }
    envelopeToUpdate->clear( );

    std::list< X_monotone_curve_2 > curves;
    for ( Edge_iterator eit = this->arr->edges_begin( ); eit != this->arr->edges_end( ); ++eit )
    {
        curves.push_back( eit->curve( ) );
    }
    Diagram_1 diagram;
    if ( lower )
    {
        CGAL::lower_envelope_x_monotone_2( curves.begin( ), curves.end( ), diagram );
    }
    else
    {
        CGAL::upper_envelope_x_monotone_2( curves.begin( ), curves.end( ), diagram );
    }

    typename Diagram_1::Edge_const_handle e = diagram.leftmost( );
    typename Diagram_1::Vertex_const_handle v;
    while ( e != diagram.rightmost( ) )
    {
        if ( ! e->is_empty( ) )
        {
            // The edge is not empty: draw a representative curve.
            // Note that the we only draw the portion of the curve
            // that overlaps the x-range defined by the two vertices
            // that are incident to this edge.

            // TODO: generate a subcurve instead of just making a segment

            Point_2 leftPoint = e->left( )->point( );
            Point_2 rightPoint = e->right( )->point( );
            X_monotone_curve_2 curve =
                this->construct_x_monotone_subcurve_2( e->curve( ), leftPoint, rightPoint );
            envelopeToUpdate->insert( curve );
        }
        v = e->right( );

        // TODO: Draw the point associated with the current vertex.
        e = v->right( );
    }
    envelopeToUpdate->modelChanged( );
}

template < class Arr_ >
void
EnvelopeCallback< Arr_ >::
showLowerEnvelope( bool show )
{
    if ( show )
    {
        std::cout << "Show lower envelope" << std::endl;
        this->scene->addItem( this->lowerEnvelope );
    }
    else
    {
        std::cout << "Hide lower envelope" << std::endl;
        this->scene->removeItem( this->lowerEnvelope );
    }
}

template < class Arr_ >
void
EnvelopeCallback< Arr_ >::
showUpperEnvelope( bool show )
{
    if ( show )
    {
        std::cout << "Show upper envelope" << std::endl;
        this->scene->addItem( this->upperEnvelope );
    }
    else
    {
        std::cout << "Hide upper envelope" << std::endl;
        this->scene->removeItem( this->upperEnvelope );
    }
}
#endif // ENVELOPE_CALLBACK_H
