#ifndef CGAL_QT_CURVE_GRAPHICS_ITEM_H
#define CGAL_QT_CURVE_GRAPHICS_ITEM_H
#include "ArrangementPainterOstream.h"
#include "Utils.h"
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <QGraphicsScene>

namespace CGAL {
namespace Qt {

template < class ArrTraits >
class CurveGraphicsItem : public GraphicsItem, public QGraphicsSceneMixin
{
public:
    // known curve types
    typedef ArrTraits Traits;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

public: // ctors
    CurveGraphicsItem( );

public: // methods
    virtual void paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget );
    virtual QRectF boundingRect( ) const;
    void insert( const X_monotone_curve_2& curve );
    void clear( );

public slots:
    void modelChanged( );

protected: // methods
    void updateBoundingBox( );

protected: // fields
    CGAL::Qt::Converter< Kernel > convert;
    ArrangementPainterOstream< Traits > painterOstream;
    std::vector< X_monotone_curve_2 > curves;
    CGAL::Bbox_2 boundingBox;  
    bool boundingBoxInitialized;
}; // class CurveGraphicsItem

template < class ArrTraits >
CurveGraphicsItem< ArrTraits >::
CurveGraphicsItem( ):
    painterOstream( 0 ),
    boundingBox( 0, 0, 0, 0 ),
    boundingBoxInitialized( false )
{
    this->setZValue( 4 );
}

template < class ArrTraits >
void
CurveGraphicsItem< ArrTraits >::
paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget )
{
    painter->setPen( QPen( ::Qt::red, 0. ) );
    QRectF clippingRectangle = this->viewportRect( );
    this->painterOstream = ArrangementPainterOstream< Traits >( painter, clippingRectangle );
    for ( int i = 0; i < this->curves.size( ); ++i )
    {
        X_monotone_curve_2 curve = this->curves[ i ];
        this->painterOstream << curve;
    }
}

template < class ArrTraits >
QRectF 
CurveGraphicsItem< ArrTraits >::
boundingRect( ) const
{
    QRectF boundingRectangle = this->convert( this->boundingBox );
    return boundingRectangle;
}

template < class ArrTraits >
void 
CurveGraphicsItem< ArrTraits >::
modelChanged( )
{
    if ( this->curves.size( ) == 0 )
    {
        this->hide( );
    }
    else
    {
        this->show( );
    }
    this->updateBoundingBox( );
    this->update( );
}

template < class ArrTraits >
void 
CurveGraphicsItem< ArrTraits >::
updateBoundingBox( )
{
    this->prepareGeometryChange( );
    if ( this->curves.size( ) == 0 )
    {
        this->boundingBox = Bbox_2( 0, 0, 0, 0 );
        this->boundingBoxInitialized = false;
        return;
    }
    else
    {
        this->boundingBox = this->curves[ 0 ].bbox( );
        this->boundingBoxInitialized = true;
    }

    for ( int i = 1; i < this->curves.size( ); ++i )
    {
        this->boundingBox = this->boundingBox + this->curves[ i ].bbox( );
    }
}

template < class ArrTraits >
void 
CurveGraphicsItem< ArrTraits >::
insert( const X_monotone_curve_2& segment )
{
    this->curves.push_back( segment );
}

template < class ArrTraits >
void 
CurveGraphicsItem< ArrTraits >::
clear( )
{
    this->curves.clear( );
}

/**
Specialization of the base template CurveGraphicsItem:

    updateBoundingBox
*/
template < class Kernel_ >
class CurveGraphicsItem< CGAL::Arr_linear_traits_2< Kernel_ > > : public GraphicsItem, public QGraphicsSceneMixin
{
public: // typedefs
    // known curve types
    typedef CGAL::Arr_linear_traits_2< Kernel_ > Traits;
    typedef Kernel_ Kernel;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Ray_2 Ray_2;

public: // ctors
    CurveGraphicsItem( ):
        painterOstream( 0 ),
        boundingBox( 0, 0, 0, 0 ),
        boundingBoxInitialized( false )
    {
        this->setZValue( 4 );
    }

public: // methods
    virtual void paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget )
    {
        QRectF clippingRectangle = this->viewportRect( );
        painter->setPen( QPen( ::Qt::red, 0. ) );
        this->painterOstream = ArrangementPainterOstream< Traits >( painter, clippingRectangle );
        for ( int i = 0; i < this->curves.size( ); ++i )
        {
            X_monotone_curve_2 curve = this->curves[ i ];
            this->painterOstream << curve;
        }
    }

    void insert( const X_monotone_curve_2& curve )
    {
        this->curves.push_back( curve );
    }

    void clear( )
    {
        this->curves.clear( );
    }

    QRectF boundingRect( ) const
    {
        QRectF res;
        if ( this->getScene( ) == NULL )
        {
            return res;
        }
    }

    virtual void setScene( QGraphicsScene* scene_ )
    {
        this->QGraphicsSceneMixin::setScene( scene_ );
        if ( this->getScene( ) == NULL )
        {
            QRectF clipRect = this->viewportRect( );
            this->convert = CGAL::Qt::Converter< Kernel >( clipRect );
        }
    }

public slots:
    void modelChanged( )
    {
        if ( this->curves.size( ) == 0 )
        {
            this->hide( );
        }
        else
        {
            this->show( );
        }
        this->updateBoundingBox( );
        this->update( );
    }

protected: // methods
    void updateBoundingBox( )
    {
        this->boundingBoxInitialized = 0;
        this->boundingBox = CGAL::Bbox_2( 0, 0, 0, 0 );
        QRectF clipRect = this->viewportRect( );
        if ( !clipRect.isValid( ) )
        {
            return;
        }
        this->convert = CGAL::Qt::Converter< Kernel >( clipRect );

        bool first = 1;
        for ( int i = 0; i < curves.size( ); ++i )
        {
            X_monotone_curve_2 curve = curves[ i ];
            if ( curve.is_segment( ) )
            {
                Segment_2 seg = curve.segment( );
                CGAL::Bbox_2 seg_bbox = seg.bbox( );
                if ( first )
                {
                    first = 0;
                    this->boundingBoxInitialized = 1;
                    this->boundingBox = seg_bbox;
                }
                else
                {
                    this->boundingBox = this->boundingBox + seg_bbox;
                }
            }
            else if ( curve.is_ray( ) )
            {
                Ray_2 ray = curve.ray( );
                QLineF qclippedRay = this->convert( ray );
                if ( qclippedRay.isNull( ) )
                    continue;
                Segment_2 clippedRay = this->convert( qclippedRay );
                if ( first )
                {
                    first = 0;
                    this->boundingBoxInitialized = 1;
                    this->boundingBox = clippedRay.bbox( );
                }
                else
                {
                    this->boundingBox = this->boundingBox + clippedRay.bbox( );
                }
            }
            else // curve.is_line( )
            {
                Line_2 line = curve.line( );
                QLineF qclippedLine = this->convert( line );
                if ( qclippedLine.isNull( ) )
                    continue;
                Segment_2 clippedLine = this->convert( qclippedLine );
                if ( first )
                {
                    first = 0;
                    this->boundingBoxInitialized = 1;
                    this->boundingBox = clippedLine.bbox( );
                }
                else
                {
                    this->boundingBox = this->boundingBox + clippedLine.bbox( );
                }
            }
        }
    }

protected: // fields
    CGAL::Qt::Converter< Kernel > convert;
    ArrangementPainterOstream< Traits > painterOstream;
    std::vector< X_monotone_curve_2 > curves;
    CGAL::Bbox_2 boundingBox;  
    bool boundingBoxInitialized;
}; // class CurveGraphicsItem

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CURVE_GRAPHICS_ITEM_H
