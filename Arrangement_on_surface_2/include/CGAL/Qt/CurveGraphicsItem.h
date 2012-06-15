#ifndef CGAL_QT_CURVE_GRAPHICS_ITEM_H
#define CGAL_QT_CURVE_GRAPHICS_ITEM_H
#include <CGAL/Qt/ArrangementPainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsItem.h>
namespace CGAL {
namespace Qt {

template < class ArrTraits >
class CurveGraphicsItem : public GraphicsItem
{
public:
    // known curve types
    typedef ArrTraits Traits;
    typedef typename Traits::Kernel Kernel;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

    CurveGraphicsItem( );

    virtual void paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget );
    void updateBoundingBox( );
    void insert( const X_monotone_curve_2& curve );
    void clear( );
    QRectF boundingRect( ) const;

public slots:
    void modelChanged( );

protected:
    CGAL::Qt::Converter< Kernel > convert;
    ArrangementPainterOstream< ArrTraits > painterOstream;
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
    painter->setPen( QPen( ::Qt::red, 1. ) );
    this->painterOstream = ArrangementPainterOstream< ArrTraits >( painter/*, clippingRectangle */ );
    for ( int i = 0; i < this->curves.size( ); ++i )
    {
        X_monotone_curve_2 curve = this->curves[ i ];
        this->painterOstream << curve;
    }
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

template < class ArrTraits >
QRectF 
CurveGraphicsItem< ArrTraits >::
boundingRect( ) const
{
    QRectF boundingRectangle = this->convert( this->boundingBox );
    return boundingRectangle;
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CURVE_GRAPHICS_ITEM_H
