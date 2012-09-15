#include "PointsGraphicsItem.h"
#include <limits>
#include <QPen>
#include <QPainter>

PointsGraphicsItem::
PointsGraphicsItem( ):
    pointPen( QPen( ::Qt::black, 3. ) )
{ }

void
PointsGraphicsItem::
paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget )
{
    double scale = painter->worldTransform( ).m11( );
    double radius = this->pointPen.width( ) / 2.0;
    radius /= scale;
    QPen savePen = painter->pen( );
    painter->setBrush( this->pointPen.brush( ) );

    for ( int i = 0; i < this->points.size( ); ++i )
    {
        QPointF pt = this->points[ i ];
        painter->drawEllipse( pt, radius, radius );
    }

    painter->setPen( savePen );
}

QRectF
PointsGraphicsItem::
boundingRect( ) const
{
    if ( this->points.size( ) == 0 )
    {
        return QRectF( );
    }
    double xmin = std::numeric_limits< double >::max( );
    double xmax = -std::numeric_limits< double >::max( );
    double ymin = std::numeric_limits< double >::max( );
    double ymax = -std::numeric_limits< double >::max( );
    for ( int i = 0; i < this->points.size( ); ++i )
    {
        QPointF pt = this->points[ i ];
        double x = pt.x( );
        double y = pt.y( );
        xmin = std::min( xmin, x );
        xmax = std::max( xmax, x );
        ymin = std::min( ymin, y );
        ymax = std::max( ymax, y );
    }
    QRectF res( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );
    res.adjust( -5, -5, 5, 5 ); // pad the borders a bit
    return res;
}

void
PointsGraphicsItem::
clear( )
{
    this->prepareGeometryChange( );

    this->points.clear( );
}

void
PointsGraphicsItem::
modelChanged( )
{
    if ( this->points.size( ) == 0 )
    {
        this->hide( );
    }
    else
    {
        this->show( );
    }
    this->update( );
}
