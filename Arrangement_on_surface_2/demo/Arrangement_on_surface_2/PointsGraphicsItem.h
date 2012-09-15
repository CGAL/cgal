#ifndef POINTS_GRAPHICS_ITEM_H
#define POINTS_GRAPHICS_ITEM_H
#include <vector>
#include <CGAL/Qt/GraphicsItem.h>
#include <QPen>

class QPainter;
class QPen;

/**
Add a set of points to the QGraphicsScene.
*/
class PointsGraphicsItem: public CGAL::Qt::GraphicsItem
{
public:
    PointsGraphicsItem( );

    virtual void paint( QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget );
    virtual QRectF boundingRect( ) const;

    template < class Point >
    void insert( const Point& point )
    {
        this->prepareGeometryChange( );

        double x = CGAL::to_double( point.x( ) );
        double y = CGAL::to_double( point.y( ) );
        this->points.push_back( QPointF( x, y ) );
    }

    void clear( );

public slots:
    virtual void modelChanged( );

protected:
    std::vector< QPointF > points;
    QPen pointPen;

}; // class PointsGraphicsItem
#endif // POINTS_GRAPHICS_ITEM_H
