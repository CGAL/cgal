#include "ArrangementDemoGraphicsView.h"
#include <iostream>
#include <QVarLengthArray>

ArrangementDemoGraphicsView::
ArrangementDemoGraphicsView( QWidget* parent ):
    QGraphicsView( parent ),
    showGrid( false ),
    gridSize( 50 )
{ }


void
ArrangementDemoGraphicsView::
setShowGrid( bool b )
{
    this->showGrid = b;
}

bool
ArrangementDemoGraphicsView::
getShowGrid( ) const
{
    return this->showGrid;
}

void 
ArrangementDemoGraphicsView::
setGridSize( int size )
{
    this->gridSize = size;
}

int 
ArrangementDemoGraphicsView::
getGridSize( ) const
{
    return this->gridSize;
}

void
ArrangementDemoGraphicsView::
drawForeground( QPainter* painter, const QRectF& rect )
{
    QRectF viewportRect = this->getViewportRect( );
    if ( this->showGrid )
    {
        QVarLengthArray< QLineF, 100 > linesX;
        QVarLengthArray< QLineF, 100 > linesY;
        qreal left = int(viewportRect.left()) - (int(viewportRect.left()) % this->gridSize);
        qreal top = int(viewportRect.top()) - (int(viewportRect.top()) % this->gridSize);
        for ( qreal x = left; x < viewportRect.right( ); x += this->gridSize )
        {
            linesX.append( QLineF( x, viewportRect.top( ), x, viewportRect.bottom( ) ) );
        }
        for ( qreal y = top; y < viewportRect.bottom( ); y += this->gridSize )
        {
            linesY.append( QLineF( viewportRect.left( ), y, viewportRect.right( ), y ) );
        }
        painter->drawLines( linesX.data( ), linesX.size( ) );
        painter->drawLines( linesY.data( ), linesY.size( ) );
    }
}

QRectF
ArrangementDemoGraphicsView::
getViewportRect( ) const
{
    QPointF p1 = this->mapToScene( 0, 0 );
    QPointF p2 = this->mapToScene( this->width( ), this->height( ) );
    QRectF res = QRectF( p1, p2 );
    return res;
}
