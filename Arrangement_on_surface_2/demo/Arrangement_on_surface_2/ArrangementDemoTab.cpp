#include "ArrangementDemoTab.h"

ArrangementDemoTabBase::
ArrangementDemoTabBase( QWidget* parent ):
    QWidget( parent ),
    layout( new QGridLayout( this ) ),
    graphicsView( new ArrangementDemoGraphicsView( this ) ),
    scene( new QGraphicsScene( -100, -100, 100, 100 ) ),
    arrangementGraphicsItem( NULL ),
    segmentInputCallback( NULL ),
    deleteCurveCallback( NULL ),
    pointLocationCallback( NULL ),
    verticalRayShootCallback( NULL ),
    mergeEdgeCallback( NULL ),
    splitEdgeCallback( NULL ),
    envelopeCallback( NULL )
{
    this->setupUi( );
}

void
ArrangementDemoTabBase::
setupUi( )
{
    this->layout->addWidget( this->graphicsView, 0, 0, 1, 1 );
    this->graphicsView->setScene( this->scene );
    this->graphicsView->setMouseTracking( true );
}

QGraphicsScene* 
ArrangementDemoTabBase::
getScene( ) const
{
    return this->scene;
}

QGraphicsView* 
ArrangementDemoTabBase::
getView( ) const
{
    return this->graphicsView;
}

CGAL::Qt::GraphicsItem*
ArrangementDemoTabBase::
getArrangementGraphicsItem( ) const
{
    return this->arrangementGraphicsItem;
}

CGAL::Qt::GraphicsViewSegmentInputBase*
ArrangementDemoTabBase::
getSegmentInputCallback( ) const
{
    return this->segmentInputCallback;
}

CGAL::Qt::Callback*
ArrangementDemoTabBase::
getDeleteCurveCallback( ) const
{
    return this->deleteCurveCallback;
}

CGAL::Qt::Callback*
ArrangementDemoTabBase::
getPointLocationCallback( ) const
{
    return this->pointLocationCallback;
}

VerticalRayShootCallbackBase*
ArrangementDemoTabBase::
getVerticalRayShootCallback( ) const
{
    return this->verticalRayShootCallback;
}

CGAL::Qt::Callback*
ArrangementDemoTabBase::
getMergeEdgeCallback( ) const
{
    return this->mergeEdgeCallback;
}

SplitEdgeCallbackBase*
ArrangementDemoTabBase::
getSplitEdgeCallback( ) const
{
    return this->splitEdgeCallback;
}

EnvelopeCallbackBase*
ArrangementDemoTabBase::
getEnvelopeCallback( ) const
{
    return this->envelopeCallback;
}
