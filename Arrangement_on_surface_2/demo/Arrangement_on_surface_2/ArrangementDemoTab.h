#ifndef ARRANGEMENT_DEMO_TAB_H
#define ARRANGEMENT_DEMO_TAB_H
#include <QWidget>
#include <QGridLayout>

#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include "ArrangementDemoGraphicsView.h"
#include "ArrangementSegmentInputCallback.h"
#include "DeleteCurveCallback.h"
#include "PointLocationCallback.h"
#include "VerticalRayShootCallback.h"
#include "MergeEdgeCallback.h"
#include "SplitEdgeCallback.h"
#include "EnvelopeCallback.h"

class ArrangementDemoTabBase : public QWidget
{
public:
    ArrangementDemoTabBase( QWidget* parent );

    virtual QGraphicsScene* getScene( ) const;
    virtual QGraphicsView* getView( ) const;

    virtual CGAL::Qt::GraphicsItem* getArrangementGraphicsItem( ) const;
    virtual CGAL::Qt::GraphicsViewSegmentInputBase* getSegmentInputCallback( ) const;
    virtual CGAL::Qt::Callback* getDeleteCurveCallback( ) const;
    virtual CGAL::Qt::Callback* getPointLocationCallback( ) const;
    virtual VerticalRayShootCallbackBase* getVerticalRayShootCallback( ) const;
    virtual CGAL::Qt::Callback* getMergeEdgeCallback( ) const;
    virtual SplitEdgeCallbackBase* getSplitEdgeCallback( ) const;
    virtual EnvelopeCallbackBase* getEnvelopeCallback( ) const;

protected:
    virtual void setupUi( );

    ArrangementDemoGraphicsView* graphicsView;
    QGraphicsScene* scene;
    QGridLayout* layout;

    CGAL::Qt::GraphicsItem* arrangementGraphicsItem;
    CGAL::Qt::GraphicsViewSegmentInputBase* segmentInputCallback;
    CGAL::Qt::Callback* deleteCurveCallback;
    CGAL::Qt::Callback* pointLocationCallback;
    VerticalRayShootCallbackBase* verticalRayShootCallback;
    CGAL::Qt::Callback* mergeEdgeCallback;
    SplitEdgeCallbackBase* splitEdgeCallback;
    EnvelopeCallbackBase* envelopeCallback;

}; // class ArrangementDemoTabBase

template < class Arr_ >
class ArrangementDemoTab : public ArrangementDemoTabBase
{
public:
    typedef ArrangementDemoTabBase Superclass;
    typedef Arr_ Arrangement;

    ArrangementDemoTab( Arrangement* arrangement_, QWidget* parent ):
        Superclass( parent ),
        arrangement( arrangement_ )
    {
        this->arrangementGraphicsItem = new CGAL::Qt::ArrangementGraphicsItem< Arrangement >( this->arrangement );
        this->segmentInputCallback = new ArrangementSegmentInputCallback< Arrangement >( this->arrangement, this );
        this->deleteCurveCallback = new DeleteCurveCallback< Arrangement >( this->arrangement, this );
        this->pointLocationCallback = new PointLocationCallback< Arrangement >( this->arrangement, this );
        this->verticalRayShootCallback = new VerticalRayShootCallback< Arrangement >( this->arrangement, this );
        this->mergeEdgeCallback = new MergeEdgeCallback< Arrangement >( this->arrangement, this );
        this->splitEdgeCallback = new SplitEdgeCallback< Arrangement >( this->arrangement, this );
        this->envelopeCallback = new EnvelopeCallback< Arrangement >( this->arrangement, this );

        this->scene->addItem( this->arrangementGraphicsItem );

        this->segmentInputCallback->setScene( this->scene );
        this->deleteCurveCallback->setScene( this->scene );
        this->pointLocationCallback->setScene( this->scene );
        this->verticalRayShootCallback->setScene( this->scene );
        this->mergeEdgeCallback->setScene( this->scene );
        this->splitEdgeCallback->setScene( this->scene );
        this->envelopeCallback->setScene( this->scene );
    }

protected:
    Arrangement* arrangement;

}; // class ArrangementDemoTab

#endif // ARRANGEMENT_DEMO_TAB_H
