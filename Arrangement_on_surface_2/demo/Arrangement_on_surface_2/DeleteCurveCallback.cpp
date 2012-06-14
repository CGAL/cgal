#include "DeleteCurveCallback.hpp"

DeleteCurveCallback::
DeleteCurveCallback( QObject* parent ):
    QObject( parent )
{ }

void
DeleteCurveCallback::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

QGraphicsScene*
DeleteCurveCallback::
getScene( ) const
{
    return this->scene;
}

