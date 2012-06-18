#include "VerticalRayShootCallback.h"

VerticalRayShootCallbackBase::
VerticalRayShootCallbackBase( QObject* parent_ ):
    CGAL::Qt::Callback( parent_ ),
    shootingUp( true )
{ }

void 
VerticalRayShootCallbackBase::
setShootingUp( bool isShootingUp )
{
    this->shootingUp = isShootingUp;
}

