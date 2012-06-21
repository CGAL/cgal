#include "SplitEdgeCallback.h"
void
SplitEdgeCallbackBase::
setSnappingEnabled( bool b )
{
    this->snappingEnabled = b;
}

void
SplitEdgeCallbackBase::
setSnapToGridEnabled( bool b )
{
    this->snapToGridEnabled = b;
}

SplitEdgeCallbackBase::
SplitEdgeCallbackBase( QObject* parent ):
    CGAL::Qt::Callback( parent ),
    snappingEnabled( false ),
    snapToGridEnabled( false )
{ }
