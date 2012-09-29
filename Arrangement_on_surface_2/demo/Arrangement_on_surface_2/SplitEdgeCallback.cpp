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

void 
SplitEdgeCallbackBase::
setColor( QColor c )
{
    this->color = c;
}

QColor 
SplitEdgeCallbackBase::
getColor( ) const
{
    return this->color;
}

SplitEdgeCallbackBase::
SplitEdgeCallbackBase( QObject* parent ):
    CGAL::Qt::Callback( parent ),
    snappingEnabled( false ),
    snapToGridEnabled( false ),
    color( ::Qt::blue )
{ }
