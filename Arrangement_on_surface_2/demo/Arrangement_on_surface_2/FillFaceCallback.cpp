#include "FillFaceCallback.h"
FillFaceCallbackBase::
FillFaceCallbackBase( QObject* parent ):
    CGAL::Qt::Callback( parent ),
    fillColor( ::Qt::black )
{

}

void
FillFaceCallbackBase::
setColor( QColor c )
{
    this->fillColor = c;
    emit modelChanged( );
}

QColor
FillFaceCallbackBase::
getColor( ) const
{
    return this->fillColor;
}
