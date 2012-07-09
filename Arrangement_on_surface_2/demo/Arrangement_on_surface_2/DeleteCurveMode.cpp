#include "DeleteCurveMode.h"

#include <QString>

DeleteCurveMode::
DeleteCurveMode( ):
    m_mode( DELETE_CURVE )
{ }
DeleteCurveMode::
DeleteCurveMode( const DeleteCurveMode& dcm ):
    m_mode( dcm.mode( ) )
{
}

DeleteCurveMode::
DeleteCurveMode( Mode mode ):
    m_mode( mode )
{ }

DeleteCurveMode::
~DeleteCurveMode( )
{ }

DeleteCurveMode::Mode 
DeleteCurveMode::
mode( ) const
{
    return this->m_mode;
}

void 
DeleteCurveMode::
setMode( Mode mode )
{
    this->m_mode = mode;
}

QString
DeleteCurveMode::
ToString( const DeleteCurveMode& mode )
{
    if ( mode.mode( ) == DELETE_CURVE )
    {
        return QString("Delete Curve");
    }
    else
    {
        return QString("Delete Edge");
    }
}
