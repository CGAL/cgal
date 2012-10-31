#ifndef DELETE_CURVE_MODE_H
#define DELETE_CURVE_MODE_H

#include <QMetaType>

class QString;

/**
An attribute describing the policy for deleting curves from the arrangement.
*/
class DeleteCurveMode
{
public:
  enum Mode {
    DELETE_CURVE,
    DELETE_EDGE
  };

  DeleteCurveMode( );
  DeleteCurveMode( const DeleteCurveMode& dcm );
  DeleteCurveMode( Mode mode );
  ~DeleteCurveMode( );

  Mode mode( ) const;
  void setMode( Mode mode );

  static QString ToString( const DeleteCurveMode& mode );

protected:
  Mode m_mode;
};

Q_DECLARE_METATYPE( DeleteCurveMode )

#endif // DELETE_CURVE_MODE_H
