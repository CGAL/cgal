#ifndef DELETE_CURVE_MODE_ITEM_EDITOR_H
#define DELETE_CURVE_MODE_ITEM_EDITOR_H

#include <QComboBox>
#include "DeleteCurveMode.h"

class QWidget;

class DeleteCurveModeItemEditor : public QComboBox
{
  Q_OBJECT
  Q_PROPERTY( DeleteCurveMode mode READ mode WRITE setMode USER true )

public:
  DeleteCurveModeItemEditor( QWidget* parent = 0 );

public:
  DeleteCurveMode mode( ) const;
  void setMode( DeleteCurveMode m );

}; // class DeleteCurveModeItemEditor

#endif // DELETE_CURVE_MODE_ITEM_EDITOR_H
