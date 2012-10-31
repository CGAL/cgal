#include "DeleteCurveModeItemEditor.h"

DeleteCurveModeItemEditor::DeleteCurveModeItemEditor( QWidget* parent ):
  QComboBox( parent )
{
  this->setFrame( false );

  QVariant deleteCurveOption = QVariant::fromValue( DeleteCurveMode( ) );
  QVariant deleteEdgeOption =
    QVariant::fromValue( DeleteCurveMode( DeleteCurveMode::DELETE_EDGE ) );
  this->insertItem( 0, "Delete Curve", deleteCurveOption );
  this->insertItem( 1, "Delete Edge", deleteEdgeOption );
}

DeleteCurveMode DeleteCurveModeItemEditor::mode( ) const
{
  return qVariantValue<DeleteCurveMode >(this->itemData(this->currentIndex( ),
                                                        Qt::UserRole ) );
}

void
DeleteCurveModeItemEditor::
setMode( DeleteCurveMode m )
{
  if ( m.mode( ) == DeleteCurveMode::DELETE_CURVE )
  {
    this->setCurrentIndex( 0 );
  }
  else
  {
    this->setCurrentIndex( 1 );
  }
}
