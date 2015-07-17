#ifndef PREFERENCE_DLG_H
#define PREFERENCE_DLG_H

#include <QDialog>
#include <QColorDialog>

#include <QGridLayout>
#include <QGroupBox> 

#include <QLabel>
#include <QLineEdit>

#include <QMessageBox>

#include <QPushButton>

#include <QSpinBox>


class QLabel;
class QSpinBox;
class QLineEdit;

class PreferenceDlg : public QDialog
{
  Q_OBJECT

  friend class Viewer;

public:
  PreferenceDlg(QWidget *parent=0);

private:
  void init(QColor, float, QColor, float, QColor, float, QColor, QColor, QColor, int);

private Q_SLOTS:
  void okClicked() { hide(); Q_EMIT( applyChanges() ); }
  void applyClicked() { Q_EMIT( applyChanges() ); }

  void setVertexColor();
  void setVertexSize(const QString&);
  void setDEdgeColor();
  void setDEdgeSize(const QString&);
  void setVEdgeColor();
  void setVEdgeSize(const QString&);
  void setFacetColor();
  void setFacetAlpha();
  void setTrackballColor();
  void setTrackballAlpha();
  void setStepLong();
  void setEmptySphereColor();
  void setEmptySphereAlpha();

  Q_SIGNALS: // Signals do not have access specifier
  void applyChanges();

private:
  QLabel *m_labelVertex;
  QLineEdit *m_editSizeV;
  QLabel *m_labelDEdge;
  QLineEdit *m_editSizeDE;
  QLabel *m_labelVEdge;
  QLineEdit *m_editSizeVE;
  QLabel *m_labelFacet;
  QSpinBox *m_spinAlphaF;
  QLabel *m_labelBall;
  QSpinBox *m_spinAlphaB;
  QSpinBox *m_spinStep;
  QLabel *m_labelSphere;
  QSpinBox *m_spinAlphaS;

  float m_fSizeVertex;
  float m_fSizeDEdge;
  float m_fSizeVEdge;
  QColor m_colorVertex;
  QColor m_colorDEdge;
  QColor m_colorVEdge;
  QColor m_colorFacet;
  QColor m_colorTrackball;
  int m_iStep;
  QColor m_colorEmptySphere;
};

#endif
