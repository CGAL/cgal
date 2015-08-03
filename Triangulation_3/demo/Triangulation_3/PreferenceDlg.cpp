#include <QtGui>
#include "PreferenceDlg.h"


PreferenceDlg::PreferenceDlg(QWidget *parent) : QDialog(parent)
{
  /* Vertex */

  // create groupbox
  QGroupBox *groupV = new QGroupBox( tr("Vertex") );
  // create buttons
  QPushButton *btnVertex = new QPushButton( tr("Set Color") );
  // create color label
  m_labelVertex = new QLabel;
  m_labelVertex->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create size label
  QLabel *labelSizeV = new QLabel( tr("Set Size") );
  // create lineedit
  m_editSizeV = new QLineEdit;

  // connect to actions
  connect( btnVertex, SIGNAL(clicked()), this, SLOT(setVertexColor()) );
  connect( m_editSizeV, SIGNAL(textChanged(const QString&)), this, SLOT(setVertexSize(const QString&)) );

  // lay out the buttons
  QGridLayout *layoutV = new QGridLayout;
  layoutV->addWidget( btnVertex, 0, 0 );
  layoutV->addWidget( m_labelVertex, 0, 1 );
  layoutV->addWidget( labelSizeV, 1, 0 );
  layoutV->addWidget( m_editSizeV, 1, 1 );
  groupV->setLayout( layoutV );

  /* Delaunau Edge */

  // create groupbox
  QGroupBox *groupDE = new QGroupBox( tr("Delaunay Edge") );
  // create button
  QPushButton *btnDEdge = new QPushButton( tr("Set Color") );
  // create color label
  m_labelDEdge = new QLabel;
  m_labelDEdge->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create size label
  QLabel *labelSizeDE = new QLabel( tr("Set Size") );
  // create lineedit
  m_editSizeDE = new QLineEdit;

  // connect to actions
  connect( btnDEdge, SIGNAL(clicked()), this, SLOT(setDEdgeColor()) );
  connect( m_editSizeDE, SIGNAL(textChanged(const QString&)), this, SLOT(setDEdgeSize(const QString&)) );

  // lay out the buttons
  QGridLayout *layoutDE = new QGridLayout;
  layoutDE->addWidget( btnDEdge, 0, 0 );
  layoutDE->addWidget( m_labelDEdge, 0, 1 );
  layoutDE->addWidget( labelSizeDE, 1, 0 );
  layoutDE->addWidget( m_editSizeDE, 1, 1 );
  groupDE->setLayout( layoutDE );

  /* Voronoi Edge */

  // create groupbox
  QGroupBox *groupVE = new QGroupBox( tr("Voronoi Edge") );
  // create button
  QPushButton *btnVEdge = new QPushButton( tr("Set Color") );
  // create color label
  m_labelVEdge = new QLabel;
  m_labelVEdge->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create size label
  QLabel *labelSizeVE = new QLabel( tr("Set Size") );
  // create lineedit
  m_editSizeVE = new QLineEdit;

  // connect to actions
  connect( btnVEdge, SIGNAL(clicked()), this, SLOT(setVEdgeColor()) );
  connect( m_editSizeVE, SIGNAL(textChanged(const QString&)), this, SLOT(setVEdgeSize(const QString&)) );

  // lay out the buttons
  QGridLayout *layoutVE = new QGridLayout;
  layoutVE->addWidget( btnVEdge, 0, 0 );
  layoutVE->addWidget( m_labelVEdge, 0, 1 );
  layoutVE->addWidget( labelSizeVE, 1, 0 );
  layoutVE->addWidget( m_editSizeVE, 1, 1 );
  groupVE->setLayout( layoutVE );

  /* Facet */

  // create groupbox
  QGroupBox *groupF = new QGroupBox( tr("Facet") );
  // create button
  QPushButton *btnFacet = new QPushButton( tr("Set Color") );
  // create color label
  m_labelFacet = new QLabel;
  m_labelFacet->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create label and spinbox
  QLabel *labelFacetA = new QLabel( tr("Transparency") );
  m_spinAlphaF = new QSpinBox;
  m_spinAlphaF->setRange(0, 255);

  // connect to actions
  connect( btnFacet, SIGNAL(clicked()), this, SLOT(setFacetColor()) );
  connect( m_spinAlphaF, SIGNAL(valueChanged(int)), this, SLOT(setFacetAlpha()) );

  // lay out the buttons
  QGridLayout *layoutF = new QGridLayout;
  layoutF->addWidget( btnFacet, 0, 0 );
  layoutF->addWidget( m_labelFacet, 0, 1 );
  layoutF->addWidget( labelFacetA, 1, 0 );
  layoutF->addWidget( m_spinAlphaF, 1, 1 );
  groupF->setLayout( layoutF );

  /* Trackball  */

  // create groupbox
  QGroupBox *groupB = new QGroupBox( tr("Trackball") );
  // create button
  QPushButton *btnBall = new QPushButton( tr("Set Color") );
  // create color label
  m_labelBall = new QLabel;
  m_labelBall->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create label and spinbox
  QLabel *labelBallA = new QLabel( tr("Transparency") );
  m_spinAlphaB = new QSpinBox;
  m_spinAlphaB->setRange(0, 255);
  // create label and spinbox
  QLabel *labelStep = new QLabel( tr("Step-long of Resizing") );
  m_spinStep = new QSpinBox;
  m_spinStep->setRange(1, 300);

  // connect to actions
  connect( btnBall, SIGNAL(clicked()), this, SLOT(setTrackballColor()) );
  connect( m_spinAlphaB, SIGNAL(valueChanged(int)), this, SLOT(setTrackballAlpha()) );
  connect( m_spinStep, SIGNAL(valueChanged(int)), this, SLOT(setStepLong()) );

  // lay out the buttons
  QGridLayout *layoutB = new QGridLayout;
  layoutB->addWidget( btnBall, 0, 0 );
  layoutB->addWidget( m_labelBall, 0, 1 );
  layoutB->addWidget( labelBallA, 1, 0 );
  layoutB->addWidget( m_spinAlphaB, 1, 1 );
  layoutB->addWidget( labelStep, 2, 0 );
  layoutB->addWidget( m_spinStep, 2, 1 );
  groupB->setLayout( layoutB );

  /* Empty Sphere  */

  // create groupbox
  QGroupBox *groupS = new QGroupBox( tr("Empty Sphere") );
  // create color label
  m_labelSphere = new QLabel;
  m_labelSphere->setFrameStyle(QFrame::Sunken | QFrame::Panel);
  // create button
  QPushButton *btnSphere = new QPushButton( tr("Set Color") );
  // create label and spinbox
  QLabel *labelSphereA = new QLabel( tr("Transparency") );
  m_spinAlphaS = new QSpinBox;
  m_spinAlphaS->setRange(0, 255);

  // connect to actions
  connect( btnSphere, SIGNAL(clicked()), this, SLOT(setEmptySphereColor()) );
  connect( m_spinAlphaS, SIGNAL(valueChanged(int)), this, SLOT(setEmptySphereAlpha()) );

  // lay out the buttons
  QGridLayout *layoutS = new QGridLayout;
  layoutS->addWidget( btnSphere, 0, 0 );
  layoutS->addWidget( m_labelSphere, 0, 1 );
  layoutS->addWidget( labelSphereA, 1, 0 );
  layoutS->addWidget( m_spinAlphaS, 1, 1 );
  groupS->setLayout( layoutS );

  /* OK buttons */
  // create groupbox
  QGroupBox *groupBtn = new QGroupBox();
  // buttons
  QPushButton *ok = new QPushButton( tr("OK") );
  QPushButton *apply = new QPushButton( tr("Apply") );
  QPushButton *cancel = new QPushButton( tr("Cancel") );
  cancel->setFocus();

  // connect to actions
  connect( ok, SIGNAL(clicked()), this, SLOT(okClicked()) );
  connect( apply, SIGNAL(clicked()), this, SLOT(applyClicked()) );
  connect( cancel, SIGNAL(clicked()), this, SLOT(reject()) );

  // lay out the buttons
  QGridLayout *layoutBtn = new QGridLayout;
  layoutBtn->addWidget( ok, 0, 0 );
  layoutBtn->addWidget( cancel, 0, 1 );
  layoutBtn->addWidget( apply, 0, 2 );
  groupBtn->setLayout( layoutBtn );

  /* dialog layout */

  // lay out the buttons
  QGridLayout *main = new QGridLayout;
  main->addWidget( groupV, 0, 1 );
  main->addWidget( groupDE, 0, 2 );
  main->addWidget( groupVE, 0, 3 );
  main->addWidget( groupF, 1, 1 );
  main->addWidget( groupB, 1, 2 );
  main->addWidget( groupS, 1, 3 );
  main->addWidget( groupBtn, 2, 2, 2, 3 );
  setLayout( main );

  // set dialog title
  setWindowTitle( tr("Preferences") );
}

void PreferenceDlg::init(QColor clrVt, float sizeV, QColor clrDE, float sizeDE,
                       QColor clrVE, float sizeVE,
                       QColor clrF, QColor clrB, QColor clrS, int iStep)
{
  // vertex color
  m_colorVertex = clrVt;
  // show the color in label
  m_labelVertex->setText(m_colorVertex.name());
  m_labelVertex->setPalette( QPalette(m_colorVertex) );
  m_labelVertex->setAutoFillBackground(true);
  // vertex size
  m_fSizeVertex = sizeV;
  m_editSizeV->setText( QString::number( m_fSizeVertex ) );

  // Delaunay edge color
  m_colorDEdge = clrDE;
  // show the color in label
  m_labelDEdge->setText( m_colorDEdge.name() );
  m_labelDEdge->setPalette( QPalette(m_colorDEdge) );
  m_labelDEdge->setAutoFillBackground(true);
  // edge size
  m_fSizeDEdge = sizeDE;
  m_editSizeDE->setText( QString::number( m_fSizeDEdge ) );

  // Voronoi edge color
  m_colorVEdge = clrVE;
  // show the color in label
  m_labelVEdge->setText( m_colorVEdge.name() );
  m_labelVEdge->setPalette( QPalette(m_colorVEdge) );
  m_labelVEdge->setAutoFillBackground(true);
  // edge size
  m_fSizeVEdge = sizeVE;
  m_editSizeVE->setText( QString::number( m_fSizeVEdge ) );

  // facet color
  m_colorFacet = clrF;
  // show the color in label
  m_labelFacet->setText( m_colorFacet.name() );
  m_labelFacet->setPalette( QPalette(m_colorFacet) );
  m_labelFacet->setAutoFillBackground(true);
  // facet transparency
  m_spinAlphaF->setValue( m_colorFacet.alpha() );

  // trackball color
  m_colorTrackball = clrB;
  // show the color in label
  m_labelBall->setText(m_colorTrackball.name());
  m_labelBall->setPalette( QPalette(m_colorTrackball) );
  m_labelBall->setAutoFillBackground(true);
  // trackball transparency
  m_spinAlphaB->setValue( m_colorTrackball.alpha() );
  // trackball resizing fineness
  m_spinStep->setValue( iStep );

  // empty sphere color
  m_colorEmptySphere = clrS;
  // show the color in label
  m_labelSphere->setText(m_colorEmptySphere.name());
  m_labelSphere->setPalette( QPalette(m_colorEmptySphere) );
  m_labelSphere->setAutoFillBackground(true);
  // trackball transparency
  m_spinAlphaS->setValue( m_colorEmptySphere.alpha() );
}

void PreferenceDlg::setVertexColor()
{
  m_colorVertex = QColorDialog::getColor(m_colorVertex, this);
  if( m_colorVertex.isValid() ) {
    m_labelVertex->setText(m_colorVertex.name());
    m_labelVertex->setPalette( QPalette(m_colorVertex) );
    m_labelVertex->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setVertexSize(const QString& str)
{
  bool ok;
  float size = str.toFloat(&ok);
  if( ok )
    m_fSizeVertex = size;
  else {
      QMessageBox mb(QMessageBox::NoIcon, tr("Error!"),
                 tr("Enter a valid floating number."),
                 QMessageBox::Ok, this);
      mb.exec();
      m_editSizeV->setFocus();
  }
}

void PreferenceDlg::setDEdgeColor()
{
  m_colorDEdge = QColorDialog::getColor(m_colorDEdge, this);
  if( m_colorDEdge.isValid() ) {
    m_labelDEdge->setText( m_colorDEdge.name() );
    m_labelDEdge->setPalette( QPalette(m_colorDEdge) );
    m_labelDEdge->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setDEdgeSize(const QString& str)
{
  bool ok;
  float size = str.toFloat(&ok);
  if( ok )
    m_fSizeDEdge = size;
  else {
      QMessageBox mb(QMessageBox::NoIcon, tr("Error!"),
                 tr("Enter a valid floating number."),
                 QMessageBox::Ok, this);
      mb.exec();
      m_editSizeDE->setFocus();
  }
}

void PreferenceDlg::setVEdgeColor()
{
  m_colorVEdge = QColorDialog::getColor(m_colorVEdge, this);
  if( m_colorVEdge.isValid() ) {
    m_labelVEdge->setText( m_colorVEdge.name() );
    m_labelVEdge->setPalette( QPalette(m_colorVEdge) );
    m_labelVEdge->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setVEdgeSize(const QString& str)
{
  bool ok;
  float size = str.toFloat(&ok);
  if( ok )
    m_fSizeVEdge = size;
  else {
      QMessageBox mb(QMessageBox::NoIcon, tr("Error!"),
                 tr("Enter a valid floating number."),
                 QMessageBox::Ok, this);
      mb.exec();
      m_editSizeVE->setFocus();
  }
}

void PreferenceDlg::setFacetColor()
{
  m_colorFacet = QColorDialog::getColor(m_colorFacet, this);
  if( m_colorFacet.isValid() ) {
    m_labelFacet->setText( m_colorFacet.name() );
    m_colorFacet.setAlpha( m_spinAlphaF->value() );
    m_labelFacet->setPalette( QPalette(m_colorFacet) );
    m_labelFacet->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setFacetAlpha()
{
  m_colorFacet.setAlpha( m_spinAlphaF->value() );
}

void PreferenceDlg::setTrackballColor()
{
  m_colorTrackball = QColorDialog::getColor(m_colorTrackball, this);
  if( m_colorTrackball.isValid() ) {
    m_labelBall->setText( m_colorTrackball.name() );
    m_colorTrackball.setAlpha( m_spinAlphaB->value() );
    m_labelBall->setPalette( QPalette(m_colorTrackball) );
    m_labelBall->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setTrackballAlpha()
{
  m_colorTrackball.setAlpha( m_spinAlphaB->value() );
}

void PreferenceDlg::setStepLong()
{
  m_iStep = m_spinStep->value();
}

void PreferenceDlg::setEmptySphereColor()
{
  m_colorEmptySphere = QColorDialog::getColor(m_colorEmptySphere, this);
  if( m_colorEmptySphere.isValid() ) {
    m_labelSphere->setText( m_colorEmptySphere.name() );
    m_colorEmptySphere.setAlpha( m_spinAlphaS->value() );
    m_labelSphere->setPalette( QPalette(m_colorEmptySphere) );
    m_labelSphere->setAutoFillBackground(true);
  }
}

void PreferenceDlg::setEmptySphereAlpha()
{
  m_colorEmptySphere.setAlpha( m_spinAlphaS->value() );
}
