#include "MainWindow.h"
#include <string>


MainWindow::MainWindow(QWidget* parent)
 : CGAL::Qt::DemosMainWindow(parent)
{
  //  Qt Automatic Connections
  //   https://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
  //  setupUi(this) automatically generates connections to the slots named
  //   "on_<action_name>_<signal_name>"
  setupUi(this);

  this->viewer->setScene(&m_scene);
  m_scene.setViewer(this->viewer);

  // setup Mode menu group actions
  QActionGroup *modeGroup = new QActionGroup(this);
  modeGroup->addAction(this->actionNormal_View);
  modeGroup->addAction(this->actionInsert_Vertex);
  modeGroup->addAction(this->actionInsert_Point);
  modeGroup->addAction(this->actionSelect_Vertex);
  modeGroup->addAction(this->actionMove_Vertex);
  modeGroup->addAction(this->actionFind_NearestNb);
  modeGroup->addAction(this->actionEmpty_Sphere);
  QObject::connect(modeGroup, SIGNAL(triggered(QAction *)), this, SLOT(setMode(QAction *)));

  // connect menu actions to viewer slots
  connectActions();

  // About menu
  // addAboutCGAL() is a function in DemoMainWindow
  //   it will add a menu action "About CGAL..." to Help menu and connect to popupAboutCGAL
  //   default popupAboutCGAL points to a fixed file directory ":/cgal/help/about_CGAL.html"
  //   here we override it with our directory
  this->addAboutCGAL();
  // addAboutDemo(QString htmlResourceName) is also a function in DemoMainWindow
  //   it will add a menu action "About Demo..." to Help menu
  //   when the action is invoked, it will popup a messageBox showing the given html
  this->addAboutDemo( ":/documentation/documentation/about.html" );

}

void MainWindow::connectActions()
{
  // Edit menu actions
  QObject::connect(this->actionIncremental_Construct, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleIncremental(bool)));
  QObject::connect(this->actionStop_Animation, SIGNAL(triggered()),
           this->viewer, SLOT(stopIncremental()));
  QObject::connect(this->viewer, SIGNAL(stopIncAnimation()),
           this, SLOT(stopAnimation()));

  // Show menu actions
  QObject::connect(this->actionShow_Axis, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleShowAxis(bool)));
  QObject::connect(this->actionShow_Vertex, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleShowVertex(bool)));
  QObject::connect(this->actionShow_DEdge, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleShowDEdge(bool)));
  QObject::connect(this->actionShow_VEdge, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleShowVEdge(bool)));
  QObject::connect(this->actionShow_Facet, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleShowFacet(bool)));
  QObject::connect(this->actionFlat, SIGNAL(toggled(bool)),
           this->viewer, SLOT(toggleFlat(bool)));

  // Preferences
  QObject::connect(this->actionPreferences, SIGNAL(triggered()),
           this->viewer, SLOT(setPreferences()));

  // Help menu actions
  QObject::connect(this->actionDemo_Help, SIGNAL(triggered()),
           this->viewer, SLOT(help()));
  QObject::connect(this->actionAbout_T3_demo, SIGNAL(triggered()),
       this, SLOT(popupAboutDemo()));

  // Quit
  QObject::connect(this->actionQuit, SIGNAL(triggered()),
           qApp, SLOT(closeAllWindows()));

  // Viewer signals
  QObject::connect(this, SIGNAL(sceneChanged()),
           this->viewer, SLOT(update()));
}

void MainWindow::closeEvent(QCloseEvent * /*event*/)
{
  viewer->writeSettings();
}

/*************************************************************/
/*  Action functions */

void MainWindow::setMode(QAction *action)
{
  if( action == this->actionNormal_View )
    this->viewer->setMode( viewer->NONE );
  else if( action == this->actionInsert_Vertex )
    this->viewer->setMode( viewer->INSERT_V );
  else if( action == this->actionInsert_Point )
    this->viewer->setMode( viewer->INSERT_PT );
  else if( action == this->actionSelect_Vertex )
    this->viewer->setMode( viewer->SELECT );
  else if( action == this->actionMove_Vertex )
    this->viewer->setMode( viewer->MOVE );
  else if( action == this->actionFind_NearestNb )
    this->viewer->setMode( viewer->FINDNB );
  else if( action == this->actionEmpty_Sphere )
    this->viewer->setMode( viewer->EMPTYSPH );
}

void MainWindow::on_actionLoad_Points_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                     tr("Open an file"),        // dialog caption
                     ".",        // initial directory
                     tr("OFF files (*.off);;XYZ files (*.xyz);;All files (*.*)"));        // selection filter
  if( fileName.isEmpty() )  return;

  // erase old data
  viewer->clear();

  // parse fileName to get the file type
  std::string fname = fileName.toLatin1().data();
  std::string ftype = fname.substr( fname.find_last_of('.')+1 );

  if ( ftype.compare("off")==0 || ftype.compare("OFF")==0 ) { // read from OFF file
    m_scene.loadPointsOFF( fname.data() );
    // set selectBuffer size (if necessary)
    viewer->setSelBuffSize();
  } else if ( ftype.compare("xyz")==0 || ftype.compare("XYZ")==0 ) { // read from XYZ file
    m_scene.loadPointsXYZ( fname.data() );
    // set selectBuffer size (if necessary)
    viewer->setSelBuffSize();
  } else {
    viewer->displayMessage( tr("Please select an OFF or XYZ file to open!") );
  }

  // update viewer
  Q_EMIT( sceneChanged() );
}

void MainWindow::on_actionSave_Points_triggered()
{
  if( m_scene.isDTEmpty() ) {
    viewer->displayMessage( tr("Error: no existing triangulation to be saved.") );
    return;
  }

  QString fileName = QFileDialog::getSaveFileName(this,
                     tr("Save an file"),        // dialog caption
                     ".",        // initial directory
                     tr("OFF files (*.off);;XYZ files (*.xyz);;All files (*.*)"));        // selection filter
  if( fileName.isEmpty() )  return;

  // parse fileName to get the file type
  std::string fname = fileName.toLatin1().data();//toAscii()
  std::string ftype = fname.substr( fname.find_last_of('.')+1 );

  if ( ftype.compare("off")==0 || ftype.compare("OFF")==0 ) { // save to OFF file
    m_scene.savePointsOFF( fname.data() );
  } else if ( ftype.compare("xyz")==0 || ftype.compare("XYZ")==0 ) { // save to XYZ file
    m_scene.savePointsXYZ( fname.data() );
  } else {
    viewer->displayMessage( tr("Please select an OFF or XYZ file to open!") );
  }
}

void MainWindow::on_actionGenerate_Points_triggered()
{
  bool isOk;

  int nPoints = QInputDialog::getInt(this,
              "3D Triangulation demo", "Number of points: ",        // caption and label
              100, // default value
              4, // min value
              2147483647, // max value
              1, // step value of arrow button
              &isOk); // if OK is pressed


  if ( isOk) {
    // erase old data
    viewer->clear();

          // generate points
    m_scene.generatePoints(nPoints);
    // set selectBuffer size (if necessary)
    viewer->setSelBuffSize();

    // update viewer
    Q_EMIT( sceneChanged() );
    viewer->changed();


  }// if(isOk)
}

void MainWindow::stopAnimation()
{
  if( this->actionIncremental_Construct->isChecked() )
    this->actionIncremental_Construct->setChecked( false );
}

void MainWindow::on_actionClear_Scene_triggered()
{
  viewer->clear();

  // update viewer
  Q_EMIT( sceneChanged() );
}

void MainWindow::popupAboutCGAL()
{
  // read contents from .html file
  QFile about_CGAL(":/documentation/documentation/about.html");
  about_CGAL.open(QIODevice::ReadOnly|QIODevice::Text);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();

  // popup a message box
  QMessageBox mb(QMessageBox::NoIcon,
                 tr("About CGAL..."),
                 about_CGAL_txt,
                 QMessageBox::Ok,
                 this);

  // set links to be accessible by mouse or keyboard
  QLabel* mb_label = mb.findChild<QLabel*>("qt_msgbox_label");
  if(mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() |
                                      ::Qt::LinksAccessibleByMouse |
                                      ::Qt::LinksAccessibleByKeyboard);
  } else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }

  mb.exec();
}
