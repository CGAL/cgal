
#include "MainWindow.h"

MainWindow::MainWindow(QWidget* parent): CGAL::Qt::DemosMainWindow(parent)
{
  setupUi(this);
  this->viewer->setScene(&scene);
  connectActions();
  this->addAboutDemo(":/cgal/help/about_Alpha_shapes_3.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::connectActions()
{
  QObject::connect(this->actionLoad_New_File, SIGNAL(triggered()), 
		   this, SLOT(open_file()));

  QObject::connect(this->alphaSlider, SIGNAL(valueChanged(int)), 
		   this, SLOT(alphaChanged(int)));

  QObject::connect(this->alphaBox, SIGNAL(valueChanged(int)),
		   this, SLOT(alphaChanged(int)));

  QObject::connect(this->alphaSlider, SIGNAL(valueChanged(int)), 
		   this->alphaBox, SLOT(setValue(int)));

  QObject::connect(this->alphaBox, SIGNAL(valueChanged(int)), 
		   this->alphaSlider, SLOT(setValue(int)));

  QObject::connect(this, SIGNAL(sceneChanged()), 
		   this->viewer, SLOT(sceneChanged()));

  QObject::connect(this, SIGNAL(alphaChanged()), 
		   this->viewer, SLOT(update()));


  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   qApp, SLOT(quit()));
}

void
MainWindow::open_file()
{

  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points File"),
						  "./data",
						  tr("pts files (*.pts)"));

  if(! fileName.isEmpty()){
    open(fileName);
  }
}


void 
MainWindow::alphaChanged(int i)
{
  if (scene.alpha_shape.number_of_alphas() > 0){
    if(i < 100){
      int n = static_cast<int>((i * scene.alpha_shape.number_of_alphas())/ 100);
      if(n == 0) n++;
      scene.alpha_shape.set_alpha(scene.alpha_shape.get_nth_alpha(n));
    } else {
      Alpha_iterator alpha_end_it = scene.alpha_shape.alpha_end();
      scene.alpha_shape.set_alpha((*(--alpha_end_it))+1);
    }
  } else {
    scene.alpha_shape.set_alpha(0);
  }
  emit (alphaChanged());
}

void
MainWindow::open(QString fileName)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  scene.alpha_shape.clear();
  scene.points.clear();
  std::ifstream ifs(qPrintable(fileName));
  int n;
  ifs >> n;
  Point_3 p;
  for(int i=0; i<n; i++){
    ifs >> p;
    scene.points.push_back(p);
  }
  timer.reset();
  timer.start();
  scene.alpha_shape.make_alpha_shape(scene.points.begin(), scene.points.end());
  scene.alpha_shape.set_alpha(16);
  timer.stop();
  

  alphaSlider->setRange(0,100);
  alphaSlider->setSliderPosition(50);

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor();
  emit (sceneChanged());
}


#include "MainWindow.moc"

