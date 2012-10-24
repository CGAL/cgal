#include <fstream>

// CGAL headers
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_hyperbolic_triangulation_2.h>
#include <CGAL/Triangulation_hyperbolic_traits_2.h>

// to be deleted
#include <CGAL/Qt/HyperbolicPainterOstream.h>
//

#include <CGAL/point_generators_2.h>

// Maintain translations
#include <TranslationInfo.h>
#include <Translations.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

// experiment
#include "GroupOfIndex2.h"
#include "OriginalDomainNeighbors.h"
#include "OriginalDomainNeighborsCommon.h"
//

// GraphicsView items and event filters (input classes)
#include "PointTranslationWithInfo.h"
#include "TriangulationCircumcircle.h"
#include "TriangulationMovingPoint.h"
#include "TriangulationConflictZone.h"
#include "TriangulationRemoveVertex.h"
#include "TriangulationPointInputAndConflictZone.h"
#include <CGAL/Qt/TriangulationGraphicsItemWithColorInfo.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_Hyperbolic_translations_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel R;
typedef CGAL::Triangulation_hyperbolic_traits_2<R> K;

typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K Gt;

// keep the name of the translation for every vertex
typedef TranslationInfo<std::wstring> Vb_info;

typedef CGAL::Triangulation_vertex_base_with_info_2< Vb_info, Gt > Vb;
typedef CGAL::Triangulation_face_base_with_info_2 <CGAL::Hyperbolic_face_info_2, Gt > Fb;

typedef CGAL::Delaunay_hyperbolic_triangulation_2< Gt, CGAL::Triangulation_data_structure_2<Vb, Fb> > Delaunay;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  
  Delaunay dt;
  QGraphicsEllipseItem* disk;
  QGraphicsScene scene;  

  CGAL::Qt::TriangulationGraphicsItem<Delaunay> * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Delaunay> * vgi;

  CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::TriangulationConflictZone<Delaunay> * cz;
  CGAL::Qt::TriangulationRemoveVertex<Delaunay> * trv;
  CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
  
  // hyperbolic translations
  typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
  typedef Hyperbolic_isometry::complex complex;
  
  // storage of static data: translations a, b, c, d
  typedef Translations<Gt> Translations;
  typedef CGAL::Qt::PointTranslationWithInfo<Delaunay, Hyperbolic_isometry> PointTranslation;
  
  PointTranslation *trs_a;
  PointTranslation *trs_b;
  PointTranslation *trs_c;
  PointTranslation *trs_d;
  
  // experiment
  typedef CGAL::Qt::OriginalDomainNeighbors<Delaunay, Hyperbolic_isometry> OriginalDomainNeighbors;
  typedef Group_of_index_2<Hyperbolic_isometry> Group_of_index_2;
  
  
  // Oct 2012
  typedef IsometryWithInfo<Gt, TranslationInfo<std::wstring> > TranslationWithInfo;
  typedef CGAL::Qt::OriginalDomainNeighborsCommon<Delaunay, TranslationWithInfo> OriginalDomainNeighborsWithInfo;
  
  std::vector<TranslationWithInfo> cloud_of_translations;
  
  OriginalDomainNeighborsWithInfo *odnwi;
  OriginalDomainNeighborsWithInfo *odnwi2;
  //
  
  // delete later, for some experiments
  Group_of_index_2::List words;
  //
  
  Group_of_index_2 g2;
  Group_of_index_2 g4;
  Group_of_index_2 g8;
  Group_of_index_2 g16;
  // extra groups, we can delete them
  Group_of_index_2 g32;
  Group_of_index_2 g64;
  Group_of_index_2 g128;
  Group_of_index_2 g256;
  Group_of_index_2 g512;
  Group_of_index_2 g1024;
  //
  OriginalDomainNeighbors *odn;
  OriginalDomainNeighbors *odn2;
  //
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionDo_translation_a_toggled(bool checked);

  void on_actionDo_translation_b_toggled(bool checked);

  void on_actionDo_translation_c_toggled(bool checked);

  void on_actionDo_translation_d_toggled(bool checked);
  
  void on_actionG_toggled(bool checked);
  
  void on_actionG2_toggled(bool checked);
  
  void on_actionG4_toggled(bool checked);
  
  void on_actionG8_toggled(bool checked);
  
  void on_actionG16_toggled(bool checked);

  void on_actionShowConflictZone_toggled(bool checked);

  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionSavePoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

signals:
  void changed();
};


// delete
bool inverses(int i, int j)
{
  if((i == 0 && j == 2) || (i == 1 && j == 3) || (i == 4 && j == 6) || (i == 5 && j == 7)) {
    return true;
  }
  if((i == 2 && j == 0) || (i == 3 && j == 1) || (i == 6 && j == 4) || (i == 7 && j == 5)) {
    return true;
  }
  
  return false;
}

template<class TList>
void GenerateWordsOfLengthLessThan4(const TList& input, TList& output)
{
  typename TList::value_type el;
  
  std::copy(input.begin(), input.end(), std::insert_iterator<TList>(output, output.end()));
  
  typename TList::const_iterator gi, gj, gk, gl;
  gi = input.begin();
  int i, j, k;
  for(i = 0, gi = input.begin(); gi != input.end(); i++, gi++) {
    for(j = 0, gj = input.begin(); gj != input.end(); j++, gj++) {
      if(inverses(i, j)) {
        continue;
      }
      
      el.g = gi->g * gj->g;
      output.push_back(el);
    }
  }
  
  for(i = 0, gi = input.begin(); gi != input.end(); i++, gi++) {
    for(j = 0, gj = input.begin(); gj != input.end(); j++, gj++) {
      for(k = 0, gk = input.begin(); gk != input.end(); k++, gk++) {
        if(inverses(i, j) || inverses(j, k)) {
          continue;
        }
        
        el.g = gi->g * gj->g * gk->g;
        output.push_back(el);
      }
    }
  }
  
}


template<class TList>
void GenerateWordsOfLengthLessThan4_2(const TList& input, TList& output)
{
  typedef typename TList::value_type Val_type;
  
  std::copy(input.begin(), input.end(), std::insert_iterator<TList>(output, output.end()));
  
  typename TList::const_iterator gi, gj, gk, gl;
  gi = input.begin();
  int i, j, k;
  for(i = 0, gi = input.begin(); gi != input.end(); i++, gi++) {
    for(j = 0, gj = input.begin(); gj != input.end(); j++, gj++) {
      if(inverses(i, j)) {
        continue;
      }
      
      Val_type el = (*gi) * (*gj);
      output.push_back(el);
    }
  }
  
  for(i = 0, gi = input.begin(); gi != input.end(); i++, gi++) {
    for(j = 0, gj = input.begin(); gj != input.end(); j++, gj++) {
      for(k = 0, gk = input.begin(); gk != input.end(); k++, gk++) {
        if(inverses(i, j) || inverses(j, k)) {
          continue;
        }
        
        Val_type el = (*gi) * (*gj) * (*gk);
        output.push_back(el);
      }
    }
  }
  
}


MainWindow::MainWindow()
  : DemosMainWindow(), dt(K(1))
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);
  
  // Add Poincaré disk
  qreal squared_radius = to_double(dt.geom_traits().unit_circle().squared_radius());
  qreal origin_x = 0, origin_y = 0, radius = CGAL::sqrt(squared_radius), diameter = 2*radius;
  qreal left_top_corner_x = origin_x - radius;
  qreal left_top_corner_y = origin_y - radius;
  qreal width = diameter, height = diameter;
  
  disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);
  scene.addItem(disk);
  
  // add the origin to the triangulation
  dt.insert(Point_2(0, 0));
  
  // Add a GraphicItem for the Delaunay triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Add a GraphicItem for the Voronoi diagram
  vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(vgi);
  vgi->hide();

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  mp = new CGAL::Qt::TriangulationMovingPoint<Delaunay>(&dt, this);
  // TriangulationMovingPoint<Delaunay> emits a modelChanged() signal each
  // time the moving point moves.
  // The following connection is for the purpose of emitting changed().
  QObject::connect(mp, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  trv = new CGAL::Qt::TriangulationRemoveVertex<Delaunay>(&dt, this);
  QObject::connect(trv, SIGNAL(modelChanged()),
		   this, SIGNAL(changed()));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Delaunay>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  cz = new CGAL::Qt::TriangulationConflictZone<Delaunay>(&scene, &dt, this);

  // initialization of translations
  trs_a = new PointTranslation(Translations::a(), scene, &dt, this);
  trs_b = new PointTranslation(Translations::b(), scene, &dt, this);
  trs_c = new PointTranslation(Translations::c(), scene, &dt, this);
  trs_d = new PointTranslation(Translations::d(), scene, &dt, this);
  
  trs_a->info().setString(L"a");
  trs_b->info().setString(L"b");
  trs_c->info().setString(L"c");
  trs_d->info().setString(L"d");
  
  // temp things
  TranslationWithInfo tr_a(Translations::a(), TranslationInfo<std::wstring>(L"a") );
  TranslationWithInfo tr_b(Translations::b(), TranslationInfo<std::wstring>(L"b") );
  TranslationWithInfo tr_c(Translations::c(), TranslationInfo<std::wstring>(L"c") );
  TranslationWithInfo tr_d(Translations::d(), TranslationInfo<std::wstring>(L"d") );
  
  TranslationWithInfo tr_inv_a(Translations::a().inverse(), TranslationInfo<std::wstring>(L"a̅") );
  TranslationWithInfo tr_inv_b(Translations::b().inverse(), TranslationInfo<std::wstring>(L"b̅") );
  TranslationWithInfo tr_inv_c(Translations::c().inverse(), TranslationInfo<std::wstring>(L"c̅") );
  TranslationWithInfo tr_inv_d(Translations::d().inverse(), TranslationInfo<std::wstring>(L"d̅") );
                                               
  std::vector<TranslationWithInfo> translations_with_info;
  translations_with_info.push_back(tr_a);
  translations_with_info.push_back(tr_b);
  translations_with_info.push_back(tr_inv_a);
  translations_with_info.push_back(tr_inv_b);
  translations_with_info.push_back(tr_c);
  translations_with_info.push_back(tr_d);
  translations_with_info.push_back(tr_inv_c);
  translations_with_info.push_back(tr_inv_d);
  
  GenerateWordsOfLengthLessThan4_2(translations_with_info, cloud_of_translations);
  std::cout << "size of cloud " << cloud_of_translations.size() << std::endl;
  
  odnwi = new OriginalDomainNeighborsWithInfo(scene, &dt, this);
  QObject::connect(odnwi, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  //
  
  QObject::connect(trs_a, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  QObject::connect(trs_b, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  QObject::connect(trs_c, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  QObject::connect(trs_d, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  
  // experiment
  std::cout << "G: " << Translations::list().size() << std::endl;
  
  g2.add_group(&Translations::list());
  std::cout << "G_2: " << g2.number_of_elements() << std::endl;
  
  g4.add_group(&g2.list());
  std::cout << "G_4: " << g4.number_of_elements() << std::endl;
  
  g8.add_group(&g4.list());
  std::cout << "G_8: " << g8.number_of_elements() << std::endl;
  
  g16.add_group(&g8.list());
  std::cout << "G_16: " <<  g16.number_of_elements() << std::endl;
  
  g32.add_group(&g16.list());
  std::cout << "G_32: " << g32.number_of_elements() << std::endl;
  
  g64.add_group(&g32.list());
  std::cout << "G_64: " << g64.number_of_elements() << std::endl;
  
  g128.add_group(&g64.list());
  std::cout << "G_128: " << g128.number_of_elements() << std::endl;
  
  g256.add_group(&g128.list());
  std::cout << "G_256: " << g256.number_of_elements() << std::endl;
  
  odn = new OriginalDomainNeighbors(scene, &dt, this);
  QObject::connect(odn, SIGNAL(modelChanged()),
                   this, SIGNAL(changed()));
  
  GenerateWordsOfLengthLessThan4(Translations::list(), words);
  //
  
  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  ag->addAction(this->actionMovingPoint);
  ag->addAction(this->actionCircumcenter);
  ag->addAction(this->actionShowConflictZone);
  ag->addAction(this->actionDo_translation_a);
  ag->addAction(this->actionDo_translation_b);
  ag->addAction(this->actionDo_translation_c);
  ag->addAction(this->actionDo_translation_d);
  ag->addAction(this->actionG);
  ag->addAction(this->actionG2);
  ag->addAction(this->actionG4);
  ag->addAction(this->actionG8);
  ag->addAction(this->actionG16);
  
  // Check two actions 
  this->actionInsertPoint->setChecked(true);
  this->actionShowDelaunay->setChecked(true);

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  
  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Delaunay_triangulation_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));
}


void
MainWindow::processInput(CGAL::Object o)
{
  Point_2 p;
  if(CGAL::assign(p, o)){
    QPointF qp(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    
    // note that if the point is on the boundary then the disk contains the point
    //uncomment!
    //if(disk->contains(qp)){
      //dt.insert(p);
    //}
    
    //delete
    /*if(disk->contains(qp)){
      std::cout << "inserted point " << p << std::endl;
      odn2 = new OriginalDomainNeighbors(scene, &dt, this, p, 1);
      QObject::connect(odn2, SIGNAL(modelChanged()), this, SIGNAL(changed()));
      
      odn2->assign(g4.begin(), g4.end());
      odn2->assign(words.begin(), words.end());
    }*/
    
    // Oct 2012
    
    if(disk->contains(qp)){
      std::cout << "inserted point " << p << std::endl;
      odnwi2 = new OriginalDomainNeighborsWithInfo(scene, &dt, this, p, 1);
      QObject::connect(odnwi2, SIGNAL(modelChanged()),
                       this, SIGNAL(changed()));
      
      //odnwi2->assign(g4.begin(), g4.end());
      odnwi2->assign(cloud_of_translations.begin(), cloud_of_translations.end());
    }
  }
  emit(changed());
}


/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
    scene.installEventFilter(trv);
  } else {
    scene.removeEventFilter(pi);
    scene.removeEventFilter(trv);
  }
}

void
MainWindow::on_actionDo_translation_a_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(trs_a);
    trs_a->show();
  } else {
    scene.removeEventFilter(trs_a);
    trs_a->hide();
  }
}

void
MainWindow::on_actionDo_translation_b_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(trs_b);
    trs_b->show();
  } else {
    scene.removeEventFilter(trs_b);
    trs_b->hide();
  }
}

void
MainWindow::on_actionDo_translation_c_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(trs_c);
    trs_c->show();
  } else {
    scene.removeEventFilter(trs_c);
    trs_c->hide();
  }
}

void
MainWindow::on_actionDo_translation_d_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(trs_d);
    trs_d->show();
  } else {
    scene.removeEventFilter(trs_d);
    trs_d->hide();
  }
}

void
MainWindow::on_actionG_toggled(bool checked)
{
  if(checked){
    odn->assign(Translations::list_begin(), Translations::list_end());
  }
}

void
MainWindow::on_actionG2_toggled(bool checked)
{
  if(checked){
    odn->assign(g2.begin(), g2.end());
  }
}

void
MainWindow::on_actionG4_toggled(bool checked)
{
  if(checked){
    odn->assign(g4.begin(), g4.end());
  }
}

void
MainWindow::on_actionG8_toggled(bool checked)
{
  if(checked){
    odn->assign(g8.begin(), g8.end());
  }
}

void
MainWindow::on_actionG16_toggled(bool checked)
{
  if(checked){
    // odn->assign(g16.begin(), g16.end());
    
    // Delete! Add all the words of length less than 4 + some words of length 4.
    //odn->assign(g4.begin(), g4.end());    
    //odn->assign(words.begin(), words.end());
    
    // Oct 2012
    odnwi->assign(cloud_of_translations.begin(), cloud_of_translations.end());
  }
}

void
MainWindow::on_actionMovingPoint_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(mp);
  } else {
    scene.removeEventFilter(mp);
  }
}

void
MainWindow::on_actionShowConflictZone_toggled(bool checked)
{

  if(checked){
    scene.installEventFilter(cz);
  } else {
    scene.removeEventFilter(cz);
  }
}

void
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(tcc);
    tcc->show();
  } else {  
    scene.removeEventFilter(tcc);
    tcc->hide();
  }
}


void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
  dgi->setVisibleEdges(checked);
}


void
MainWindow::on_actionShowVoronoi_toggled(bool checked)
{
  vgi->setVisible(checked);
}


void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  QRectF rect = CGAL::Qt::viewportsBbox(&scene);
  CGAL::Qt::Converter<K> convert;  
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg(isor.min(), isor.max());
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInteger(this, 
                             tr("Number of random points"),
                             tr("Enter number of random points"),
			     100,
			     0,
			     std::numeric_limits<int>::max(),
			     1,
			     &ok);

  if(!ok) {
    return;
  }

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    points.push_back(*pg++);
  }
  dt.insert(points.begin(), points.end());
  // default cursor
  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points file"),
						  ".");
  if(! fileName.isEmpty()){
    open(fileName);
  }
}


void
MainWindow::open(QString fileName)
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::ifstream ifs(qPrintable(fileName));
  
  K::Point_2 p;
  std::vector<K::Point_2> points;
  while(ifs >> p) {
    points.push_back(p);
  }
  dt.insert(points.begin(), points.end());

  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  actionRecenter->trigger();
  emit(changed());
    
}

void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save points"),
						  ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    for(Delaunay::Finite_vertices_iterator 
          vit = dt.finite_vertices_begin(),
          end = dt.finite_vertices_end();
        vit!= end; ++vit)
    {
      ofs << vit->point() << std::endl;
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Hyperbolic_translations_2_demo.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Delaunay_triangulation_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }

  return app.exec();
}
