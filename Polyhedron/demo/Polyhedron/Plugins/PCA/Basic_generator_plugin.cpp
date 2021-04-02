#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QMessageBox>
#include <QBitmap>
#include <QTabBar>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include "ui_Basic_generator_widget.h"

struct Face : public std::array<int,3>
{
  Face(int i, int j, int k)
  {
    (*this)[0] = i;
    (*this)[1] = j;
    (*this)[2] = k;
  }
};

class GeneratorWidget :
    public QDockWidget,
    public Ui::BasicGenerator
{
  Q_OBJECT
public:
  GeneratorWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};

typedef Kernel::Point_3 Point;

namespace euler =  CGAL::Euler;
using namespace CGAL::Three;
namespace params = CGAL::parameters;

class Q_DECL_EXPORT Basic_generator_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "basic_generator_plugin.json")
public :
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    for(int i=0; i<=POLYLINE; ++i)
      nbs[i]=0;


    QMenu* menuFile = mw->findChild<QMenu*>("menuFile");

    QMenu* menu = menuFile->findChild<QMenu*>("menuGenerateObject");
    if(!menu){
      QAction* actionLoad = mw->findChild<QAction*>("actionLoadPlugin");
      menu = new QMenu(tr("Generate &Object"), menuFile);
      menu->setObjectName("menuGenerateObject");
      menuFile->insertMenu(actionLoad, menu);
    }


    QAction* actionPrism       = new QAction("P&rism", mw);
    QAction* actionSphere      = new QAction("&Sphere", mw);
    QAction* actionPyramid     = new QAction("Py&ramid", mw);
    QAction* actionHexahedron  = new QAction("&Hexahedron", mw);
    QAction* actionTetrahedron = new QAction("&Tetrahedron", mw);
    QAction* actionGrid        = new QAction("&Grid", mw);
    QAction* actionPointSet    = new QAction("Po&int Set", mw);
    QAction* actionPolyline    = new QAction("Po&lyline", mw);




    connect(actionPrism, SIGNAL(triggered()),
            this, SLOT(on_actionPrism_triggered()));
    _actions << actionPrism;

    connect(actionSphere, SIGNAL(triggered()),
            this, SLOT(on_actionSphere_triggered()));
    _actions << actionSphere;

    connect(actionPyramid, SIGNAL(triggered()),
            this, SLOT(on_actionPyramid_triggered()));
    _actions << actionPyramid;

    connect(actionHexahedron, SIGNAL(triggered()),
            this, SLOT(on_actionHexahedron_triggered()));
    _actions << actionHexahedron;

    connect(actionTetrahedron, SIGNAL(triggered()),
            this, SLOT(on_actionTetrahedron_triggered()));
    _actions << actionTetrahedron;

    connect(actionGrid, SIGNAL(triggered()),
            this, SLOT(on_actionGrid_triggered()));
    _actions << actionGrid;

    connect(actionPointSet, SIGNAL(triggered()),
            this, SLOT(on_actionPointSet_triggered()));
    _actions << actionPointSet;

    connect(actionPolyline, SIGNAL(triggered()),
            this, SLOT(on_actionPolyline_triggered()));
    _actions << actionPolyline;

    Q_FOREACH(QAction* action, _actions)
    {
      menu->addAction(action);
    }
    dock_widget = new GeneratorWidget("Basic Objects", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
    connect(dock_widget->generateButton, &QAbstractButton::clicked,
            this, &Basic_generator_plugin::on_generate_clicked);
    connect(dock_widget->selector_tabWidget, &QTabWidget::currentChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget, &GeneratorWidget::visibilityChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget->prismCheckBox, &QCheckBox::stateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget->pyramidCheckBox, &QCheckBox::stateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
  }

  bool applicable(QAction*) const
  {
    //not in the Operations menu
    return false;
  }
  QList<QAction*> actions() const {
    return QList<QAction*>();
  }
public Q_SLOTS:

  void on_actionPrism_triggered();
  void on_actionSphere_triggered();
  void on_actionPyramid_triggered();
  void on_actionHexahedron_triggered();
  void on_actionTetrahedron_triggered();
  void on_actionGrid_triggered();
  void on_actionPointSet_triggered();
  void on_actionPolyline_triggered();
  void on_generate_clicked();
  void on_tab_changed();
  void closure(){ dock_widget->hide(); }
private:
  QList<QAction*> _actions;
  GeneratorWidget* dock_widget;
  enum VolumeTab
  {
    PRISM=0,
    SPHERE,
    PYRAMID,
    HEXAHEDRON,
    TETRAHEDRON,
    GRID,
    POINT_SET,
    POLYLINE
  };

  int nbs[POLYLINE+1];
  template<class Facegraph_item>
  void generateCube();
  template<class Facegraph_item>
  void generatePrism();
  template<class Facegraph_item>
  void generatePyramid();
  template<class Facegraph_item>
  void generateSphere();
  template<class Facegraph_item>
  void generateTetrahedron();
  void generatePoints();
  void generateLines();
  template<class Facegraph_item>
  void generateGrid();
}; //end of class Basic_generator_plugin

//show the widget
void Basic_generator_plugin::on_actionPrism_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(PRISM);
}
void Basic_generator_plugin::on_actionSphere_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(SPHERE);
}
void Basic_generator_plugin::on_actionPyramid_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(PYRAMID);
}
void Basic_generator_plugin::on_actionHexahedron_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(HEXAHEDRON);
}
void Basic_generator_plugin::on_actionTetrahedron_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(TETRAHEDRON);
}
void Basic_generator_plugin::on_actionGrid_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(GRID);
}
void Basic_generator_plugin::on_actionPointSet_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(POINT_SET);
}
void Basic_generator_plugin::on_actionPolyline_triggered()
{
  dock_widget->show(); dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(POLYLINE);
}
void Basic_generator_plugin::on_tab_changed()
{
  QString name;
  int nb = 0;
  switch(dock_widget->selector_tabWidget->currentIndex())
  {
  case PRISM:
  {
    name = QString("Prism");
    nb = nbs[PRISM];
    QPixmap pic;
    if(dock_widget->prismCheckBox->isChecked())
      pic = QPixmap(":/cgal/Polyhedron_3/resources/prism.png");
    else
      pic = QPixmap(":/cgal/Polyhedron_3/resources/prism-open.png");
    dock_widget->prism_picLabel->setPixmap(pic);
    dock_widget->prism_picLabel->show();
  }
    break;
  case SPHERE:
    name = QString("Sphere");
    nb = nbs[SPHERE];
    break;

  case PYRAMID:
  {
    name = QString("Pyramid");
    nb = nbs[PYRAMID];
    QPixmap pic;
    if(dock_widget->pyramidCheckBox->isChecked())
      pic = QPixmap(":/cgal/Polyhedron_3/resources/pyramid.png");
    else
      pic = QPixmap(":/cgal/Polyhedron_3/resources/pyramid-open.png");
    dock_widget->pyramid_picLabel->setPixmap(pic);
    dock_widget->pyramid_picLabel->show();
  }
    break;

  case HEXAHEDRON:
    name = QString("Hexahedron");
    nb = nbs[HEXAHEDRON];
    break;

  case TETRAHEDRON:
    name = QString("Tetrahedron");
    nb = nbs[TETRAHEDRON];
    break;
  case GRID:
    name = QString("Grid");
    nb = nbs[GRID];
    break;
  case POINT_SET:
    name = QString("Point_set");
    nb = nbs[POINT_SET];
    break;

  case POLYLINE:
    name = QString("Polyline");
    nb = nbs[POLYLINE];
    break;
  default:
    break;
  };
  dock_widget->name_lineEdit->setText((nb==0)?name:QString(name).append(QString("%1").arg(nb)));
}
//generate
void Basic_generator_plugin::on_generate_clicked()
{
  switch(dock_widget->selector_tabWidget->currentIndex())
  {
  case PRISM:
    generatePrism<Scene_surface_mesh_item>();
    ++nbs[PRISM];
    break;

  case SPHERE:
    generateSphere<Scene_surface_mesh_item>();
    ++nbs[SPHERE];
    break;

  case PYRAMID:
    generatePyramid<Scene_surface_mesh_item>();
    ++nbs[PYRAMID];
    break;

  case HEXAHEDRON:
    generateCube<Scene_surface_mesh_item>();
    ++nbs[HEXAHEDRON];
    break;

  case TETRAHEDRON:
    generateTetrahedron<Scene_surface_mesh_item>();
    ++nbs[TETRAHEDRON];
    break;
  case GRID:
    generateGrid<Scene_surface_mesh_item>();
    ++nbs[GRID];
    break;
  case POINT_SET:
    generatePoints();
    ++nbs[POINT_SET];
    break;

  case POLYLINE:
    generateLines();
    ++nbs[POLYLINE];
    break;
  default:
    break;
  };
  on_tab_changed();
}
//make a non triangle hexahedron
template<class Facegraph_item>
void Basic_generator_plugin::generateCube()
{
  typename Facegraph_item::Face_graph cube;
  if(dock_widget->tabWidget_2->currentIndex() == 0)
  {
    QString point_texts[8];
    Point points[8];
    point_texts[0] = dock_widget->cubeP0->text(); point_texts[4] = dock_widget->cubeP4->text();
    point_texts[1] = dock_widget->cubeP1->text(); point_texts[5] = dock_widget->cubeP5->text();
    point_texts[2] = dock_widget->cubeP2->text(); point_texts[6] = dock_widget->cubeP6->text();
    point_texts[3] = dock_widget->cubeP3->text(); point_texts[7] = dock_widget->cubeP7->text();

    for(int i=0; i<8; ++i)
    {
      QStringList list = point_texts[i].split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
      if (list.isEmpty()) return;
      if (list.size()!=3){
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR : Input should consists of 3 doubles.");
        msgBox->exec();
        return;
      }
      double coords[3];
      for(int j=0; j<3; ++j)
      {
        bool ok;
        coords[j] = list.at(j).toDouble(&ok);
        if(!ok)
        {
            QMessageBox *msgBox = new QMessageBox;
            msgBox->setWindowTitle("Error");
            msgBox->setText("ERROR : Coordinates are invalid.");
            msgBox->exec();
            return;
        }
      }
      points[i] = Point(coords[0], coords[1], coords[2]);
    }

    CGAL::make_hexahedron(
          points[0],
        points[1],
        points[2],
        points[3],

        points[4],
        points[5],
        points[6],
        points[7],
        cube);
  }
  else
  {
    QString text = dock_widget->extremaEdit->text();
    QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (list.isEmpty()) return;
    if (list.size()!=6){
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Input should consists of 6 doubles.");
      msgBox->exec();
      return;
    }

    for(int i=0; i<6; ++i)
    {
      bool ok;
      list.at(i).toDouble(&ok);
      if(!ok)
      {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR : Coordinates are invalid.");
        msgBox->exec();
        return;
      }
    }
    CGAL::make_hexahedron(
          Point(list.at(0).toDouble(),list.at(1).toDouble(),list.at(5).toDouble()),
          Point(list.at(3).toDouble(),list.at(1).toDouble(),list.at(5).toDouble()),
          Point(list.at(3).toDouble(),list.at(1).toDouble(),list.at(2).toDouble()),
          Point(list.at(0).toDouble(),list.at(1).toDouble(),list.at(2).toDouble()),

          Point(list.at(0).toDouble(),list.at(4).toDouble(),list.at(2).toDouble()),
          Point(list.at(0).toDouble(),list.at(4).toDouble(),list.at(5).toDouble()),
          Point(list.at(3).toDouble(),list.at(4).toDouble(),list.at(5).toDouble()),
          Point(list.at(3).toDouble(),list.at(4).toDouble(),list.at(2).toDouble()),
          cube);
  }
  Facegraph_item* cube_item = new Facegraph_item(cube);
  cube_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(cube_item);
  scene->setSelectedItem(id);
}
//make a prism
template<class Facegraph_item>
void Basic_generator_plugin::generatePrism()
{
  //gets the precision parameter
  int nb_vertices = dock_widget->prismSpinBox->value();
  double height(dock_widget->prismHeightSpinBox->value()),
      radius(dock_widget->prismBaseSpinBox->value());
  bool is_closed = dock_widget->prismCheckBox->isChecked();

  QString text = dock_widget->prism_lineEdit->text();
  QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty()) return;
  if (list.size()!=3){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of 3 doubles.");
    msgBox->exec();
    return;
  }
  double coords[3];
  for(int i=0; i<3; ++i)
  {
    bool ok;
    coords[i] = list.at(i).toDouble(&ok);
    if(!ok)
    {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Coordinates are invalid.");
      msgBox->exec();
      return;
    }
  }

  typename Facegraph_item::Face_graph prism;
  make_regular_prism(nb_vertices,
                     prism,
                     Point(coords[0], coords[1], coords[2]),
                     height,
                     radius,
                     is_closed);

  Facegraph_item* prism_item = new Facegraph_item(prism);
  prism_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(prism_item);
  scene->setSelectedItem(id);
}
//make a pyramid
template<class Facegraph_item>
void Basic_generator_plugin::generatePyramid()
{
  int nb_vertices = dock_widget->pyramidSpinBox->value();
  double height(dock_widget->pyramidHeightSpinBox->value()),
      radius(dock_widget->pyramidBaseSpinBox->value());
  bool is_closed = dock_widget->pyramidCheckBox->isChecked();

  QString text = dock_widget->pyramid_lineEdit->text();
  QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty()) return;
  if (list.size()!=3){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of 3 doubles.");
    msgBox->exec();
    return;
  }
  double coords[3];
  for(int i=0; i<3; ++i)
  {
    bool ok;
    coords[i] = list.at(i).toDouble(&ok);
    if(!ok)
    {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Coordinates are invalid.");
      msgBox->exec();
      return;
    }
  }

  typename Facegraph_item::Face_graph pyramid;
  make_pyramid(nb_vertices,
               pyramid,
               Point(coords[0], coords[1], coords[2]),
               height,
               radius,
               is_closed);

  Facegraph_item* pyramid_item = new Facegraph_item(pyramid);
  pyramid_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(pyramid_item);
  scene->setSelectedItem(id);
}
//make a sphere
template<class Facegraph_item>
void Basic_generator_plugin::generateSphere()
{
  int precision = dock_widget->SphereSpinBox->value();
  QString text = dock_widget->center_radius_lineEdit->text();
  QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty()) return;
  if (list.size()!=4){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of four doubles.");
    msgBox->exec();
    return;
  }

  double radius = list.at(3).toDouble();
  Point center = Point(list.at(0).toDouble(), list.at(1).toDouble(), list.at(2).toDouble());
  typename Facegraph_item::Face_graph sphere;
  make_icosahedron(sphere, center, radius);
  typedef typename boost::property_map<typename Facegraph_item::Face_graph, CGAL::vertex_point_t>::type VPMap;
  if(precision !=0)
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sphere,
                                                  params::number_of_iterations(precision));
  VPMap vpmap = get(CGAL::vertex_point, sphere);
  //emplace the points back on the sphere
  for(typename boost::graph_traits<typename Facegraph_item::Face_graph>::vertex_descriptor vd : vertices(sphere))
  {
    Kernel::Vector_3 vec(center, get(vpmap, vd));
    vec = radius*vec/CGAL::sqrt(vec.squared_length());
    put(vpmap, vd, Kernel::Point_3(center.x() + vec.x(),
                                   center.y() + vec.y(),
                                   center.z() + vec.z()));
  }

  Facegraph_item* sphere_item = new Facegraph_item(sphere);
  sphere_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(sphere_item);
  scene->setSelectedItem(id);
}
//make a tetrahedron
template<class Facegraph_item>
void Basic_generator_plugin::generateTetrahedron()
{
  typename Facegraph_item::Face_graph tetrahedron;

  QString point_texts[4];
  Point points[4];
  point_texts[0] = dock_widget->tetP0->text();
  point_texts[1] = dock_widget->tetP1->text();
  point_texts[2] = dock_widget->tetP2->text();
  point_texts[3] = dock_widget->tetP3->text();

  for(int i=0; i<4; ++i)
  {
    QStringList list = point_texts[i].split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (list.isEmpty()) return;
    if (list.size()!=3){
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Input should consists of 3 doubles.");
      msgBox->exec();
      return;
    }
    double coords[3];
    for(int j=0; j<3; ++j)
    {
      bool ok;
      coords[j] = list.at(j).toDouble(&ok);
      if(!ok)
      {
          QMessageBox *msgBox = new QMessageBox;
          msgBox->setWindowTitle("Error");
          msgBox->setText("ERROR : Coordinates are invalid.");
          msgBox->exec();
          return;
      }
    }
    points[i] = Point(coords[0], coords[1], coords[2]);
  }
  CGAL::make_tetrahedron(points[0],
                         points[1],
                         points[2],
                         points[3],
                         tetrahedron);

  Facegraph_item* tet_item = new Facegraph_item(tetrahedron);
  tet_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(tet_item);
  scene->setSelectedItem(id);
}
//make a point set
void Basic_generator_plugin::generatePoints()
{
  QString text = dock_widget->point_textEdit->toPlainText();
  Scene_points_with_normal_item* item = new Scene_points_with_normal_item();
  QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  int counter = 0;
  double coord[3];
  bool ok = true;
  if (list.isEmpty()) return;
  if (list.size()%3!=0){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of triplets.");
    msgBox->exec();
    return;
  }

  Q_FOREACH(QString s, list)
  {
    if(!s.isEmpty())
    {
      double res = s.toDouble(&ok);
      if(!ok)
      {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR : Coordinates are invalid.");
        msgBox->exec();
        break;
      }
      else
      {
        coord[counter] = res;
        counter++;
      }
    }
    if(counter == 3)
    {
      const Kernel::Point_3 p(coord[0], coord[1], coord[2]);
      item->point_set()->insert(p);
      counter =0;
    }
  }
  if(ok)
  {
    dock_widget->point_textEdit->clear();
    item->point_set()->unselect_all();
    item->setName(dock_widget->name_lineEdit->text());
    item->setColor(Qt::black);
    item->invalidateOpenGLBuffers();
    Scene_interface::Item_id id = scene->addItem(item);
    scene->setSelectedItem(id);
  }
}
//make a polyline
void Basic_generator_plugin::generateLines()
{
  QString text = dock_widget->line_textEdit->toPlainText();
  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;
  polylines.resize(polylines.size()+1);
  std::vector<Scene_polylines_item::Point_3>& polyline = *(polylines.rbegin());
  QStringList polylines_metadata;

  QStringList list = text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  int counter = 0;
  double coord[3];
  bool ok = true;
  if (list.isEmpty()) return;
  if(!dock_widget->polygon_checkBox->isChecked() && list.size()%3!=0){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of triplets.");
    msgBox->exec();
    return;
  }
  else if(dock_widget->polygon_checkBox->isChecked()&& list.size()%2!=0){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of pairs.");
    msgBox->exec();
    return;
  }
  Q_FOREACH(QString s, list)
  {
    if(!s.isEmpty())
    {
      double res = s.toDouble(&ok);
      if(!ok)
      {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR : Coordinates are invalid.");
        msgBox->exec();
        break;
      }
      else
      {
        coord[counter] = res;
        counter++;
      }
    }
    if(!dock_widget->polygon_checkBox->isChecked() && counter == 3)
    {
      Scene_polylines_item::Point_3 p(coord[0], coord[1], coord[2]);
      polyline.push_back(p);
      counter =0;
    }
    else if(dock_widget->polygon_checkBox->isChecked() && counter == 2)
    {
      Scene_polylines_item::Point_3 p(coord[0], coord[1], 0);
      polyline.push_back(p);
      counter = 0;
    }
  }
  if(dock_widget->polygon_checkBox->isChecked())
  {
    polyline.push_back(polyline.front()); //polygon_2 are not closed.
  }
  if(ok)
  {
    dock_widget->line_textEdit->clear();
    if(dock_widget->fill_checkBox->isChecked())
    {
      CGAL::Three::Three::CursorScopeGuard guard(Qt::WaitCursor);
      QApplication::processEvents();
      if(polyline.front() != polyline.back()) {
        polyline.push_back(polyline.front());
      }
      if(polyline.size() < 4) { // no triangle, skip it (needs at least 3 + 1 repeat)
        QMessageBox::warning(mw, "Warning", "Needs at least 3 points to triangulate. Aborting.");
        return;
      }
      std::vector<Face> patch;
      CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polyline,
                                                               std::back_inserter(patch),
                                                               CGAL::Polygon_mesh_processing::parameters::use_delaunay_triangulation(true));

      if(patch.empty()) {
          QMessageBox::warning(mw, "Warning", "Triangulation failed.");
          return;
      }
      SMesh* poly = new SMesh;
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(polyline,
                                                                  patch,
                                                                  *poly);

      Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item(poly);
      poly_item->setName(dock_widget->name_lineEdit->text());
      poly_item->setRenderingMode(FlatPlusEdges);
      scene->setSelectedItem(scene->addItem(poly_item));
    }
    else
    {
      Scene_polylines_item* item = new Scene_polylines_item();
      item->polylines = polylines;
      item->invalidateOpenGLBuffers();
      item->setName(dock_widget->name_lineEdit->text());
      item->setColor(Qt::black);
      item->setProperty("polylines metadata", polylines_metadata);
      Scene_interface::Item_id id = scene->addItem(item);
      scene->setSelectedItem(id);
    }
  }
}

struct Point_generator
{
  std::size_t w,h;
  Point bl, ur;
  Point_generator(const std::size_t& w,
                  const std::size_t& h,
                  const Point& bl,
                  const Point& ur)
    :w(w), h(h), bl(bl), ur(ur)
  {}
  Point operator()(std::size_t i, std::size_t j) const
  {
    return Point(bl.x() + i*(ur.x()-bl.x())/(w-1),
                 bl.y() + j*(ur.y()-bl.y())/(h-1),
                 bl.z() + j*(ur.z()-bl.z())/(h-1));
  }
};
template<class Facegraph_item>
void Basic_generator_plugin::generateGrid()
{
  typedef typename Facegraph_item::Face_graph Face_graph;
  Face_graph grid;

  QString points_text;
  Point extrema[2];
  using size_type = typename boost::graph_traits<Face_graph>::vertices_size_type;
  size_type nb_cells[2];
  bool triangulated = dock_widget->grid_checkBox->isChecked();
  points_text= dock_widget->grid_lineEdit->text();

  QStringList list = points_text.split(QRegExp("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty()) return;
  if (list.size()!=6){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of 6 doubles.");
    msgBox->exec();
    return;
  }
  double coords[6];
  for(int j=0; j<6; ++j)
  {
    bool ok;
    coords[j] = list.at(j).toDouble(&ok);
    if(!ok)
    {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Coordinates are invalid.");
      msgBox->exec();
      return;
    }
  }
  extrema[0] = Point(coords[0], coords[1], coords[2]);
  extrema[1] = Point(coords[3], coords[4], coords[5]);
  nb_cells[0] = static_cast<size_type>(dock_widget->gridX_spinBox->value());
  nb_cells[1] = static_cast<size_type>(dock_widget->gridY_spinBox->value());

  //nb_points = nb_cells+1
  Point_generator point_gen(nb_cells[0]+1, nb_cells[1]+1, extrema[0], extrema[1]);

  CGAL::make_grid(nb_cells[0], nb_cells[1], grid, point_gen, triangulated);
  Facegraph_item* grid_item = new Facegraph_item(grid);
  grid_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(grid_item);
  scene->setSelectedItem(id);
}
#include "Basic_generator_plugin.moc"
