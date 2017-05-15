#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QMessageBox>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Scene_polyhedron_item.h"
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Kernel_traits.h>
#include "ui_Basic_volumes_generator_dialog.h"

class VolumeDialog :
    public QDialog,
    public Ui::VolumeDialog
{
  Q_OBJECT
public:
  VolumeDialog(QWidget* =0)
  {
    setupUi(this);
  }
};

typedef Kernel::Point_3 Point;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type VPMap;

namespace euler =  CGAL::Euler;
using namespace CGAL::Three;
class Q_DECL_EXPORT Basic_volumes_generator_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    QAction* actionCube = new QAction("Generate Hexahedron", mw);
    QAction* actionPrism = new QAction("Generate Regular Prism", mw);
    QAction* actionPyramid = new QAction("Generate Pyramid", mw);
    QAction* actionSphere = new QAction("Generate Sphere", mw);
    QAction* actionTetrahedron = new QAction("Generate Tetrahedron", mw);
    connect(actionCube, SIGNAL(triggered()),
            this, SLOT(on_actionCube_triggered()));
    connect(actionPrism, SIGNAL(triggered()),
            this, SLOT(on_actionPrism_triggered()));
    connect(actionPyramid, SIGNAL(triggered()),
            this, SLOT(on_actionPyramid_triggered()));
    connect(actionSphere, SIGNAL(triggered()),
            this, SLOT(on_actionSphere_triggered()));
    connect(actionTetrahedron, SIGNAL(triggered()),
            this, SLOT(on_actionTetrahedron_triggered()));
    _actions << actionCube
             << actionPrism
             << actionPyramid
             << actionSphere
             << actionTetrahedron;
    Q_FOREACH(QAction* action, _actions)
      action->setProperty("subMenuName", "Basic Volume Generation");
  }

  bool applicable(QAction*) const
  {
    //always applicable
    return true;
  }
  QList<QAction*> actions() const {
    return _actions;
  }
public Q_SLOTS:
  void on_actionCube_triggered();
  void on_actionPrism_triggered();
  void on_actionPyramid_triggered();
  void on_actionSphere_triggered();
  void on_actionTetrahedron_triggered();
private:
  QList<QAction*> _actions;
}; //end of class Basic_volumes_generator_plugin

//make a non triangle hexahedron
void Basic_volumes_generator_plugin::on_actionCube_triggered()
{
  VolumeDialog *dialog = new VolumeDialog();
  dialog->stackedWidget->setCurrentIndex(3);
  //opens the dialog
  if(!dialog->exec())
    return;

  Polyhedron cube;
  if(dialog->tabWidget_2->currentIndex() == 0)
  {
    CGAL::make_hexahedron(Point(dialog->p0xBox->value(), dialog->p0yBox->value(), dialog->p0zBox->value()),
                          Point(dialog->p1xBox->value(), dialog->p1yBox->value(), dialog->p1zBox->value()),
                          Point(dialog->p2xBox->value(), dialog->p2yBox->value(), dialog->p2zBox->value()),
                          Point(dialog->p3xBox->value(), dialog->p3yBox->value(), dialog->p3zBox->value()),

                          Point(dialog->p4xBox->value(), dialog->p4yBox->value(), dialog->p4zBox->value()),
                          Point(dialog->p5xBox->value(), dialog->p5yBox->value(), dialog->p5zBox->value()),
                          Point(dialog->p6xBox->value(), dialog->p6yBox->value(), dialog->p6zBox->value()),
                          Point(dialog->p7xBox->value(), dialog->p7yBox->value(), dialog->p7zBox->value()),
                          cube);
  }
  else
  {
    QString text = dialog->extremaEdit->text();
    QStringList list = text.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    if (list.isEmpty()) return;
    if (list.size()!=6){
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Input should consists of 6 doubles.");
      msgBox->exec();
      return;
    }

    CGAL::make_hexahedron(Point(list.at(0).toDouble(),list.at(1).toDouble(),list.at(2).toDouble()),
                          Point(list.at(3).toDouble(),list.at(1).toDouble(),list.at(2).toDouble()),
                          Point(list.at(3).toDouble(),list.at(1).toDouble(),list.at(5).toDouble()),
                          Point(list.at(0).toDouble(),list.at(1).toDouble(),list.at(5).toDouble()),

                          Point(list.at(0).toDouble(),list.at(4).toDouble(),list.at(5).toDouble()),
                          Point(list.at(0).toDouble(),list.at(4).toDouble(),list.at(2).toDouble()),
                          Point(list.at(3).toDouble(),list.at(4).toDouble(),list.at(2).toDouble()),
                          Point(list.at(3).toDouble(),list.at(4).toDouble(),list.at(5).toDouble()),
                          cube);
  }
  Scene_polyhedron_item* cube_item = new Scene_polyhedron_item(cube);
  cube_item->setName(QString("Hexahedron"));
  scene->addItem(cube_item);
}
//make a prism
void Basic_volumes_generator_plugin::on_actionPrism_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->stackedWidget->setCurrentIndex(0);
  //opens the dialog
  if(!dialog->exec())
    return;
  int nb_vertices = dialog->prismSpinBox->value();
  double height(dialog->prismHeightSpinBox->value()),
      radius(dialog->prismBaseSpinBox->value()),
      center_x(dialog->prismXSpinBox->value()),
      center_y(dialog->prismYSpinBox->value()),
      center_z(dialog->prismZSpinBox->value());
  bool is_closed = dialog->prismCheckBox->isChecked();

  Polyhedron prism;
  make_regular_prism(nb_vertices,
                     prism,
                     Point(center_x,
                           center_y,
                           center_z),
                     height,
                     radius,
                     is_closed);

  Scene_polyhedron_item* prism_item = new Scene_polyhedron_item(prism);
  prism_item->setName(QString("Prism"));
  scene->addItem(prism_item);
}

//make a pyramid
void Basic_volumes_generator_plugin::on_actionPyramid_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->stackedWidget->setCurrentIndex(2);

  //opens the dialog
  if(!dialog->exec())
    return;
  int nb_vertices = dialog->pyramidSpinBox->value();
  double height(dialog->pyramidHeightSpinBox->value()),
      radius(dialog->pyramidBaseSpinBox->value()),
      center_x(dialog->pyramidXSpinBox->value()),
      center_y(dialog->pyramidYSpinBox->value()),
      center_z(dialog->pyramidZSpinBox->value());
  bool is_closed = dialog->pyramidCheckBox->isChecked();

  Polyhedron pyramid;
  make_pyramid(nb_vertices,
                     pyramid,
                     Point(center_x,
                           center_y,
                           center_z),
                     height,
                     radius,
                     is_closed);

  Scene_polyhedron_item* pyramid_item = new Scene_polyhedron_item(pyramid);
  pyramid_item->setName(QString("Pyramid"));
  scene->addItem(pyramid_item);
}

//make a sphere
void Basic_volumes_generator_plugin::on_actionSphere_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->stackedWidget->setCurrentIndex(1);
  //opens the dialog
  if(!dialog->exec())
    return;
  int precision = dialog->SphereSpinBox->value();
  double radius;
  Point center;
  if(dialog->tabWidget->currentIndex() == 0)
  {
    radius = dialog->sphereRadiusSpinBox->value();
    center = Point(dialog->sphereXSpinBox->value(),
               dialog->sphereYSpinBox->value(),
               dialog->sphereZSpinBox->value());
  }
  else
  {
    QString text = dialog->center_radius_lineEdit->text();
    QStringList list = text.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    if (list.isEmpty()) return;
    if (list.size()!=4){
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR : Input should consists of four doubles.");
      msgBox->exec();
      return;
    }
    radius = list.at(3).toDouble();
    center = Point(list.at(0).toDouble(), list.at(1).toDouble(), list.at(2).toDouble());
  }
  Polyhedron sphere;
  make_icosahedron(sphere, center, radius);
  if(precision !=0)
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sphere,
                                                  precision);
  VPMap vpmap = get(CGAL::vertex_point, sphere);
  //emplace the points back on the sphere
  BOOST_FOREACH(vertex_descriptor vd, vertices(sphere))
  {
    Kernel::Vector_3 vec(get(vpmap, vd), center);
    vec = radius*vec/CGAL::sqrt(vec.squared_length());
    put(vpmap, vd, Kernel::Point_3(center.x() + vec.x(),
                                   center.y() + vec.y(),
                                   center.z() + vec.z()));
  }

  Scene_polyhedron_item* sphere_item = new Scene_polyhedron_item(sphere);
  sphere_item->setName(QString("Sphere"));
  scene->addItem(sphere_item);
}
//make a tetrahedron
void Basic_volumes_generator_plugin::on_actionTetrahedron_triggered()
{
  Polyhedron tetrahedron;
  CGAL::make_tetrahedron(Point(0.0, 0.0, 0.0),
                         Point(1.0, 0.0, 0.0),
                         Point(0.0, 1.0, 0.0),
                         Point(0.0, 0.0, 1.0),
                         tetrahedron);

  Scene_polyhedron_item* tet_item = new Scene_polyhedron_item(tetrahedron);
  tet_item->setName(QString("Tetrahedron"));
  scene->addItem(tet_item);
}
#include "Basic_volumes_generator_plugin.moc"
