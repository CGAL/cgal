#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Scene_polyhedron_item.h"

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
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    QAction* actionCube = new QAction("Generate Cube", mw);
    QAction* actionCylinder = new QAction("Generate Cylinder", mw);
    QAction* actionSphere = new QAction("Generate Sphere", mw);
    QAction* actionTetrahedron = new QAction("Generate Tetrahedron", mw);
    connect(actionCube, SIGNAL(triggered()),
            this, SLOT(on_actionCube_triggered()));
    connect(actionCylinder, SIGNAL(triggered()),
            this, SLOT(on_actionCylinder_triggered()));
    connect(actionSphere, SIGNAL(triggered()),
            this, SLOT(on_actionSphere_triggered()));
    connect(actionTetrahedron, SIGNAL(triggered()),
            this, SLOT(on_actionTetrahedron_triggered()));
    _actions << actionCube
             << actionCylinder
             << actionSphere
             << actionTetrahedron;
    Q_FOREACH(QAction* action, _actions)
      action->setProperty("subMenuName", "Basic Volume Generation");
  }

  bool applicable(QAction*) const
  {
    return true;
  }
  QList<QAction*> actions() const {
    return _actions;
  }
public Q_SLOTS:
  void on_actionCube_triggered();
  void on_actionCylinder_triggered();
  void on_actionSphere_triggered();
  void on_actionTetrahedron_triggered();
private:
  QList<QAction*> _actions;
}; //end of class Basic_volumes_generator_plugin

//make a triangulated cube
void Basic_volumes_generator_plugin::on_actionCube_triggered()
{
  Polyhedron cube;
  VPMap vpmap = get(CGAL::vertex_point, cube);
  vertex_descriptor vertices[8];
  for(int i=0; i<8; ++i)
    vertices[i] = add_vertex(cube);

  //vertices
  put(vpmap, vertices[0],Point(-0.5,-0.5,-0.5));
  put(vpmap, vertices[1],Point(0.5,-0.5,-0.5));
  put(vpmap, vertices[2],Point(0.5,0.5,-0.5));
  put(vpmap, vertices[3],Point(-0.5,0.5,-0.5));
  put(vpmap, vertices[4],Point(-0.5,-0.5,0.5));
  put(vpmap, vertices[5],Point(0.5,-0.5,0.5));
  put(vpmap, vertices[6],Point(0.5,0.5,0.5));
  put(vpmap, vertices[7],Point(-0.5,0.5,0.5));
  //faces
  std::vector<vertex_descriptor> face;
  face.resize(3);

  face[0] = vertices[0];
  face[1] = vertices[1];
  face[2] = vertices[2];
  euler::add_face(face, cube);
  face[0] = vertices[0];
  face[1] = vertices[2];
  face[2] = vertices[3];
  euler::add_face(face, cube);
  face[0] = vertices[1];
  face[1] = vertices[5];
  face[2] = vertices[6];
  euler::add_face(face, cube);
  face[0] = vertices[1];
  face[1] = vertices[6];
  face[2] = vertices[2];
  euler::add_face(face, cube);
  face[0] = vertices[5];
  face[1] = vertices[4];
  face[2] = vertices[7];
  euler::add_face(face, cube);
  face[0] = vertices[5];
  face[1] = vertices[7];
  face[2] = vertices[6];
  euler::add_face(face, cube);
  face[0] = vertices[4];
  face[1] = vertices[0];
  face[2] = vertices[3];
  euler::add_face(face, cube);
  face[0] = vertices[4];
  face[1] = vertices[3];
  face[2] = vertices[7];
  euler::add_face(face, cube);
  face[0] = vertices[2];
  face[1] = vertices[6];
  face[2] = vertices[7];
  euler::add_face(face, cube);
  face[0] = vertices[2];
  face[1] = vertices[7];
  face[2] = vertices[3];
  euler::add_face(face, cube);
  face[0] = vertices[4];
  face[1] = vertices[5];
  face[2] = vertices[1];
  euler::add_face(face, cube);
  face[0] = vertices[4];
  face[1] = vertices[1];
  face[2] = vertices[0];
  euler::add_face(face, cube);

  Scene_polyhedron_item* cube_item = new Scene_polyhedron_item(cube);
  cube_item->setName(QString("Cube"));
  scene->addItem(cube_item);
}
//make a cylinder
void Basic_volumes_generator_plugin::on_actionCylinder_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->sphereGroupBox->setVisible(false);
  dialog->cylinderGroupBox->setVisible(true);
  //opens the dialog
  if(!dialog->exec())
    return;
  int precision = dialog->cylinderSpinBox->value();
  bool is_closed = dialog->cylinderCheckBox->isChecked();
  const float to_rad = static_cast<float>(CGAL_PI / 180.0);

  Polyhedron cylinder;
  VPMap vpmap = get(CGAL::vertex_point, cylinder);
  const int nb_vertices = 360/precision * 2;
  vertex_descriptor vertices[nb_vertices];
  for(int i=0; i<nb_vertices; ++i)
    vertices[i] = add_vertex(cylinder);

  //fill vertices
  for(int i=0; i < nb_vertices/2; ++i)
  {
    put(vpmap, vertices[i], Point(0.5*cos(i*precision*to_rad),0.5,-0.5*sin(i*precision*to_rad) ));
    put(vpmap, vertices[i+nb_vertices/2], Point(0.5*cos(i*precision*to_rad),-0.5,-0.5*sin(i*precision*to_rad)));
  }
  std::vector<vertex_descriptor> face;
  face.resize(3);
  //fill faces
  for(int i=0; i<nb_vertices/2; ++i)
  {
    face[0] = vertices[(i+1)%(nb_vertices/2)];
    face[1] = vertices[i];
    face[2] = vertices[(i+1)%(nb_vertices/2) + nb_vertices/2];
    euler::add_face(face, cylinder);

    face[0] = vertices[(i+1)%(nb_vertices/2) + nb_vertices/2];
    face[1] = vertices[i];
    face[2] = vertices[i + nb_vertices/2];
    euler::add_face(face, cylinder);
  }

  //close
  if(is_closed)
  {
    //add the center of the fans
    vertex_descriptor top = add_vertex(cylinder);
    vertex_descriptor bot = add_vertex(cylinder);
    put(vpmap, top, Point(0,0.5,0));
    put(vpmap, bot, Point(0,-0.5,0));

    //add the faces
    for(int i=0; i<nb_vertices/2; ++i)
    {
      face[0] = vertices[i];
      face[1] = vertices[(i+1)%(nb_vertices/2)];
      face[2] = top;
      euler::add_face(face, cylinder);

      face[0] = bot;
      face[1] = vertices[(i+1)%(nb_vertices/2) + nb_vertices/2];
      face[2] = vertices[i + nb_vertices/2];
      euler::add_face(face, cylinder);
    }
  }
  Scene_polyhedron_item* cylinder_item = new Scene_polyhedron_item(cylinder);
  cylinder_item->setName(QString("Cylinder"));
  scene->addItem(cylinder_item);
}
//make a sphere
void Basic_volumes_generator_plugin::on_actionSphere_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->cylinderGroupBox->setVisible(false);
  dialog->sphereGroupBox->setVisible(true);
  //opens the dialog
  if(!dialog->exec())
    return;
  int precision = dialog->SphereSpinBox->value();
  const float to_rad = static_cast<float>(CGAL_PI / 180.0);

  Polyhedron sphere;
  VPMap vpmap = get(CGAL::vertex_point, sphere);

  const int nb_vertices = 179/precision * 360/precision;
  vertex_descriptor vertices[nb_vertices];
  for(std::size_t i=0; i<nb_vertices; ++i)
    vertices[i] = add_vertex(sphere);

  //fill vertices
  int id=0;
  const int per_meridian = 180/precision;
  const int per_parallel = 360/precision;
  for(int t=1; t<per_meridian; ++t)
  {
    for(int p=0; p<per_parallel; ++p)
    {
      put(vpmap, vertices[id++], Point(0.5*sin(t*to_rad*precision)*sin(p*to_rad*precision), 0.5*cos(t*to_rad*precision), 0.5*sin(t*to_rad*precision)*cos(p*to_rad*precision)));
    }
  };
  vertex_descriptor top = add_vertex(sphere);
  vertex_descriptor bot = add_vertex(sphere);
  put(vpmap, top, Point(0,0.5,0));
  put(vpmap, bot, Point(0,-0.5,0));
  std::vector<vertex_descriptor> face;
  face.resize(3);
  //fill faces
  for(int i=0; i<per_meridian-2; ++i)
  {
    for(int j=0; j < per_parallel; ++j)
    {
      face[0] = vertices[i*per_parallel+(j+1)%per_parallel];
      face[1] = vertices[i*per_parallel+j];
      face[2] = vertices[(i+1)*per_parallel+(j+1)%per_parallel];
      euler::add_face(face, sphere);
      face[0] = vertices[(i+1)*per_parallel+(j+1)%per_parallel];
      face[1] = vertices[i*per_parallel+j];
      face[2] = vertices[(i+1)*per_parallel+j%per_parallel];
      euler::add_face(face, sphere);
    }
  }

  //fill poles
  for(int i=0; i < per_parallel; ++i)
  {
    face[0] = top;
    face[1] = vertices[i];
    face[2] = vertices[(i+1)%(per_parallel)];
    euler::add_face(face, sphere);

    face[0] = bot;
    face[1] = vertices[(i+1)%(per_parallel) + (per_meridian-2)*per_parallel];
    face[2] = vertices[i + (per_meridian-2)*per_parallel];
    euler::add_face(face, sphere);
  }
  Scene_polyhedron_item* sphere_item = new Scene_polyhedron_item(sphere);
  sphere_item->setName(QString("Sphere"));
  scene->addItem(sphere_item);
}
//make a tetrahedron
void Basic_volumes_generator_plugin::on_actionTetrahedron_triggered()
{
  Polyhedron tetrahedron;
  VPMap vpmap = get(CGAL::vertex_point, tetrahedron);
  vertex_descriptor vertices[4];
  for(int i=0; i<4; ++i)
    vertices[i] = add_vertex(tetrahedron);

  //vertices
  put(vpmap, vertices[0],Point(-0.5,-0.5,-0.5));
  put(vpmap, vertices[1],Point(0.5,-0.5,-0.5));
  put(vpmap, vertices[2],Point(0,0.5,-0.5));
  put(vpmap, vertices[3],Point(0,0,0.5));
  //faces
  std::vector<vertex_descriptor> face;
  face.resize(3);

  face[0] = vertices[1];
  face[1] = vertices[0];
  face[2] = vertices[2];
  euler::add_face(face, tetrahedron);
  face[0] = vertices[0];
  face[1] = vertices[1];
  face[2] = vertices[3];
  euler::add_face(face, tetrahedron);
  face[0] = vertices[0];
  face[1] = vertices[3];
  face[2] = vertices[2];
  euler::add_face(face, tetrahedron);
  face[0] = vertices[1];
  face[1] = vertices[2];
  face[2] = vertices[3];
  euler::add_face(face, tetrahedron);


  Scene_polyhedron_item* tet_item = new Scene_polyhedron_item(tetrahedron);
  tet_item->setName(QString("Tetrahedron"));
  scene->addItem(tet_item);
}
#include "Basic_volumes_generator_plugin.moc"
