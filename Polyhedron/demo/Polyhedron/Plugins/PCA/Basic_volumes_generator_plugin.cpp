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

//make a non triangle cube
void Basic_volumes_generator_plugin::on_actionCube_triggered()
{
  Polyhedron cube;
  CGAL::make_hexahedron(Point(-0.5,-0.5,-0.5),
                        Point(0.5,-0.5,-0.5),
                        Point(0.5,0.5,-0.5),
                        Point(-0.5,0.5,-0.5),

                        Point(-0.5,0.5,0.5),
                        Point(-0.5,-0.5,0.5),
                        Point(0.5,-0.5,0.5),
                        Point(0.5,0.5,0.5),
                        cube);


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
  int nb_vertices = dialog->cylinderSpinBox->value();
  float precision = 360*2.0/nb_vertices;
  //if(360%precision != 0)
  //  QMessageBox::warning(mw, tr("WARNING"),
  //                       tr("Precision should divide 360 for a better result"),
  //                       QMessageBox::Ok);
  bool is_closed = dialog->cylinderCheckBox->isChecked();
  const float to_rad = static_cast<float>(CGAL_PI / 180.0);

  Polyhedron cylinder;
  VPMap vpmap = get(CGAL::vertex_point, cylinder);
  //const int nb_vertices = 360/precision * 2;
  std::vector<vertex_descriptor> vertices;
  vertices.resize(nb_vertices);
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
  Polyhedron sphere;
  VPMap vpmap = get(CGAL::vertex_point, sphere);
  // create the initial icosahedron
  std::vector<vertex_descriptor> v_vertices;
  v_vertices.resize(12);
  for(int i=0; i<12; ++i)
    v_vertices[i] = add_vertex(sphere);
  float t = (1.0 + CGAL::sqrt(5.0)) / 2.0;

  put(vpmap, v_vertices[0],Kernel::Point_3(-1,  t,  0));
  put(vpmap, v_vertices[1],Kernel::Point_3( 1,  t,  0));
  put(vpmap, v_vertices[2],Kernel::Point_3(-1, -t,  0));
  put(vpmap, v_vertices[3],Kernel::Point_3( 1, -t,  0));

  put(vpmap, v_vertices[4],Kernel::Point_3( 0, -1,  t));
  put(vpmap, v_vertices[5],Kernel::Point_3( 0,  1,  t));
  put(vpmap, v_vertices[6],Kernel::Point_3( 0, -1, -t));
  put(vpmap, v_vertices[7],Kernel::Point_3( 0,  1, -t));

  put(vpmap, v_vertices[8],Kernel::Point_3( t,  0, -1));
  put(vpmap, v_vertices[9],Kernel::Point_3( t,  0,  1));
  put(vpmap, v_vertices[10],Kernel::Point_3(-t, 0, -1));
  put(vpmap, v_vertices[11],Kernel::Point_3(-t, 0,  1));

  std::vector<vertex_descriptor> face;
  face.resize(3);
  face[1] = v_vertices[0]; face[0] = v_vertices[11]; face[2] = v_vertices[5];
  euler::add_face(face, sphere);
  face[1] = v_vertices[0]; face[0] = v_vertices[5]; face[2] = v_vertices[1];
  euler::add_face(face, sphere);
  face[1] = v_vertices[0]; face[0] = v_vertices[1]; face[2] = v_vertices[7];
  euler::add_face(face, sphere);
  face[1] = v_vertices[0]; face[0] = v_vertices[7]; face[2] = v_vertices[10];
  euler::add_face(face, sphere);
  face[1] = v_vertices[0]; face[0] = v_vertices[10]; face[2] = v_vertices[11];
  euler::add_face(face, sphere);

  face[1] = v_vertices[1] ; face[0] = v_vertices[5] ; face[2] = v_vertices[9];
  euler::add_face(face, sphere);
  face[1] = v_vertices[5] ; face[0] = v_vertices[11]; face[2] = v_vertices[4];
  euler::add_face(face, sphere);
  face[1] = v_vertices[11]; face[0] = v_vertices[10]; face[2] = v_vertices[2];
  euler::add_face(face, sphere);
  face[1] = v_vertices[10]; face[0] = v_vertices[7] ; face[2] = v_vertices[6];
  euler::add_face(face, sphere);
  face[1] = v_vertices[7] ; face[0] = v_vertices[1] ; face[2] = v_vertices[8];
  euler::add_face(face, sphere);

  face[1] = v_vertices[3] ; face[0] = v_vertices[9] ; face[2] = v_vertices[4];
  euler::add_face(face, sphere);
  face[1] = v_vertices[3] ; face[0] = v_vertices[4] ; face[2] = v_vertices[2];
  euler::add_face(face, sphere);
  face[1] = v_vertices[3] ; face[0] = v_vertices[2] ; face[2] = v_vertices[6];
  euler::add_face(face, sphere);
  face[1] = v_vertices[3] ; face[0] = v_vertices[6] ; face[2] = v_vertices[8];
  euler::add_face(face, sphere);
  face[1] = v_vertices[3] ; face[0] = v_vertices[8] ; face[2] = v_vertices[9];
  euler::add_face(face, sphere);

  face[1] = v_vertices[4] ; face[0] = v_vertices[9] ; face[2] = v_vertices[5] ;
  euler::add_face(face, sphere);
  face[1] = v_vertices[2] ; face[0] = v_vertices[4] ; face[2] = v_vertices[11];
  euler::add_face(face, sphere);
  face[1] = v_vertices[6] ; face[0] = v_vertices[2] ; face[2] = v_vertices[10];
  euler::add_face(face, sphere);
  face[1] = v_vertices[8] ; face[0] = v_vertices[6] ; face[2] = v_vertices[7] ;
  euler::add_face(face, sphere);
  face[1] = v_vertices[9] ; face[0] = v_vertices[8] ; face[2] = v_vertices[1] ;
  euler::add_face(face, sphere);
  if(precision !=0)
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sphere, precision);
  BOOST_FOREACH(vertex_descriptor vd, vertices(sphere))
  {
   Kernel::Vector_3 vec(get(vpmap, vd), Kernel::Point_3(0,0,0));
   vec = vec/CGAL::sqrt(vec.squared_length());
   put(vpmap, vd, Kernel::Point_3(vec.x(), vec.y(), vec.z()));
  }
  Scene_polyhedron_item* sphere_item = new Scene_polyhedron_item(sphere);
  sphere_item->setName(QString("Sphere"));
  scene->addItem(sphere_item);
}
//make a tetrahedron
void Basic_volumes_generator_plugin::on_actionTetrahedron_triggered()
{
  Polyhedron tetrahedron;
  CGAL::make_tetrahedron(Point(-0.5,-0.5,-0.5),
                         Point(0.5,-0.5,-0.5),
                         Point(0,0.5,-0.5),
                         Point(0,0,0.5),
                         tetrahedron);

  Scene_polyhedron_item* tet_item = new Scene_polyhedron_item(tetrahedron);
  tet_item->setName(QString("Tetrahedron"));
  scene->addItem(tet_item);
}
#include "Basic_volumes_generator_plugin.moc"
