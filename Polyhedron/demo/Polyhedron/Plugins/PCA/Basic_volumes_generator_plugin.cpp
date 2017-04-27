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
    QAction* actionPrism = new QAction("Generate Regular Prism", mw);
    QAction* actionSphere = new QAction("Generate Sphere", mw);
    QAction* actionTetrahedron = new QAction("Generate Tetrahedron", mw);
    connect(actionCube, SIGNAL(triggered()),
            this, SLOT(on_actionCube_triggered()));
    connect(actionPrism, SIGNAL(triggered()),
            this, SLOT(on_actionPrism_triggered()));
    connect(actionSphere, SIGNAL(triggered()),
            this, SLOT(on_actionSphere_triggered()));
    connect(actionTetrahedron, SIGNAL(triggered()),
            this, SLOT(on_actionTetrahedron_triggered()));
    _actions << actionCube
             << actionPrism
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
  void on_actionPrism_triggered();
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
//make a prism
template<class Mesh, class Traits>
void make_regular_prism(int nb_vertices,
                        Mesh& prism,
                        typename CGAL::Point_3<Traits> center = CGAL::Point_3<Traits>(0,0,0),
                        double height = 1.0,
                        double radius = 1.0,
                        bool is_closed = true)
{
  const double to_rad = static_cast<float>(CGAL_PI / 180.0);
  const double precision = 360/nb_vertices;
  double diameter = 2*radius;
  VPMap vpmap = get(CGAL::vertex_point, prism);
  //const int nb_vertices = 360/precision * 2;
  std::vector<vertex_descriptor> vertices;
  vertices.resize(nb_vertices*2);
  for(int i=0; i<nb_vertices*2; ++i)
    vertices[i] = add_vertex(prism);

  //fill vertices
  for(int i=0; i < nb_vertices; ++i)
  {
    put(vpmap, vertices[i], Point(0.5*diameter*cos(i*precision*to_rad)+center.x(),0.5*height+center.y(),-0.5*diameter*sin(i*precision*to_rad) + center.z()));
    put(vpmap, vertices[i+nb_vertices], Point(0.5*diameter*cos(i*precision*to_rad)+center.x(),-0.5*height+center.y(),-0.5*diameter*sin(i*precision*to_rad)+center.z()));
  }
  std::vector<vertex_descriptor> face;
  face.resize(3);
  //fill faces
  for(int i=0; i<nb_vertices; ++i)
  {
    face[0] = vertices[(i+1)%(nb_vertices)];
    face[1] = vertices[i];
    face[2] = vertices[(i+1)%(nb_vertices) + nb_vertices];
    euler::add_face(face, prism);

    face[0] = vertices[(i+1)%(nb_vertices) + nb_vertices];
    face[1] = vertices[i];
    face[2] = vertices[i + nb_vertices];
    euler::add_face(face, prism);
  }

  //close
  if(is_closed)
  {
    //add the center of the fans
    vertex_descriptor top = add_vertex(prism);
    vertex_descriptor bot = add_vertex(prism);
    put(vpmap, top, Point(center.x(),0.5*height+center.y(),center.z()));
    put(vpmap, bot, Point(center.x(),-0.5*height+center.y(),center.z()));

    //add the faces
    for(int i=0; i<nb_vertices; ++i)
    {
      face[0] = vertices[i];
      face[1] = vertices[(i+1)%(nb_vertices)];
      face[2] = top;
      euler::add_face(face, prism);

      face[0] = bot;
      face[1] = vertices[(i+1)%(nb_vertices) + nb_vertices];
      face[2] = vertices[i + nb_vertices];
      euler::add_face(face, prism);
    }
  }
}

template<class Mesh, class Traits>
void make_icosahedron(Mesh& mesh,
                      typename CGAL::Point_3<Traits> center = CGAL::Point_3<Traits>(0,0,0),
                      double radius = 1.0)
{
  VPMap vpmap = get(CGAL::vertex_point, mesh);
  // create the initial icosahedron
  std::vector<vertex_descriptor> v_vertices;
  v_vertices.resize(12);
  for(int i=0; i<12; ++i)
    v_vertices[i] = add_vertex(mesh);
  double t = (radius + radius*CGAL::sqrt(5.0)) / 2.0;

  put(vpmap, v_vertices[0],typename Traits::Point_3(-radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[1],typename Traits::Point_3( radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[2],typename Traits::Point_3(-radius + center.x(), -t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[3],typename Traits::Point_3( radius + center.x(), -t + center.y(), 0.0 + center.z()));

  put(vpmap, v_vertices[4],typename Traits::Point_3( 0.0 + center.x(), -radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[5],typename Traits::Point_3( 0.0 + center.x(),  radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[6],typename Traits::Point_3( 0.0 + center.x(), -radius + center.y(), -t + center.z()));
  put(vpmap, v_vertices[7],typename Traits::Point_3( 0.0 + center.x(),  radius + center.y(), -t + center.z()));

  put(vpmap, v_vertices[8],typename Traits::Point_3(  t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[9],typename Traits::Point_3(  t + center.x(), 0.0 + center.y(),  radius + center.z()));
  put(vpmap, v_vertices[10],typename Traits::Point_3(-t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[11],typename Traits::Point_3(-t + center.x(), 0.0 + center.y(),  radius + center.z()));

  std::vector<vertex_descriptor> face;
  face.resize(3);
  face[1] = v_vertices[0]; face[0] = v_vertices[11]; face[2] = v_vertices[5];
  euler::add_face(face, mesh);
  face[1] = v_vertices[0]; face[0] = v_vertices[5]; face[2] = v_vertices[1];
  euler::add_face(face, mesh);
  face[1] = v_vertices[0]; face[0] = v_vertices[1]; face[2] = v_vertices[7];
  euler::add_face(face, mesh);
  face[1] = v_vertices[0]; face[0] = v_vertices[7]; face[2] = v_vertices[10];
  euler::add_face(face, mesh);
  face[1] = v_vertices[0]; face[0] = v_vertices[10]; face[2] = v_vertices[11];
  euler::add_face(face, mesh);

  face[1] = v_vertices[1] ; face[0] = v_vertices[5] ; face[2] = v_vertices[9];
  euler::add_face(face, mesh);
  face[1] = v_vertices[5] ; face[0] = v_vertices[11]; face[2] = v_vertices[4];
  euler::add_face(face, mesh);
  face[1] = v_vertices[11]; face[0] = v_vertices[10]; face[2] = v_vertices[2];
  euler::add_face(face, mesh);
  face[1] = v_vertices[10]; face[0] = v_vertices[7] ; face[2] = v_vertices[6];
  euler::add_face(face, mesh);
  face[1] = v_vertices[7] ; face[0] = v_vertices[1] ; face[2] = v_vertices[8];
  euler::add_face(face, mesh);

  face[1] = v_vertices[3] ; face[0] = v_vertices[9] ; face[2] = v_vertices[4];
  euler::add_face(face, mesh);
  face[1] = v_vertices[3] ; face[0] = v_vertices[4] ; face[2] = v_vertices[2];
  euler::add_face(face, mesh);
  face[1] = v_vertices[3] ; face[0] = v_vertices[2] ; face[2] = v_vertices[6];
  euler::add_face(face, mesh);
  face[1] = v_vertices[3] ; face[0] = v_vertices[6] ; face[2] = v_vertices[8];
  euler::add_face(face, mesh);
  face[1] = v_vertices[3] ; face[0] = v_vertices[8] ; face[2] = v_vertices[9];
  euler::add_face(face, mesh);

  face[1] = v_vertices[4] ; face[0] = v_vertices[9] ; face[2] = v_vertices[5] ;
  euler::add_face(face, mesh);
  face[1] = v_vertices[2] ; face[0] = v_vertices[4] ; face[2] = v_vertices[11];
  euler::add_face(face, mesh);
  face[1] = v_vertices[6] ; face[0] = v_vertices[2] ; face[2] = v_vertices[10];
  euler::add_face(face, mesh);
  face[1] = v_vertices[8] ; face[0] = v_vertices[6] ; face[2] = v_vertices[7] ;
  euler::add_face(face, mesh);
  face[1] = v_vertices[9] ; face[0] = v_vertices[8] ; face[2] = v_vertices[1] ;
  euler::add_face(face, mesh);
}
void Basic_volumes_generator_plugin::on_actionPrism_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->sphereGroupBox->setVisible(false);
  dialog->prismGroupBox->setVisible(true);
  //opens the dialog
  if(!dialog->exec())
    return;
  int nb_vertices = dialog->prismSpinBox->value();
  double height(dialog->heightSpinBox->value()),
      radius(dialog->prismBaseSpinBox->value()),
      center_x(dialog->prismXSpinBox->value()),
      center_y(dialog->prismYSpinBox->value()),
      center_z(dialog->prismZSpinBox->value());
  bool is_closed = dialog->prismCheckBox->isChecked();

  Polyhedron prism;
  make_regular_prism(nb_vertices, prism, Point(center_x,center_y,center_z), height, radius, is_closed);
  Scene_polyhedron_item* prism_item = new Scene_polyhedron_item(prism);
  prism_item->setName(QString("Prism"));
  scene->addItem(prism_item);
}
//make a sphere
void Basic_volumes_generator_plugin::on_actionSphere_triggered()
{
  //gets the precision parameter
  VolumeDialog *dialog = new VolumeDialog();
  dialog->prismGroupBox->setVisible(false);
  dialog->sphereGroupBox->setVisible(true);
  //opens the dialog
  if(!dialog->exec())
    return;
  int precision = dialog->SphereSpinBox->value();
  double radius = dialog->sphereRadiusSpinBox->value();
  Point center(dialog->sphereXSpinBox->value(),
               dialog->sphereYSpinBox->value(),
               dialog->sphereZSpinBox->value());
  Polyhedron sphere;
  make_icosahedron(sphere, center, radius);
  if(precision !=0)
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sphere, precision);
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
