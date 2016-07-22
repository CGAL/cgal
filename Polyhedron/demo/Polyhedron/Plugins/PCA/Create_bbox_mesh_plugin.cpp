#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/gl.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>


const int cube[][3] = { { 0, 1, 3 },
                        { 3, 1, 2 },
                        { 0, 4, 1 },
                        { 1, 4, 5 },
                        { 3, 2, 7 },
                        { 7, 2, 6 },
                        { 4, 0, 3 },
                        { 7, 4, 3 },
                        { 6, 4, 7 },
                        { 6, 5, 4 },
                        { 1, 5, 6 },
                        { 2, 1, 6 } };

using namespace CGAL::Three;

struct Build_bbox_mesh : 
  public CGAL::Modifier_base<Polyhedron::HalfedgeDS>
{
  Scene_interface::Bbox bbox;
  typedef Polyhedron::HalfedgeDS HDS;

public:
  Build_bbox_mesh(Scene_interface::Bbox b) : bbox(b) {}

  void operator()( HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( 8, 12, 24);
    typedef HDS::Vertex   Vertex;
    typedef Vertex::Point Point;
    B.add_vertex( Point( bbox.xmin(), bbox.ymin(), bbox.zmin())); // -1 -1 -1
    B.add_vertex( Point( bbox.xmin(), bbox.ymax(), bbox.zmin())); // -1 1 -1
    B.add_vertex( Point( bbox.xmax(), bbox.ymax(), bbox.zmin())); // 1 1 -1
    B.add_vertex( Point( bbox.xmax(), bbox.ymin(), bbox.zmin())); // 1 -1 -1
    B.add_vertex( Point( bbox.xmin(), bbox.ymin(), bbox.zmax())); // -1 -1 1
    B.add_vertex( Point( bbox.xmin(), bbox.ymax(), bbox.zmax())); // -1 1 1
    B.add_vertex( Point( bbox.xmax(), bbox.ymax(), bbox.zmax())); // 1 1 1
    B.add_vertex( Point( bbox.xmax(), bbox.ymin(), bbox.zmax())); // 1 -1 1
    for(int i = 0; i < 12; ++i) {
      B.begin_facet();
      B.add_vertex_to_facet( cube[i][0]);
      B.add_vertex_to_facet( cube[i][1]);
      B.add_vertex_to_facet( cube[i][2]);
      B.end_facet();
    }
    B.end_surface();
  }
};

class Create_bbox_mesh_plugin : 
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const;
  bool applicable(QAction*) const {
    if(scene->mainSelectionIndex() != -1)
      return true;
  return false;}

protected:
  void bbox(bool extended = false);

public Q_SLOTS:
  void createBbox() {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    bbox();
    QApplication::restoreOverrideCursor();
  }
  void createExtendedBbox() {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    bbox(true);
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface* scene;
  QAction* actionBbox;
  QAction* actionExtendedBbox;

}; // end Create_bbox_mesh_plugin

void Create_bbox_mesh_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  actionBbox = new QAction(tr("Create &Bbox Mesh"), mainWindow);
  actionBbox->setObjectName("createBboxMeshAction");
  connect(actionBbox, SIGNAL(triggered()),
          this, SLOT(createBbox()));
  actionExtendedBbox = new QAction(tr("Create &Extended Bbox Mesh"), mainWindow);
  actionExtendedBbox->setObjectName("createExtendedBboxMeshAction");
  connect(actionExtendedBbox, SIGNAL(triggered()),
          this, SLOT(createExtendedBbox()));
}

QList<QAction*> Create_bbox_mesh_plugin::actions() const {
  return QList<QAction*>() << actionBbox << actionExtendedBbox;
}

void Create_bbox_mesh_plugin::bbox(bool extended)
{
  Scene_interface::Bbox bbox;
  bool initialized = false;

  Q_FOREACH(int index, scene->selectionIndices()) {
    Scene_item* item = scene->item(index);
    if(item->isFinite() && ! item->isEmpty()) {
      if(initialized) {
        bbox = bbox + item->bbox();
      } else {
        bbox = item->bbox();
        initialized = true;
      }
    }
  }
  std::cerr << "bbox dimensions: " << bbox.xmax() - bbox.xmin()
            << "\n                 " << bbox.ymax() - bbox.ymin()
            << "\n                 " << bbox.zmax() - bbox.zmin()
            << std::endl;

  if(extended) {
    const double delta_x = ( bbox.xmax() - bbox.xmin() ) / 20.;
    const double delta_y = ( bbox.ymax() - bbox.ymin() ) / 20.;
    const double delta_z = ( bbox.zmax() - bbox.zmin() ) / 20.;
bbox = Scene_interface::Bbox(
    bbox.xmin() - delta_x,
    bbox.ymin() - delta_y,
    bbox.zmin() - delta_z,
    bbox.xmax() + delta_x,
    bbox.ymax() + delta_y,
    bbox.zmax() + delta_z);
  }

  Polyhedron p;
  
  Build_bbox_mesh b(bbox);
  p.delegate(b);

  Scene_item* item = new Scene_polyhedron_item(p);
  item->setName("Scene bbox mesh");
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
}

#include "Create_bbox_mesh_plugin.moc"
