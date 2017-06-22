#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/gl.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>


using namespace CGAL::Three;

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
  QMainWindow* mw;
  QAction* actionBbox;
  QAction* actionExtendedBbox;

}; // end Create_bbox_mesh_plugin

void Create_bbox_mesh_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
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
  
  Scene_item* item;
  Kernel::Iso_cuboid_3 ic(bbox);
  if(mw->property("is_polyhedorn_mode").toBool()){
    Polyhedron* p = new Polyhedron;
    CGAL::make_hexahedron(ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], ic[6], ic[7],*p);
    item = new Scene_polyhedron_item(p);
  } else {
    SMesh* p = new SMesh;
    CGAL::make_hexahedron(ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], ic[6], ic[7],*p);

    item = new Scene_surface_mesh_item(p);
  }

  item->setName("Scene bbox mesh");
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
}

#include "Create_bbox_mesh_plugin.moc"
