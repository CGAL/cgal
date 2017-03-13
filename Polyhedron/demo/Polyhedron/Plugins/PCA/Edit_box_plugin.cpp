#include <QtCore/qglobal.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Scene_edit_box_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/boost/graph/helpers.h>
#include <QAction>
#include <QMainWindow>
#include <QApplication>


#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
using namespace CGAL::Three;
class Edit_box_plugin :
    public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionBbox
                             << actionExport;
  }

  bool applicable(QAction* a) const {
    if(a==actionBbox &&scene->numberOfEntries() > 0)
      return true;
    else if(a==actionExport )
    {
      for(int i = 0, end = scene->numberOfEntries();
          i < end; ++i)
      {
        if(qobject_cast<Scene_edit_box_item*>(scene->item(i)))
        {
          return true;
        }
      }
    }
    return false;}
public Q_SLOTS:

  void bbox();
  void enableAction();
  void exportToPoly();

private:
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionBbox;
  QAction* actionExport;


}; // end Edit_box_plugin

void Edit_box_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionBbox = new QAction(tr("Create Bbox"), mainWindow);
  connect(actionBbox, SIGNAL(triggered()),
          this, SLOT(bbox()));
  actionExport = new QAction(tr("Export to Face_graph item"), mainWindow);
  connect(actionExport, SIGNAL(triggered()),
          this, SLOT(exportToPoly()));
}

void Edit_box_plugin::bbox()
{
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(qobject_cast<Scene_edit_box_item*>(scene->item(i)))
      return;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);
  Scene_edit_box_item* item = new Scene_edit_box_item(scene);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName("Edit box");
  item->setRenderingMode(FlatPlusEdges);
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(item);
  scene->addItem(item);
  actionBbox->setEnabled(false);

  QApplication::restoreOverrideCursor();
}

void Edit_box_plugin::enableAction() {
  actionBbox->setEnabled(true);
}

void Edit_box_plugin::exportToPoly()
{
  int id =0;
  const qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
  Kernel::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);
  Scene_edit_box_item* item = NULL;
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    item = qobject_cast<Scene_edit_box_item*>(scene->item(i));
    if(item)
    {
      id = i;
      break;
    }
  }

  Polyhedron::Point_3 points[8];
  for(int i=0; i<8; ++i)
  {
    points[i] = Polyhedron::Point_3(item->point(i,0),item->point(i,2), item->point(i,1))-offset;
  }

  if(mw->property("is_polyhedron_mode").toBool()){
    Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
    CGAL::make_hexahedron(points[0],
                          points[1],
                          points[2],
                          points[3],
                          points[4],
                          points[5],
                          points[6],
                          points[7],
                          *poly_item->polyhedron());
    CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
    item->setName("Edit box");
    item->setRenderingMode(FlatPlusEdges);
    scene->replaceItem(id, poly_item, true);
    item->deleteLater();
    actionBbox->setEnabled(true);
  }else{
 Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item();
    CGAL::make_hexahedron(points[0],
                          points[1],
                          points[2],
                          points[3],
                          points[4],
                          points[5],
                          points[6],
                          points[7],
                          *poly_item->polyhedron());
    CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
    item->setName("Edit box");
    item->setRenderingMode(FlatPlusEdges);
    scene->replaceItem(id, poly_item, true);
    item->deleteLater();
    actionBbox->setEnabled(true);
  }
}
#include "Edit_box_plugin.moc"
