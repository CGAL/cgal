#include <QtCore/qglobal.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Locally_shortest_path_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/boost/graph/helpers.h>
#include <QAction>
#include <QMainWindow>
#include <QApplication>


#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
using namespace CGAL::Three;
class Locally_shortest_path_plugin :
    public QObject,
    public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTracePath << actionTraceBezier;
  }

  bool applicable(QAction*) const {
    if(scene->numberOfEntries() > 0)
    {
      int item_id = scene->mainSelectionIndex();
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(item_id));
    }
    return false;
  }
public Q_SLOTS:

  void trace_path();
  void trace_bezier();
  void enableAction();
  void connectNewViewer(QObject* o)
  {
    for(int i=0; i<scene->numberOfEntries(); ++i)
    {
      Locally_shortest_path_item* item = qobject_cast<Locally_shortest_path_item*>(
            scene->item(i));
      if(item)
        o->installEventFilter(item);
    }
  }

private:
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionTracePath;
  QAction* actionTraceBezier;
}; // end Locally_shortest_path_plugin

void Locally_shortest_path_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionTracePath = new QAction(tr("Create Locally Shortest Path"), mainWindow);
  actionTraceBezier= new QAction(tr("Trace Bezier Curve"), mainWindow);
  connect(actionTracePath, SIGNAL(triggered()),
          this, SLOT(trace_path()));
  connect(actionTraceBezier, SIGNAL(triggered()),
          this, SLOT(trace_bezier()));
  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
}

void Locally_shortest_path_plugin::trace_path()
{
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(qobject_cast<Locally_shortest_path_item*>(scene->item(i)))
      return;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

  Scene_polylines_item* polyline_item = new Scene_polylines_item();

  polyline_item->setName(tr("Locally Shortest Path"));
  polyline_item->setColor(Qt::red);
  scene->addItem(polyline_item);
  polyline_item->invalidateOpenGLBuffers();

  Locally_shortest_path_item* item = new Locally_shortest_path_item(scene, sm_item, polyline_item,2);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName("Source and Target Points");
  item->setRenderingMode(FlatPlusEdges);
  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(item);

  scene->addItem(item);
  actionTracePath->setEnabled(false);
  actionTraceBezier->setEnabled(false);

  QApplication::restoreOverrideCursor();
}

void Locally_shortest_path_plugin::trace_bezier()
{
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(qobject_cast<Locally_shortest_path_item*>(scene->item(i)))
      return;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

  Scene_polylines_item* polyline_item = new Scene_polylines_item();

  polyline_item->setName(tr("Locally Shortest Path"));
  polyline_item->setColor(Qt::red);
  scene->addItem(polyline_item);
  polyline_item->invalidateOpenGLBuffers();

  Locally_shortest_path_item* item = new Locally_shortest_path_item(scene, sm_item, polyline_item,4);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName("Control Points");
  item->setRenderingMode(FlatPlusEdges);
  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(item);

  scene->addItem(item);
  actionTracePath->setEnabled(false);
  actionTraceBezier->setEnabled(false);

  QApplication::restoreOverrideCursor();
}

void Locally_shortest_path_plugin::enableAction() {
  actionTracePath->setEnabled(true);
  actionTraceBezier->setEnabled(true);
}

#include "Locally_shortest_path_plugin.moc"
