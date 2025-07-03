#include <QtCore/qglobal.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Variational_medial_axis_skeleton_item.h"
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
class Variational_medial_axis_skeleton_plugin :
    public QObject,
    public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTraceSkeleton;
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

  void trace_skeleton();
  void enableAction();
  void connectNewViewer(QObject* o)
  {
    for(int i=0; i<scene->numberOfEntries(); ++i)
    {
      Variational_medial_axis_skeleton_item* item = qobject_cast<Variational_medial_axis_skeleton_item*>(
            scene->item(i));
      if(item)
        o->installEventFilter(item);
    }
  }

private:
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionTraceSkeleton;
}; // end Variational_medial_axis_skeleton_plugin

void Variational_medial_axis_skeleton_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionTraceSkeleton = new QAction(tr("Variational Medial Axis Skeleton"), mainWindow);
  connect(actionTraceSkeleton, SIGNAL(triggered()),
          this, SLOT(trace_skeleton()));
  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
}

void Variational_medial_axis_skeleton_plugin::trace_skeleton()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

  //todo: add dialog

  Variational_medial_axis_skeleton_item* item = new Variational_medial_axis_skeleton_item(scene, sm_item, 200);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName(sm_item->name() + " medial axis");
  item->setRenderingMode(FlatPlusEdges);
  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(item);

  scene->addItem(item);
  item->fill_subitems();

  sm_item->setVisible(false);

  QApplication::restoreOverrideCursor();
}


void Variational_medial_axis_skeleton_plugin::enableAction() {
  actionTraceSkeleton->setEnabled(true);
}

#include "Variational_medial_axis_skeleton_plugin.moc"
