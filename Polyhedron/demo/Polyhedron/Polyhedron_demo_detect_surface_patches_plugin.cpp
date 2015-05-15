#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_type.h"

#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

#include <set>
#include <algorithm>

class Polyhedron_demo_detect_surface_patches_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionSurfacePatches = new QAction("Detect surface patches", mw);
    actionSurfacePatches->setObjectName("detectSurfacePatches");
    if (actionSurfacePatches )
    {
      connect(actionSurfacePatches, SIGNAL(triggered()),
              this, SLOT(detectSurfacePatches()));
    }
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionSurfacePatches;
  }

  bool applicable(QAction*) const
  {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())))
        return true;
    }
    return false;
  }

public slots:
  void detectSurfacePatches();
  
private:
  QAction* actionSurfacePatches;
}; // end Polyhedron_demo_detect_surface_patches_plugin


void
Polyhedron_demo_detect_surface_patches_plugin::
detectSurfacePatches()
{
  typedef std::pair<int,Polyhedron*> Poly_tuple;
  
  // Get selected items
  QList<Poly_tuple> polyhedrons;
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if ( NULL == item ) { continue; }
    
    Polyhedron* pMesh = item->polyhedron();
    if ( NULL == pMesh ) { continue; }
    
    polyhedrons << std::make_pair(index, pMesh);
  }
  
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
  // Detect surface patches
  CGAL::Mesh_3::Detect_features_in_polyhedra<Polyhedron> detect_features;
  Q_FOREACH(Poly_tuple tuple, polyhedrons)
  {
    Polyhedron* pMesh = tuple.second;
    if ( NULL == pMesh ) { continue; }
    
    // Detect patches in current polyhedron
    detect_features.detect_surface_patches(*pMesh);

    // update scene
    scene->itemChanged(tuple.first);
  }
  
  // default cursor
  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_detect_surface_patches_plugin,
                 Polyhedron_demo_detect_surface_patches_plugin)

#include "Polyhedron_demo_detect_surface_patches_plugin.moc"
