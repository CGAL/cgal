#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
using namespace CGAL::Three;
class Polyhedron_demo_triangulate_facets_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;
    actionTriangulateFacets = new QAction("Triangulate Facets", mw);
    actionTriangulateFacets->setProperty("subMenuName","Polygon Mesh Processing");
    if(actionTriangulateFacets) {
      connect(actionTriangulateFacets, SIGNAL(triggered()),
              this, SLOT(triangulate())); 
    }
  };

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTriangulateFacets;
  }

  bool applicable(QAction*) const {
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices()){
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) ||
           qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) )
        return true;
    }
    return false;
  }

public Q_SLOTS:
   void triangulate() {
      QApplication::setOverrideCursor(Qt::WaitCursor);
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())  {

    Scene_polyhedron_item* item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Polyhedron* pMesh = item->polyhedron();
      if(!pMesh) continue;
      if(pMesh->is_pure_triangle()) {
        messages->warning(tr("The polyhedron \"%1\" is already triangulated.")
                          .arg(item->name()));
        continue;
      }


      if(!CGAL::Polygon_mesh_processing::triangulate_faces(*pMesh))
        messages->warning(tr("Some facets could not be triangulated."));

      CGAL_assertion_code(pMesh->normalize_border());
      CGAL_assertion(pMesh->is_valid(false, 3));

      item->invalidateOpenGLBuffers();
      scene->itemChanged(item);
   
    } else {
      Scene_surface_mesh_item* sm_item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      SMesh* pMesh = sm_item->polyhedron();
      if(!pMesh) continue;
      if(is_triangle_mesh(*pMesh)) {
      messages->warning(tr("The polyhedron  \"%1\"  is already triangulated.")
                          .arg(sm_item->name()) );
        continue;
      }
      if(!CGAL::Polygon_mesh_processing::triangulate_faces(*pMesh))
        messages->warning(tr("Some facets could not be triangulated."));

      sm_item->invalidateOpenGLBuffers();
      scene->itemChanged(sm_item);
    } // end of if(item)
    } // end of the loop on the selected items   

    // default cursor
    QApplication::restoreOverrideCursor();
  }
  
private:
  QAction* actionTriangulateFacets;
  Messages_interface* messages;
};

#include "Triangulate_facets_plugin.moc"
