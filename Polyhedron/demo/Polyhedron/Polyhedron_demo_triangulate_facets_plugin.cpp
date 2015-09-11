#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Polyhedron_demo_plugin_interface.h"
#include "Messages_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

class Polyhedron_demo_triangulate_facets_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // To silent a warning -Woverloaded-virtual
  // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
  using Polyhedron_demo_plugin_helper::init;

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;
    actionTriangulateFacets = new QAction("Triangulate facets", mw);
    if(actionTriangulateFacets) {
      connect(actionTriangulateFacets, SIGNAL(triggered()),
              this, SLOT(triangulate())); 
    }
    actionUnTriangulateFacets = new QAction("Untriangulate facets", mw);
    if(actionUnTriangulateFacets) {
      connect(actionUnTriangulateFacets, SIGNAL(triggered()),
              this, SLOT(untriangulate())); 
    }
  };

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTriangulateFacets
                             << actionUnTriangulateFacets;
  }

  bool applicable(QAction*) const { 
    Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())  {
      Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
      if(!item) return false;
    }
    return true;
  }


public Q_SLOTS:
  void untriangulate() {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
    Scene_polyhedron_item* item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Polyhedron* pMesh = item->polyhedron();
      if(!pMesh) return;

      QApplication::setOverrideCursor(Qt::WaitCursor);

      for(Polyhedron::Edge_iterator 
            eit = pMesh->edges_begin(),
            end = pMesh->edges_end();
          eit != end; /*increment is done manually*/)
      {
        // std::cerr << (void*)&*eit << std::endl;
        Polyhedron::Edge_iterator eit_copy = eit++;
        if(!eit_copy->is_border()) {
          Polyhedron::Facet_handle fh1 = eit_copy->facet();
          Polyhedron::Facet_handle fh2 = eit_copy->opposite()->facet();
          if( fh1 != fh2 &&  
              !eit_copy->vertex()->is_bivalent() && 
              !eit_copy->opposite()->vertex()->is_bivalent())
          {
            Kernel::Vector_3 v1 =
              CGAL::Polygon_mesh_processing::compute_face_normal(fh1, *pMesh);
            Kernel::Vector_3 v2 =
              CGAL::Polygon_mesh_processing::compute_face_normal(fh2, *pMesh);
            if(v1 * v2 > 0.99) {
              // std::cerr << "join\n";
              // pMesh->is_valid(true);
              pMesh->join_facet(eit_copy);
            }
          }
        }
      }
      CGAL_assertion_code(pMesh->normalize_border());
      // CGAL_assertion(pMesh->is_valid(true, 3));
      item->invalidate_buffers();
      scene->itemChanged(item);
      // default cursor
      QApplication::restoreOverrideCursor();
    }
  }

  void triangulate() {
    Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())  {

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

      QApplication::setOverrideCursor(Qt::WaitCursor);

      CGAL::Polygon_mesh_processing::triangulate_faces(*pMesh);

      CGAL_assertion_code(pMesh->normalize_border());
      CGAL_assertion(pMesh->is_valid(false, 3));

      item->invalidate_buffers();
      scene->itemChanged(item);
      // default cursor
      QApplication::restoreOverrideCursor();
    } // end of if(item)

    } // end of the loop on the selected items
  }
  
private:
  QAction* actionTriangulateFacets;
  QAction* actionUnTriangulateFacets;  
  Messages_interface* messages;
};

#include "Polyhedron_demo_triangulate_facets_plugin.moc"
