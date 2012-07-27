#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Polyhedron_demo_plugin_interface.h"
#include "Messages_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include "CGAL/compute_normal.h"
#include "CGAL/triangulate_polyhedron.h"

class Polyhedron_demo_triangulate_facets_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
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

  bool applicable() const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }


public slots:
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
        std::cerr << (void*)&*eit << std::endl;
        Polyhedron::Edge_iterator eit_copy = eit++;
        if(!eit_copy->is_border()) {
          Polyhedron::Facet_handle fh1 = eit_copy->facet();
          Polyhedron::Facet_handle fh2 = eit_copy->opposite()->facet();
          typedef Polyhedron::Facet Facet;
          if( fh1 != fh2 &&  
              !eit_copy->vertex()->is_bivalent() && 
              !eit_copy->opposite()->vertex()->is_bivalent())
          {
            Kernel::Vector_3 v1 = compute_facet_normal<Facet, Kernel>(*fh1);
            Kernel::Vector_3 v2 = compute_facet_normal<Facet, Kernel>(*fh2);
            if(v1 * v2 > 0.99) {
              std::cerr << "join\n";
              // pMesh->is_valid(true);
              pMesh->join_facet(eit_copy);
            }
          }
        }
      }
      CGAL_assertion_code(pMesh->normalize_border());
      // CGAL_assertion(pMesh->is_valid(true, 3));
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

      CGAL::triangulate_polyhedron(*pMesh);

      CGAL_assertion_code(pMesh->normalize_border());
      CGAL_assertion(pMesh->is_valid(false, 3));

      scene->itemChanged(item);

    } // end of if(item)

    } // end of the loop on the selected items
    // default cursor
    QApplication::restoreOverrideCursor();
  }
  
private:
  QAction* actionTriangulateFacets;
  QAction* actionUnTriangulateFacets;  
  Messages_interface* messages;
};

Q_EXPORT_PLUGIN2(Polyhedron_demo_triangulate_facets_plugin, Polyhedron_demo_triangulate_facets_plugin)

#include "Polyhedron_demo_triangulate_facets_plugin.moc"
