#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Polyhedron_demo_plugin_interface.h"
#include "Messages_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>

#include "CGAL/compute_normal.h"
      // Vector n = compute_facet_normal<Facet,Kernel>(*f);
      // ::glNormal3d(n.x(),n.y(),n.z());

typedef CGAL::Triangulation_2_filtered_projection_traits_3<Kernel>   Traits;

typedef CGAL::Triangulation_vertex_base_with_info_2<Polyhedron::Halfedge_handle,
                                                    Traits>      Vb;

typedef CGAL::Constrained_triangulation_face_base_2<Traits>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::No_intersection_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS,
                                                   Itag>         CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>          CDT;

class Polyhedron_demo_triangulate_facets_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;
    actionTriangulationFacets = new QAction("Triangulate facets", mw);
    if(actionTriangulationFacets) {
      connect(actionTriangulationFacets, SIGNAL(triggered()),
              this, SLOT(triangulate())); 
    }
  };

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTriangulationFacets;
  }

public slots:
  void triangulate() {
    CGAL::set_error_behaviour(CGAL::ABORT);
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
    Scene_polyhedron_item* item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Polyhedron* pMesh = item->polyhedron();
      if(!pMesh) return;
      if(pMesh->is_pure_triangle()) {
        messages->warning(tr("The polyhedron \"%1\" is already triangulated.")
                          .arg(item->name()));
        return;
      }

      typedef Polyhedron::Facet Facet;
      typedef Polyhedron::Facet_iterator Facet_iterator;
      for(Facet_iterator 
            fit = pMesh->facets_begin(),
            end = pMesh->facets_end();
          fit != end; ++fit)
      {
        Kernel::Vector_3 normal = compute_facet_normal<Facet,Kernel>(*fit);

        Traits cdt_traits(normal);
        CDT cdt(cdt_traits);

        Facet::Halfedge_around_facet_circulator 
          he_circ = fit->facet_begin(),
          he_circ_end(he_circ);
        CDT::Vertex_handle previous, first;
        do {
          CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
          if(first == 0) {
            first = vh;
          }
          vh->info() = he_circ;
          if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
          }
          previous = vh;
        } while( ++he_circ != he_circ_end );
        cdt.insert_constraint(previous, first);

        std::cerr << cdt.number_of_vertices() << std::endl;
        for(CDT::Finite_edges_iterator
              eit = cdt.finite_edges_begin(),
              end = cdt.finite_edges_end();
            eit != end; ++eit)
        {
          if(!cdt.is_constrained(*eit)) {
            const CDT::Vertex_handle va = eit->first->vertex(cdt. cw(eit->second));
            const CDT::Vertex_handle vb = eit->first->vertex(cdt.ccw(eit->second));
            std::cerr << va->point() << " " 
                      << (void*)&*va->info() << " "
                      << (void*)&*va->info()->next() << "\n";
            std::cerr << vb->point() << " " 
                      << (void*)&*vb->info() << " "
                      << (void*)&*vb->info()->next() << "\n";
            pMesh->split_facet(va->info(),vb->info());
          }
        }
      }
      scene->itemChanged(item);
    }
  }
  
private:
  QAction* actionTriangulationFacets;
  Messages_interface* messages;
};

Q_EXPORT_PLUGIN2(Polyhedron_demo_triangulate_facets_plugin, Polyhedron_demo_triangulate_facets_plugin);

#include "Polyhedron_demo_triangulate_facets_plugin.moc"
