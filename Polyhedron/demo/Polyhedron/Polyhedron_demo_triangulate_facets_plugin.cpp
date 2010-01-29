#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Polyhedron_demo_plugin_interface.h"
#include "Messages_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/HalfedgeDS_decorator.h>

#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>

#include "CGAL/compute_normal.h"

#include <queue>

typedef Polyhedron::Halfedge_handle Halfedge_handle;

typedef CGAL::Triangulation_2_filtered_projection_traits_3<Kernel>     Traits;

typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
                                                    Traits>            Vb;

struct Face_info {
  Halfedge_handle e[3];
  bool is_external;
};

typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
                                                  Traits>              Fb1;

typedef CGAL::Constrained_triangulation_face_base_2<Traits, Fb1>       Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    TDS;
typedef CGAL::No_intersection_tag                                      Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS,Itag>   CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>                CDT;

typedef Polyhedron::HalfedgeDS HDS;

class Triangulate_modifier : public CGAL::Modifier_base<HDS> {
  CDT* cdt;
  Polyhedron::Facet_handle fh;

public:
  Triangulate_modifier(CDT* cdt,
                       Polyhedron::Facet_handle fh) 
    : cdt(cdt), fh(fh)
  {
  }

  bool is_external(CDT::Face_handle fh) const {
    return fh->info().is_external;
  }

  void operator()(HDS& hds) {
    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);
    typedef Polyhedron::Halfedge Halfedge;

    decorator.make_hole(fh->halfedge());
    for(CDT::Finite_edges_iterator
          eit = cdt->finite_edges_begin(),
          end = cdt->finite_edges_end();
        eit != end; ++eit)
    {
      CDT::Face_handle fh = eit->first;
      const int index = eit->second;
      CDT::Face_handle opposite_fh = fh->neighbor(eit->second);
      const int opposite_index = opposite_fh->index(fh);
      const CDT::Vertex_handle va = fh->vertex(cdt-> cw(index));
      const CDT::Vertex_handle vb = fh->vertex(cdt->ccw(index));

      if( ! (is_external(fh) && is_external(opposite_fh)) && 
          ! cdt->is_constrained(*eit) ) 
      {
        // strictly internal edge
        Halfedge_handle h = hds.edges_push_back(Halfedge(),
                                                Halfedge());
        fh->info().e[index] = h;
        opposite_fh->info().e[opposite_index] = h->opposite();

        decorator.set_vertex(h, va->info()->vertex());
        decorator.set_vertex(h->opposite(), vb->info()->vertex());
      }
      if( cdt->is_constrained(*eit) )
      {
        if(!is_external(fh)) {
          fh->info().e[index] = va->info();
        }
        if(!is_external(opposite_fh)) {
          opposite_fh->info().e[opposite_index] = vb->info();
        }
      }
    }
    for(CDT::Finite_faces_iterator
          fit = cdt->finite_faces_begin(),
          end = cdt->finite_faces_end();
        fit != end; ++fit)
    {
      if(!is_external(fit)) 
      {
        Halfedge_handle h0 = fit->info().e[0];
        Halfedge_handle h1 = fit->info().e[1];
        Halfedge_handle h2 = fit->info().e[2];
        CGAL_assertion( h0 != Halfedge_handle() );
        CGAL_assertion( h1 != Halfedge_handle() );
        CGAL_assertion( h2 != Halfedge_handle() );

        typedef Halfedge::Base HBase;
        h0->HBase::set_next(h1);
        decorator.set_prev(h1, h0);
        h1->HBase::set_next(h2);
        decorator.set_prev(h2, h1);
        h2->HBase::set_next(h0);
        decorator.set_prev(h0, h2);

        decorator.fill_hole(h0);        
      }
    }
  }
};

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

public slots:
  void untriangulate() {
    CGAL::set_error_behaviour(CGAL::ABORT);
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

      QApplication::setOverrideCursor(Qt::WaitCursor);

      typedef Polyhedron::Facet Facet;
      typedef Polyhedron::Facet_iterator Facet_iterator;
      typedef Polyhedron::Facet_handle Facet_handle;

      // One need to store facet handles into a vector, because the list of
      // facets of the polyhedron will be modified during the loop, and
      // that invalidates the range [facets_begin(), facets_end()[.
      std::vector<Facet_handle> facets;
      facets.reserve(pMesh->size_of_facets());
      for(Facet_iterator 
            fit = pMesh->facets_begin(),
            end = pMesh->facets_end();
          fit != end; ++fit) {
        facets.push_back(fit);
      }

      // Iterates on the vector of facet handles
      for(std::vector<Facet_handle>::iterator 
            fit_it = facets.begin(),
            end = facets.end();
          fit_it != end; ++fit_it)
      {
        Facet_handle fit = *fit_it;
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

        // sets mark is_external
        for(CDT::Finite_faces_iterator
              fit = cdt.finite_faces_begin(),
              end = cdt.finite_faces_end();
            fit != end; ++fit)
        {
          fit->info().is_external = false;
        }
        std::queue<CDT::Face_handle> face_queue;
        face_queue.push(cdt.infinite_vertex()->face());
        while(! face_queue.empty() ) {
          CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          CGAL_assertion(cdt.is_infinite(fh));
          if(fh->info().is_external) continue;
          std::cerr << (void*)(&*fh) << std::endl;
          fh->info().is_external = true;
          for(int i = 0; i <3; ++i) {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
              face_queue.push(fh->neighbor(i));
            }
          }
        }
        // then modify the polyhedron
        Triangulate_modifier modifier(&cdt, fit);
        pMesh->delegate(modifier);
      }
      CGAL_assertion_code(pMesh->normalize_border());
      CGAL_assertion(pMesh->is_valid(false, 3));
      scene->itemChanged(item);
      // default cursor
      QApplication::restoreOverrideCursor();
    }
  }
  
private:
  QAction* actionTriangulateFacets;
  QAction* actionUnTriangulateFacets;  
  Messages_interface* messages;
};

Q_EXPORT_PLUGIN2(Polyhedron_demo_triangulate_facets_plugin, Polyhedron_demo_triangulate_facets_plugin);

#include "Polyhedron_demo_triangulate_facets_plugin.moc"
