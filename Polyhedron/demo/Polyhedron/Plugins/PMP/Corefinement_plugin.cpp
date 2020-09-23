#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include "Scene_surface_mesh_item.h"

#include <CGAL/Polygon_mesh_processing/corefinement.h>

using namespace CGAL::Three;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

class Polyhedron_demo_corefinement_sm_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "corefinement_plugin.json")

  enum bool_op {CRF_UNION, CRF_INTER, CRF_MINUS, CRF_MINUS_OP};

public:

  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;
    actionCorefine = new QAction("Corefine", mw);
    actionCorefine->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    if(actionCorefine)
      connect(actionCorefine, SIGNAL(triggered()),  this, SLOT(corefine()));

    actionUnion = new QAction("Compute Union", mw);
    actionUnion->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    if(actionUnion)
      connect(actionUnion, SIGNAL(triggered()),  this, SLOT(corefine_and_union()));

    actionInter = new QAction("Compute Intersection", mw);
    actionInter->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    if(actionInter)
      connect(actionInter, SIGNAL(triggered()),  this, SLOT(corefine_and_inter()));

    actionDiff = new QAction("Compute Difference", mw);
    actionDiff->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    if(actionDiff)
      connect(actionDiff, SIGNAL(triggered()),  this, SLOT(corefine_and_diff()));

    actionDiffRev = new QAction("Compute Opposite Difference", mw);
    actionDiffRev->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    if(actionDiffRev)
      connect(actionDiffRev, SIGNAL(triggered()),  this, SLOT(corefine_and_diff_rev()));
  };

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionCorefine
                             << actionUnion
                             << actionInter
                             << actionDiff
                             << actionDiffRev;
  }

  bool applicable(QAction*) const {
    if(scene->selectionIndices().size() != 2)
      return false;
    CGAL::Three::Scene_item* item1 = scene->item(scene->selectionIndices().first());
    CGAL::Three::Scene_item* item2 = scene->item(scene->selectionIndices().last());

    if( qobject_cast<Scene_surface_mesh_item*>(item1))
    {
      if(!qobject_cast<Scene_surface_mesh_item*>(item2))
        return false;
    }
    else
      return false;
    return true;
  }

public Q_SLOTS:
   void corefine() {
     if(scene->selectionIndices().size() != 2)
       return;

     CGAL::Three::Scene_item* item1 = scene->item(scene->selectionIndices().first());
     CGAL::Three::Scene_item* item2 = scene->item(scene->selectionIndices().last());
     if( qobject_cast<Scene_surface_mesh_item*>(item1))
     {
       apply_corefine(qobject_cast<Scene_surface_mesh_item*>(item1),
                      qobject_cast<Scene_surface_mesh_item*>(item2));
     }
  }

  void corefine_and_bool_op(bool_op op)
  {
    if(scene->selectionIndices().size() != 2)
      return;

    CGAL::Three::Scene_item* item1 = scene->item(scene->selectionIndices().first());
    CGAL::Three::Scene_item* item2 = scene->item(scene->selectionIndices().last());
    if( qobject_cast<Scene_surface_mesh_item*>(item1))
    {
      apply_corefine_and_bool_op(qobject_cast<Scene_surface_mesh_item*>(item1),
                     qobject_cast<Scene_surface_mesh_item*>(item2),
                     op);
    }
  }

  void corefine_and_union()
  {
    corefine_and_bool_op(CRF_UNION);
  }
  void corefine_and_inter()
  {
    corefine_and_bool_op(CRF_INTER);
  }
  void corefine_and_diff()
  {
    corefine_and_bool_op(CRF_MINUS);
  }
  void corefine_and_diff_rev()
  {
    corefine_and_bool_op(CRF_MINUS_OP);
  }


private:
  QAction* actionCorefine;
  QAction* actionUnion;
  QAction* actionInter;
  QAction* actionDiff;
  QAction* actionDiffRev;
  Messages_interface* messages;
  template<class Item>
  void apply_corefine(Item* item1, Item* item2)
  {
    if(! CGAL::is_triangle_mesh(*item1->face_graph())) {
      CGAL::Three::Three::warning(tr("The face graph \"%1\" is not triangulated.")
                        .arg(item1->name()));
      return;
    }
    if(! CGAL::is_triangle_mesh(*item2->face_graph())) {
      CGAL::Three::Three::warning(tr("The face graph \"%1\" is not triangulated.")
                        .arg(item2->name()));
      return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    try{
      PMP::corefine(*item1->face_graph(), *item2->face_graph(), params::throw_on_self_intersection(true));
      item1->resetColors();
      item2->resetColors();
      item1->invalidateOpenGLBuffers();
      item2->invalidateOpenGLBuffers();
      scene->itemChanged(item2);
      scene->itemChanged(item1);
    }
    catch(CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception)
    {
      CGAL::Three::Three::warning(tr("The requested operation is not possible due to the presence of self-intersections in the neighborhood of the intersection."));
    }
    // default cursor
    QApplication::restoreOverrideCursor();
  }

  template< class Item>
  void apply_corefine_and_bool_op(Item* first_item, Item* item,bool_op op )
  {
    typedef typename Item::Face_graph FaceGraph;
    if(! CGAL::is_triangle_mesh(*first_item->face_graph())) {
      CGAL::Three::Three::warning(tr("The polyhedron \"%1\" is not triangulated.")
                        .arg(first_item->name()));
      return;
    }
    if(! CGAL::is_triangle_mesh(*item->face_graph())) {
      CGAL::Three::Three::warning(tr("The polyhedron \"%1\" is not triangulated.")
                        .arg(item->name()));
      return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FaceGraph* new_poly = new FaceGraph();
    QString str_op;
    FaceGraph P, Q;
    try{
      switch(op)
      {
        case CRF_UNION:
          P = *first_item->face_graph(), Q = *item->face_graph();
          if (! PMP::corefine_and_compute_union(P, Q, *new_poly, params::throw_on_self_intersection(true)) )
          {
            delete new_poly;
            CGAL::Three::Three::warning(tr("The result of the requested operation is not manifold and has not been computed."));
            // default cursor
            QApplication::restoreOverrideCursor();
            return;
          }
          str_op = "Union";
        break;
        case CRF_INTER:
          P = *first_item->polyhedron(), Q = *item->polyhedron();
          if (! PMP::corefine_and_compute_intersection(P, Q, *new_poly, params::throw_on_self_intersection(true)) )
          {
            delete new_poly;
            CGAL::Three::Three::warning(tr("The result of the requested operation is not manifold and has not been computed."));
            // default cursor
            QApplication::restoreOverrideCursor();
            return;
          }
          str_op = "Intersection";
        break;
        case CRF_MINUS_OP:
          std::swap(first_item, item);
          CGAL_FALLTHROUGH;
        case CRF_MINUS:
          P = *first_item->polyhedron(), Q = *item->polyhedron();
          if (! PMP::corefine_and_compute_difference(P, Q, *new_poly, params::throw_on_self_intersection(true)) )
          {
            delete new_poly;
            CGAL::Three::Three::warning(tr("The result of the requested operation is not manifold and has not been computed."));
            // default cursor
            QApplication::restoreOverrideCursor();
            return;
          }
          str_op = "Difference";
      }
    }
    catch(CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception)
    {
      CGAL::Three::Three::warning(tr("The requested operation is not possible due to the presence of self-intersections in the neighborhood of the intersection."));
      QApplication::restoreOverrideCursor();
    }

    first_item->invalidateOpenGLBuffers();
    item->invalidateOpenGLBuffers();
    scene->itemChanged(item);
    scene->itemChanged(first_item);

    Item* new_item = new Item(new_poly);
    new_item->setName(QString("%1 of %2 and %3").arg(str_op).arg(first_item->name()).arg(item->name()));
    new_item->setColor(first_item->color());
    new_item->setRenderingMode(first_item->renderingMode());
    new_item->setVisible(first_item->visible());
    scene->addItem(new_item);
    new_item->invalidateOpenGLBuffers();

    // default cursor
    QApplication::restoreOverrideCursor();
  }

};

#include "Corefinement_plugin.moc"
