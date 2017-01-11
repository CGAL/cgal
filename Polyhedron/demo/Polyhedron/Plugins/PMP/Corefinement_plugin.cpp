#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/corefinement.h>

using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;
class Polyhedron_demo_corefinement_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

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
    int nb_selected=0;
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) )
        ++nb_selected;
    return nb_selected==2;
  }

public Q_SLOTS:
   void corefine() {
    Scene_polyhedron_item* first_item = NULL;
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
    {
      Scene_polyhedron_item* item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));

      if(item)
      {
        if (first_item==NULL)
        {
          first_item = item;
          continue;
        }

        if(!first_item->polyhedron()->is_pure_triangle()) {
          messages->warning(tr("The polyhedron \"%1\" is not triangulated.")
                            .arg(first_item->name()));
          return;
        }
        if(!item->polyhedron()->is_pure_triangle()) {
          messages->warning(tr("The polyhedron \"%1\" is not triangulated.")
                            .arg(item->name()));
          return;
        }

        QApplication::setOverrideCursor(Qt::WaitCursor);
        PMP::corefine(*first_item->polyhedron(), *item->polyhedron());
        first_item->invalidateOpenGLBuffers();
        item->invalidateOpenGLBuffers();
        scene->itemChanged(item);
        scene->itemChanged(first_item);
        // default cursor
        QApplication::restoreOverrideCursor();
        break;
      } // end of if(item)
    } // end of the loop on the selected items
  }


  void corefine_and_bool_op(bool_op op)
  {
    Scene_polyhedron_item* first_item = NULL;
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
    {
      Scene_polyhedron_item* item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));

      if(item)
      {
        if (first_item==NULL)
        {
          first_item = item;
          continue;
        }

        if(!first_item->polyhedron()->is_pure_triangle()) {
          messages->warning(tr("The polyhedron \"%1\" is not triangulated.")
                            .arg(first_item->name()));
          return;
        }
        if(!item->polyhedron()->is_pure_triangle()) {
          messages->warning(tr("The polyhedron \"%1\" is not triangulated.")
                            .arg(item->name()));
          return;
        }

        QApplication::setOverrideCursor(Qt::WaitCursor);
        Polyhedron* new_poly = new Polyhedron();
        QString str_op;
        Polyhedron P, Q;
        switch(op)
        {
          case CRF_UNION:
            P = *first_item->polyhedron(), Q = *item->polyhedron();
            if (! PMP::corefine_and_compute_union(P, Q, *new_poly) )
            {
              delete new_poly;
              messages->warning(tr("The result of the requested operation is not manifold and has not been computed."));
              return;
            }
            str_op = "Union";
          break;
          case CRF_INTER:
            P = *first_item->polyhedron(), Q = *item->polyhedron();
            if (! PMP::corefine_and_compute_intersection(P, Q, *new_poly) )
            {
              delete new_poly;
              messages->warning(tr("The result of the requested operation is not manifold and has not been computed."));
              return;
            }
            str_op = "Intersection";
          break;
          case CRF_MINUS_OP:
            std::swap(first_item, item);
          case CRF_MINUS:
            P = *first_item->polyhedron(), Q = *item->polyhedron();
            if (! PMP::corefine_and_compute_difference(P, Q, *new_poly) )
            {
              delete new_poly;
              messages->warning(tr("The result of the requested operation is not manifold and has not been computed."));
              return;
            }
            str_op = "Difference";
          break;
        }

        first_item->invalidateOpenGLBuffers();
        item->invalidateOpenGLBuffers();
        scene->itemChanged(item);
        scene->itemChanged(first_item);

        Scene_polyhedron_item* new_item = new Scene_polyhedron_item(new_poly);
        new_item->setName(QString("%1 of %2 and %3").arg(str_op).arg(first_item->name()).arg(item->name()));
        new_item->setColor(first_item->color());
        new_item->setRenderingMode(first_item->renderingMode());
        new_item->setVisible(first_item->visible());
        scene->addItem(new_item);
        new_item->invalidateOpenGLBuffers();

        // default cursor
        QApplication::restoreOverrideCursor();
        break;
      } // end of if(item)
    } // end of the loop on the selected items
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
};

#include "Corefinement_plugin.moc"
