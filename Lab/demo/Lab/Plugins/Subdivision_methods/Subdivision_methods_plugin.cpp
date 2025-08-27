#include <QElapsedTimer>
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QInputDialog>

#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include "Messages_interface.h"
#include "Scene_surface_mesh_item.h"
#include "SMesh_type.h"
#include <CGAL/subdivision_method_3.h>

using namespace CGAL::Three;
namespace params = CGAL::parameters;

class CGAL_Lab_subdivision_methods_plugin :
  public QObject,
  public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")
public:
  // used by CGAL_Lab_plugin_helper
  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface* m)
  {
      this->mw = mw;
      messages = m;
      scene = scene_interface;
      QAction *actionLoop = new QAction("Loop", mw);
      QAction *actionUpsample = new QAction("Upsample", mw);
      QAction *actionCatmullClark = new QAction("Catmull Clark", mw);
      QAction *actionSqrt3 = new QAction("Sqrt3", mw);
      QAction *actionDooSabin = new QAction("Doo Sabin", mw);
      actionLoop->setObjectName("actionLoop");
      actionUpsample->setObjectName("actionUpsample");
      actionCatmullClark->setObjectName("actionCatmullClark");
      actionSqrt3->setObjectName("actionSqrt3");
      actionDooSabin->setObjectName("actionDooSabin");
      _actions << actionLoop
               << actionUpsample
               << actionCatmullClark
               << actionSqrt3
               << actionDooSabin;
      for(QAction* action : _actions)
        action->setProperty("subMenuName", "3D Surface Subdivision Methods");
      autoConnectActions();

  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionLoop_triggered();
  void on_actionUpsample_triggered();
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();
  void on_actionDooSabin_triggered();
private :
  Messages_interface* messages;
  QList<QAction*> _actions;
  template<class FaceGraphItem>
  void apply_loop(FaceGraphItem* item, int nb_steps);
  template<class FaceGraphItem>
  void apply_upsample(FaceGraphItem* item, int nb_steps);
  template<class FaceGraphItem>
  void apply_catmullclark(FaceGraphItem* item, int nb_steps);
  template<class FaceGraphItem>
  void apply_sqrt3(FaceGraphItem* item, int nb_steps);
  template<class FaceGraphItem>
  void apply_doosabin(FaceGraphItem* item, int nb_steps);
}; // end CGAL_Lab_subdivision_methods_plugin


template<class FaceGraphItem>
void CGAL_Lab_subdivision_methods_plugin::apply_loop(FaceGraphItem* item, int nb_steps)
{
  typename FaceGraphItem::Face_graph* graph = item->face_graph();
  QElapsedTimer time;
  time.start();
  CGAL::Three::Three::information("Loop subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Loop_subdivision(*graph, params::number_of_iterations(nb_steps));
  CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

template<class FaceGraphItem>
void CGAL_Lab_subdivision_methods_plugin::apply_upsample(FaceGraphItem* item, int nb_steps)
{
  typename FaceGraphItem::Face_graph* graph = item->face_graph();
  if(!graph) return;
  QElapsedTimer time;
  time.start();
  CGAL::Three::Three::information("Upsample subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (is_triangle_mesh(*graph))
    CGAL::Subdivision_method_3::Loop_subdivision(*graph, params::number_of_iterations(nb_steps)
                                                                .do_not_modify_geometry(true));
  else
    CGAL::Subdivision_method_3::CatmullClark_subdivision(*graph, params::number_of_iterations(nb_steps)
                                                                        .do_not_modify_geometry(true));

  CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void CGAL_Lab_subdivision_methods_plugin::on_actionLoop_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
    return;

  bool ok = true;
  int nb_steps = QInputDialog::getInt(mw,
                                      QString("Number of Iterations"),
                                      QString("Choose number of iterations"),
                                      1 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                      &ok);
  if(!ok)
    return;

  apply_loop(sm_item, nb_steps);
}

void CGAL_Lab_subdivision_methods_plugin::on_actionUpsample_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
    return;

  bool ok = true;
  int nb_steps = QInputDialog::getInt(mw,
                                      QString("Number of Iterations"),
                                      QString("Choose number of iterations"),
                                      1 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                      &ok);
  if(!ok)
    return;

  apply_upsample(sm_item, nb_steps);
}

template<class FaceGraphItem>
void CGAL_Lab_subdivision_methods_plugin::apply_catmullclark(FaceGraphItem* item, int nb_steps)
{
  typename FaceGraphItem::Face_graph* graph = item->face_graph();
  if(!graph) return;
  QElapsedTimer time;
  time.start();
  CGAL::Three::Three::information("Catmull-Clark subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::CatmullClark_subdivision(*graph, params::number_of_iterations(nb_steps));
  CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void CGAL_Lab_subdivision_methods_plugin::on_actionCatmullClark_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
    return;

  bool ok = true;
  int nb_steps = QInputDialog::getInt(mw,
                                      QString("Number of Iterations"),
                                      QString("Choose number of iterations"),
                                      1 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                      &ok);
  if(!ok)
    return;

  apply_catmullclark(sm_item, nb_steps);
}

template<class FaceGraphItem>
void CGAL_Lab_subdivision_methods_plugin::apply_sqrt3(FaceGraphItem* item, int nb_steps)
{
  typename FaceGraphItem::Face_graph* graph = item->face_graph();
  if(!graph) return;
  QElapsedTimer time;
  time.start();
  CGAL::Three::Three::information("Sqrt-3 subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Sqrt3_subdivision(*graph, params::number_of_iterations(nb_steps));
  CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void CGAL_Lab_subdivision_methods_plugin::on_actionSqrt3_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
    return;

  bool ok = true;
  int nb_steps = QInputDialog::getInt(mw,
                                      QString("Number of Iterations"),
                                      QString("Choose number of iterations"),
                                      1 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                      &ok);
  if(!ok)
    return;

  apply_sqrt3(sm_item, nb_steps);

}

template<class FaceGraphItem>
void CGAL_Lab_subdivision_methods_plugin::apply_doosabin(FaceGraphItem* item, int nb_steps)
{
  typename FaceGraphItem::Face_graph* graph = item->face_graph();
  if(!graph) return;
  QElapsedTimer time;
  time.start();
  CGAL::Three::Three::information("Doo-Sabin subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::DooSabin_subdivision(*graph, params::number_of_iterations(nb_steps));
  CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void CGAL_Lab_subdivision_methods_plugin::on_actionDooSabin_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
    return;

  bool ok = true;
  int nb_steps = QInputDialog::getInt(mw,
                                      QString("Number of Iterations"),
                                      QString("Choose number of iterations"),
                                      1 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                      &ok);
  if(!ok)
    return;

  apply_doosabin(sm_item, nb_steps);

}

#include "Subdivision_methods_plugin.moc"
