#include "Scene_surface_mesh_item.h"
#include "Scene_nef_polyhedron_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QElapsedTimer>
#include <QMessageBox>
using namespace CGAL::Three;
class Polyhedron_demo_nef_plugin :
  public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

void init(QMainWindow* mw,
          Scene_interface* scene_interface,
          Messages_interface* )
{
    scene = scene_interface;
    this->mw = mw;
    QAction *actionConvexDecomposition = new QAction("Convex Decomposition", mw);
    actionConvexDecomposition->setProperty("subMenuName",
                                           "Convex Decomposition");

    actionToNef = new QAction("Convert to Nef Polyhedron", mw);
    actionToPoly = new QAction("Convert to Normal FaceGraph", mw);
    QAction *actionUnion = new QAction("Union (A/B)", mw);
    QAction *actionIntersection = new QAction("Intersection (A/B)", mw);
    QAction *actionDifference = new QAction("&Difference (A/B)", mw);
    QAction *actionMinkowskiSum = new QAction("&Minkowski Sum (A/B)", mw);

    connect(actionToNef       , SIGNAL(triggered()), this, SLOT(on_actionToNef_triggered()));
    connect(actionToPoly      , SIGNAL(triggered()), this, SLOT(on_actionToPoly_triggered()));
    connect(actionUnion       , SIGNAL(triggered()), this, SLOT(on_actionUnion_triggered()));
    connect(actionIntersection, SIGNAL(triggered()), this, SLOT(on_actionIntersection_triggered()));
    connect(actionDifference  , SIGNAL(triggered()), this, SLOT(on_actionDifference_triggered()));
    connect(actionMinkowskiSum, SIGNAL(triggered()), this, SLOT(on_actionMinkowskiSum_triggered()));


    _actions<< actionToNef
            << actionToPoly
            << actionUnion
            << actionIntersection
            << actionDifference
            << actionMinkowskiSum;
    Q_FOREACH(QAction * action, _actions)
      action->setProperty("subMenuName",
                          "Boolean Operations");
    _actions<< actionConvexDecomposition;
    connect(actionConvexDecomposition, SIGNAL(triggered()), this, SLOT(on_actionConvexDecomposition_triggered()));

}

bool applicable(QAction* action) const {
  const int indexA = scene->selectionAindex();
  const int indexB = scene->selectionBindex();

  if( action == actionToNef)
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

  else if( action == actionToPoly)
    return qobject_cast<Scene_nef_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));

  return qobject_cast<Scene_nef_polyhedron_item*>(scene->item(indexA))
      || qobject_cast<Scene_nef_polyhedron_item*>(scene->item(indexB));
}

QList<QAction*> actions() const {
  return _actions;
}

private:
  Scene_interface* scene;
  QAction *actionToNef;
  QAction *actionToPoly;
  QList<QAction*> _actions;
  QMainWindow* mw;

  enum  Boolean_operation { BOOLEAN_UNION,
                            BOOLEAN_INTERSECTION,
                            BOOLEAN_DIFFERENCE,
                            MINKOWSKI_SUM
  };
  void boolean_operation(const Boolean_operation operation);

public Q_SLOTS:
  void on_actionToNef_triggered();
  void on_actionToPoly_triggered();
  void on_actionUnion_triggered();
  void on_actionIntersection_triggered();
  void on_actionDifference_triggered();
  void on_actionMinkowskiSum_triggered();
  void on_actionConvexDecomposition_triggered();
}; // end class Polyhedron_demo_nef_plugin

void
Polyhedron_demo_nef_plugin::on_actionToNef_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* sm_item =
      qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
  {
    return;
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);
  QElapsedTimer time;
  time.start();
  std::cerr << "Convert facegraph to nef polyhedron...";

  Scene_nef_polyhedron_item* new_nef_item;
  new_nef_item = Scene_nef_polyhedron_item::from_polygon_mesh(sm_item);
  new_nef_item->setName(tr("%1 (Nef_polyhedron_3)").arg(sm_item->name()));
  new_nef_item->setRenderingMode(sm_item->renderingMode());
  sm_item->setVisible(false);
  scene->itemChanged(index);
  new_nef_item->invalidateOpenGLBuffers();
  scene->addItem(new_nef_item);
  std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
}


void
Polyhedron_demo_nef_plugin::on_actionConvexDecomposition_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* smitem =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  Scene_nef_polyhedron_item* item =
    (smitem)? Scene_nef_polyhedron_item::from_polygon_mesh(smitem)
           : qobject_cast<Scene_nef_polyhedron_item*>(scene->item(index));
  QApplication::restoreOverrideCursor();
  if(item) {
    QElapsedTimer time;
    time.start();
    std::cerr << "Convex decomposition...";

    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::list<Scene_surface_mesh_item*> convex_parts;
    item->convex_decomposition(convex_parts);
    int i = 0;
    for(std::list<Scene_surface_mesh_item*>::iterator it = convex_parts.begin();
        it != convex_parts.end();
        ++it){
      (*it)->setName(tr("part %1 of %2").arg(i++).arg(item->name()));
      (*it)->setRenderingMode(item->renderingMode());
      scene->addItem(*it);
    }

    if(smitem){
      delete item;
      smitem->setVisible(false);
    } else {
      item->setVisible(false);
    }

    std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
    QApplication::restoreOverrideCursor();
  } else {
    std::cerr << "Only a Facegraph or a Nef Polyhedron can be decomposed in convex parts" << std::endl;
  }
}

void
Polyhedron_demo_nef_plugin::on_actionToPoly_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_nef_polyhedron_item* item =
    qobject_cast<Scene_nef_polyhedron_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QElapsedTimer time;
    time.start();
    std::cerr << "Convert nef polyhedron to facegraph...";

    if(!item->is_simple())
    {
      QMessageBox::warning(mw,
                           tr("Cannot convert"),
                           tr("The nef polyhedron \"%1\" is not simple, "
                              "and thus cannot be converted!")
                              .arg(item->name()));
      QApplication::restoreOverrideCursor();
      return;
    }
        Scene_surface_mesh_item* new_item = item->convert_to_surface_mesh();
        new_item->setName(tr("%1 (from nef)").arg(item->name()));
        new_item->setRenderingMode(item->renderingMode());
        item->setVisible(false);
        scene->itemChanged(index);
        scene->addItem(new_item);
        std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
        QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_nef_plugin::on_actionUnion_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  boolean_operation(BOOLEAN_UNION);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_nef_plugin::on_actionIntersection_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  boolean_operation(BOOLEAN_INTERSECTION);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_nef_plugin::on_actionDifference_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  boolean_operation(BOOLEAN_DIFFERENCE);
  QApplication::restoreOverrideCursor();
}

void
Polyhedron_demo_nef_plugin::on_actionMinkowskiSum_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  boolean_operation(MINKOWSKI_SUM);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_nef_plugin::boolean_operation(const Boolean_operation operation)
{
  const int indexA = scene->selectionAindex();
  const int indexB = scene->selectionBindex();

  if(indexA < 0 || indexB < 0) return;
  if(indexA == indexB) return;

  if(qobject_cast<Scene_surface_mesh_item*>(scene->item(indexA)) ||
     qobject_cast<Scene_surface_mesh_item*>(scene->item(indexB))) {
    QMenu* menu = mw->findChild<QMenu*>("menu_Boolean_operations");
    if(!menu) qWarning("Do not find object named \"menu_Boolean_operations\"!");
    QMessageBox::warning(mw,
                         tr("Boolean operation cannot be applied on normal facegraphs"),
                         tr("You need to call the operation \"%1\" in the menu \"%2\".")
                         .arg(QString("Convert to Nef Polyhedron"))
                         .arg(menu ? menu->title() : "Boolean Operations"));
  }
  Scene_nef_polyhedron_item* itemA =
    qobject_cast<Scene_nef_polyhedron_item*>(scene->item(indexA));
  Scene_nef_polyhedron_item* itemB =
    qobject_cast<Scene_nef_polyhedron_item*>(scene->item(indexB));
  if(!itemA || !itemB)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  // copy itemA
  Scene_nef_polyhedron_item* new_item = 0;
  if(operation != MINKOWSKI_SUM) {
    new_item = new Scene_nef_polyhedron_item(*itemA->nef_polyhedron());
  };

 // perform Boolean operation
  std::cout << "Boolean operation...";
  QElapsedTimer time;
  time.start();
  switch(operation)
  {
  case BOOLEAN_UNION:
    (*new_item) += (*itemB);
    break;
  case BOOLEAN_INTERSECTION:
    (*new_item) *= (*itemB);
    break;
  case BOOLEAN_DIFFERENCE:
    (*new_item) -= (*itemB);
    break;
  case MINKOWSKI_SUM:
    new_item = Scene_nef_polyhedron_item::sum(*itemA,
                                              *itemB);
  }
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  QString name;
  switch(operation)
  {
  case BOOLEAN_UNION:
    name = tr("%1 union %2");
    break;
  case BOOLEAN_INTERSECTION:
    name = tr("%1 intersection %2");
    break;
  case BOOLEAN_DIFFERENCE:
    name = tr("%1 minus %2");
    break;
  case MINKOWSKI_SUM:
    name = tr("Minkowski sum of %1 and %2");
  }

  new_item->setName(name.arg(itemA->name(), itemB->name()));
  new_item->setColor(Qt::green);
  new_item->setRenderingMode(FlatPlusEdges);
  itemA->setRenderingMode(Wireframe);
  itemB->setRenderingMode(Wireframe);
  scene->addItem(new_item);
  scene->itemChanged(indexA);
  scene->itemChanged(indexB);

  QApplication::restoreOverrideCursor();
}

#include "Nef_plugin.moc"
