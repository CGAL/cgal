#include "Scene_polyhedron_item.h"
#include "Scene_nef_polyhedron_item.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

class Polyhedron_demo_nef_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  QStringList actionsNames() const {
    return QStringList() << "actionToNef" 
                         << "actionToPoly"
                         << "actionUnion"
                         << "actionIntersection"
                         << "actionDifference"
                         << "actionConvexDecomposition"
                         << "actionMinkowskiSum";
  }

  bool applicable() const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
      || qobject_cast<Scene_nef_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

//   QList<QAction*> actions() const {
//     QMenu* menu = new QMenu(tr("Boolean operations"), mw);
//     QAction* action_to_nef = new QAction(tr("Convert to nef polyhedron"), mw);
//     connect(action_to_nef, SIGNAL(triggered()),
//             this, SLOT(convert_to_nef_polyhedron()));
//     QAction* action_to_poly = new QAction(tr("Convert nef polyhedron to polyhedron"), mw);
//     connect(action_to_poly, SIGNAL(triggered()),
//             this, SLOT(convert_to_polyhedron()));
//     menu->addAction(action_to_nef);
//     menu->addAction(action_to_poly);

//     return QList<QAction*>() << menu->menuAction();
//   }

private:
  enum  Boolean_operation { BOOLEAN_UNION,
                            BOOLEAN_INTERSECTION,
                            BOOLEAN_DIFFERENCE,
                            MINKOWSKI_SUM
  };
  void boolean_operation(const Boolean_operation operation);

public slots:
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
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    std::cerr << "Convert polyhedron to nef polyhedron...";

    Scene_nef_polyhedron_item* new_nef_item = 
      Scene_nef_polyhedron_item::from_polyhedron(item);

    new_nef_item->setName(tr("%1 (to nef)").arg(item->name()));
    new_nef_item->setRenderingMode(item->renderingMode());
    item->setVisible(false);
    scene->itemChanged(index);
    scene->addItem(new_nef_item);
    std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
    QApplication::restoreOverrideCursor();
  }
}


void
Polyhedron_demo_nef_plugin::on_actionConvexDecomposition_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* pitem = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_nef_polyhedron_item* item =   
    (pitem)? Scene_nef_polyhedron_item::from_polyhedron(pitem)
           : qobject_cast<Scene_nef_polyhedron_item*>(scene->item(index));
  
  if(item) {
    QTime time;
    time.start();
    std::cerr << "Convex decomposition...";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    std::list<Scene_polyhedron_item*> convex_parts;
    item->convex_decomposition(convex_parts);
    int i = 0;
    for(std::list<Scene_polyhedron_item*>::iterator it = convex_parts.begin();
        it != convex_parts.end();
        ++it){
      (*it)->setName(tr("part %1 of %2").arg(i++).arg(item->name()));
      (*it)->setRenderingMode(item->renderingMode());
      scene->addItem(*it);
    }

    if(pitem){
      delete item;
      pitem->setVisible(false);
    } else {
      item->setVisible(false);
    }

    std::cerr << "ok (" << time.elapsed() << " ms)" << std::endl;
    QApplication::restoreOverrideCursor();
  } else {
    std::cerr << "Only a Polyhedron or a Nef Polyhedron can be decomposed in convex parts" << std::endl; 
  }
}

void
Polyhedron_demo_nef_plugin::on_actionToPoly_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_nef_polyhedron_item* item = 
    qobject_cast<Scene_nef_polyhedron_item*>(scene->item(index));

  if(item)
  {
    QTime time;
    time.start();
    std::cerr << "Convert nef polyhedron to polyhedron...";

    if(!item->is_simple())
    {
      QMessageBox::warning(mw,
                           tr("Cannot convert"),
                           tr("The nef polyhedron \"%1\" is not simple, "
                              "and thus cannot be converted!")
                              .arg(item->name()));
      return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    Scene_polyhedron_item* new_item = item->convert_to_polyhedron();
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
  boolean_operation(BOOLEAN_UNION);
}

void Polyhedron_demo_nef_plugin::on_actionIntersection_triggered()
{
  boolean_operation(BOOLEAN_INTERSECTION);
}

void Polyhedron_demo_nef_plugin::on_actionDifference_triggered()
{
  boolean_operation(BOOLEAN_DIFFERENCE);
}

void
Polyhedron_demo_nef_plugin::on_actionMinkowskiSum_triggered()
{
  boolean_operation(MINKOWSKI_SUM);
}

void Polyhedron_demo_nef_plugin::boolean_operation(const Boolean_operation operation)
{
  const int indexA = scene->selectionAindex();
  const int indexB = scene->selectionBindex();

  if(indexA < 0 || indexB < 0) return;
  if(indexA == indexB) return;

  if(qobject_cast<Scene_polyhedron_item*>(scene->item(indexA)) ||
     qobject_cast<Scene_polyhedron_item*>(scene->item(indexB))) {
    QMenu* menu = mw->findChild<QMenu*>("menu_Boolean_operations");
    if(!menu) qWarning("Do not find object named \"menu_Boolean_operations\"!");
    QMessageBox::warning(mw,
                         tr("Boolean operation cannot be applied on normal polyhedron"),
                         tr("You need to call the operation \"%1\" in the menu \"%2\".")
                         .arg(actions_map["actionToNef"]->text())
                         .arg(menu ? menu->title() : "Boolean operations"));
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
  QTime time;
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

Q_EXPORT_PLUGIN2(Polyhedron_demo_nef_plugin, Polyhedron_demo_nef_plugin)

#include "Polyhedron_demo_nef_plugin.moc"
