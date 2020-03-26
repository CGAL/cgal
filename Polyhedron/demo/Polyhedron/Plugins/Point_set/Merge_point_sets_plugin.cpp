#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>


using namespace CGAL::Three;
class Polyhedron_demo_merge_point_sets_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionMergePointSets;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    actionMergePointSets = new QAction(tr("Merge"), mainWindow);
    actionMergePointSets->setObjectName("actionMergePointSets");
    connect(actionMergePointSets, SIGNAL(triggered()), this, SLOT(on_actionMergePointSets_triggered()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMergePointSets;
  }

  bool applicable(QAction*) const {
    if (scene->selectionIndices().size() < 2)
      return false;

    Q_FOREACH(int index, scene->selectionIndices())
      {
        if ( qobject_cast<Scene_points_with_normal_item*>(scene->item(index)) )
          return true;
      }
    return false;
  }

public Q_SLOTS:
  void on_actionMergePointSets_triggered();
private :
  Scene_interface *scene;
}; // end Polyhedron_demo_merge_point_sets_plugin

void Polyhedron_demo_merge_point_sets_plugin::on_actionMergePointSets_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Three::Scene_interface::Item_id mainSelectionIndex
    = scene->selectionIndices().first();
  Scene_points_with_normal_item* mainSelectionItem
    = qobject_cast<Scene_points_with_normal_item*>(scene->item(mainSelectionIndex));

  QList<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices())
    {
      if (index == mainSelectionIndex)
        continue;

      Scene_points_with_normal_item* item =
        qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
      if(item)
        {
          indices_to_remove.push_front(index);
          mainSelectionItem->point_set()->merge_with (*(item->point_set()));
          mainSelectionItem->setName(tr("%1 + %2").arg(mainSelectionItem->name()).arg(item->name()));
        }
    }


  mainSelectionItem->invalidateOpenGLBuffers();
  scene->itemChanged(mainSelectionIndex);

  //remove the other items
  Q_FOREACH(int index, indices_to_remove)
  {
    scene->erase(index);
  }
  QApplication::restoreOverrideCursor();
}


#include "Merge_point_sets_plugin.moc"
