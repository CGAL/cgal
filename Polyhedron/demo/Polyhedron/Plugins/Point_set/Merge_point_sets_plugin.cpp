#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
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
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionMergePointSets;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface) {
    actionMergePointSets = new QAction(tr("Merge"), mainWindow);
    actionMergePointSets->setObjectName("actionMergePointSets");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
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
}; // end Polyhedron_demo_merge_point_sets_plugin

void Polyhedron_demo_merge_point_sets_plugin::on_actionMergePointSets_triggered()
{
  CGAL::Three::Scene_interface::Item_id mainSelectionIndex
    = scene->mainSelectionIndex();
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
          std::copy (item->point_set()->begin(), item->point_set()->end(),
                     std::back_inserter (*(mainSelectionItem->point_set())));
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
}


#include "Merge_point_sets_plugin.moc"
