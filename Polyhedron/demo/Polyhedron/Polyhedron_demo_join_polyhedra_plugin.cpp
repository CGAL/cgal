#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Polyhedron_copy_3.h>


class Polyhedron_demo_join_polyhedra_plugin:
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

  #if QT_VERSION >= 0x050000
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")//New for Qt5 version !
  #endif

  QAction* actionJoinPolyhedra;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionJoinPolyhedra; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* /* m */)
  {
    actionJoinPolyhedra= new QAction(tr("Join selected polyhedra"), mainWindow);
    actionJoinPolyhedra->setObjectName("actionJoinPolyhedra");
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    int nb_poly=0;
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) )
      {
        if (nb_poly==1) return true;
        ++nb_poly;
      }
    }
    return false;
  }

public slots:
  void on_actionJoinPolyhedra_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin

void Polyhedron_demo_join_polyhedra_plugin::on_actionJoinPolyhedra_triggered()
{
  Scene_interface::Item_id mainSelectionIndex = -1;
  Scene_polyhedron_item* mainSelectionItem = NULL;


  QList<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices()) {
    if (mainSelectionIndex==-1){
      mainSelectionItem =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
      std::cout << mainSelectionItem << std::endl;
      if(mainSelectionItem!=NULL) mainSelectionIndex=index;
      continue;
    }

    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item)
    {
      indices_to_remove.push_front(index);
      CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron::HalfedgeDS, false>
        modifier( *(item->polyhedron()) );
      mainSelectionItem->polyhedron()->delegate(modifier);
    }
  }

  scene->itemChanged(mainSelectionIndex);

  //remove the other items
  Q_FOREACH(int index, indices_to_remove)
  {
    scene->erase(index);
  }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(Polyhedron_demo_join_polyhedra_plugin, Polyhedron_demo_join_polyhedra_plugin)
#endif

#include "Polyhedron_demo_join_polyhedra_plugin.moc"
