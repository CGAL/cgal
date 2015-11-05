#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>

class Polyhedron_demo_join_and_split_polyhedra_plugin:
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionJoinPolyhedra, *actionSplitPolyhedra, *actionColorConnectedComponents;
  Messages_interface* msg_interface;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionJoinPolyhedra << actionSplitPolyhedra << actionColorConnectedComponents; }
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m)
  {
    msg_interface = m;
    actionJoinPolyhedra= new QAction(tr("Join selected polyhedra"), mainWindow);
    actionJoinPolyhedra->setObjectName("actionJoinPolyhedra");
    actionSplitPolyhedra= new QAction(tr("Split selected polyhedra"), mainWindow);
    actionSplitPolyhedra->setObjectName("actionSplitPolyhedra");
    actionColorConnectedComponents = new QAction(tr("Color each connected component of selected polyhedra"), mainWindow);
    actionColorConnectedComponents->setObjectName("actionColorConnectedComponents");
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable(QAction*) const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) )
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void on_actionJoinPolyhedra_triggered();
  void on_actionSplitPolyhedra_triggered();
  void on_actionColorConnectedComponents_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionJoinPolyhedra_triggered()
{
  Scene_interface::Item_id mainSelectionIndex = -1;
  Scene_polyhedron_item* mainSelectionItem = NULL;


  QList<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices()) {
    if (mainSelectionIndex==-1){
      mainSelectionItem =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
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

  mainSelectionItem->invalidate_buffers();
  scene->itemChanged(mainSelectionIndex);

  //remove the other items
  Q_FOREACH(int index, indices_to_remove)
  {
    scene->erase(index);
  }
}

struct Polyhedron_appender{
  Polyhedron_appender(std::list<Polyhedron*>& new_polyhedra):
    m_new_polyhedra(new_polyhedra) {}
  void operator()(const Polyhedron& p){
    m_new_polyhedra.push_back( new Polyhedron(p) );
  }
  std::list<Polyhedron*>& m_new_polyhedra;
};

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionSplitPolyhedra_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()) {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item)
    {
      std::list<Polyhedron*> new_polyhedra;
      CGAL::internal::corefinement::extract_connected_components(
        *item->polyhedron(),
        boost::make_function_output_iterator(Polyhedron_appender(new_polyhedra))
      );

      if (new_polyhedra.size()==1)
      {
        delete new_polyhedra.front();
        msg_interface->information( tr("%1 has only one connected component").arg(item->name()) );
        continue;
      }

      int cc=0;
      BOOST_FOREACH(Polyhedron* polyhedron_ptr, new_polyhedra)
      {
        Scene_polyhedron_item* new_item=new Scene_polyhedron_item(polyhedron_ptr);
        new_item->setName(tr("%1 - CC %2").arg(item->name()).arg(cc));
        ++cc;
        scene->addItem(new_item);
      }
      item->setVisible(false);
    }
  }
}

struct Polyhedron_cc_marker{
  int cc_index;
  Polyhedron_cc_marker() : cc_index(0) {}
  void start_new_connected_component(){
    ++cc_index;
  }

  template <class Facet_iterator>
  void mark(Facet_iterator begin, Facet_iterator end)
  {
    for(;begin!=end; ++begin)
      (*begin)->set_patch_id(cc_index-1);
  }
};

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionColorConnectedComponents_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()) {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item)
    {
      std::list<Polyhedron*> new_polyhedra;
      Polyhedron_cc_marker marker;
      CGAL::internal::corefinement::mark_connected_components(
        *item->polyhedron(),
        CGAL::internal::corefinement::Dummy_true(),
        marker
      );
      item->invalidate_buffers();
      scene->itemChanged(item);
    }
  }
}


#include "Polyhedron_demo_join_and_split_polyhedra_plugin.moc"
