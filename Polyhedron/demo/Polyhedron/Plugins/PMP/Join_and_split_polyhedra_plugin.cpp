#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/internal/corefinement/connected_components.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_map.hpp>
#include "Color_map.h"
using namespace CGAL::Three;
class Polyhedron_demo_join_and_split_polyhedra_plugin:
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  QAction* actionJoinPolyhedra, *actionSplitPolyhedra, *actionColorConnectedComponents;
  Messages_interface* msg_interface;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionJoinPolyhedra << actionSplitPolyhedra << actionColorConnectedComponents; }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m)
  {
    mw = mainWindow;
    scene = scene_interface;
    msg_interface = m;

    actionJoinPolyhedra= new QAction(tr("Join Selected Polyhedra"), mainWindow);
    actionJoinPolyhedra->setProperty("subMenuName", "Operations on Polyhedra");
    actionJoinPolyhedra->setObjectName("actionJoinPolyhedra");

    actionSplitPolyhedra= new QAction(tr("Split Selected Polyhedra"), mainWindow);
    actionSplitPolyhedra->setProperty("subMenuName", "Operations on Polyhedra");
    actionSplitPolyhedra->setObjectName("actionSplitPolyhedra");

    actionColorConnectedComponents = new QAction(tr("Color Each Connected Component"), mainWindow);
    actionColorConnectedComponents ->setProperty("subMenuName", "Polygon Mesh Processing");
    actionColorConnectedComponents->setObjectName("actionColorConnectedComponents");

    autoConnectActions();
  }

  bool applicable(QAction* a) const
  {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
      else if (a == actionColorConnectedComponents
            && qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void on_actionJoinPolyhedra_triggered();
  void on_actionSplitPolyhedra_triggered();
  void on_actionColorConnectedComponents_triggered();

private :
  CGAL::Three::Scene_interface* scene;
  std::vector<QColor> colors_;
}; // end Polyhedron_demo_polyhedron_stitching_plugin

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionJoinPolyhedra_triggered()
{
  CGAL::Three::Scene_interface::Item_id mainSelectionIndex
    = scene->mainSelectionIndex();
  Scene_polyhedron_item* mainSelectionItem
    = qobject_cast<Scene_polyhedron_item*>(scene->item(mainSelectionIndex));
  if(!mainSelectionItem)
  {
    std::cerr<<"No selected polyhedron_item"<<std::endl;
    return;
  }
  QList<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices()) {
    if (index == mainSelectionIndex)
      continue;

    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item)
    {
      indices_to_remove.push_front(index);
      CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron::HalfedgeDS, false>
        modifier( *(item->polyhedron()) );
      mainSelectionItem->polyhedron()->delegate(modifier);
    }
    else
      std::cerr<<"No selected polyhedron_item"<<std::endl;
  }
  mainSelectionItem->invalidateOpenGLBuffers();
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
    colors_.clear();
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item)
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      std::list<Polyhedron*> new_polyhedra;
      CGAL::internal::corefinement::extract_connected_components(
        *item->polyhedron(),
        boost::make_function_output_iterator(Polyhedron_appender(new_polyhedra))
      );

      if (new_polyhedra.size()==1)
      {
        delete new_polyhedra.front();
        msg_interface->information( tr("%1 has only one connected component").arg(item->name()) );
        QApplication::restoreOverrideCursor();
        continue;
      }

      int cc=0;
      compute_color_map(item->color(), item->isItemMulticolor() ? static_cast<unsigned int>(new_polyhedra.size()) : 1,
                        std::back_inserter(colors_));
      Scene_group_item *group = new Scene_group_item("CC");
       scene->addItem(group);
      BOOST_FOREACH(Polyhedron* polyhedron_ptr, new_polyhedra)
      {
        Scene_polyhedron_item* new_item=new Scene_polyhedron_item(polyhedron_ptr);
        new_item->setName(tr("%1 - CC %2").arg(item->name()).arg(cc));
        new_item->setColor(colors_[item->isItemMulticolor()? cc : 0]);
        ++cc;
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      item->setVisible(false);
      QApplication::restoreOverrideCursor();
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
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::set<Scene_polyhedron_item*> to_skip;

  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(item && to_skip.find(item) == to_skip.end())
    {
      item->setItemIsMulticolor(true);
      Polyhedron_cc_marker marker;
      CGAL::internal::corefinement::mark_connected_components(
        *item->polyhedron(),
        CGAL::internal::corefinement::Dummy_true(),
        marker
      );
      item->invalidateOpenGLBuffers();
      scene->itemChanged(item);
    }
    else
    {
      Scene_polyhedron_selection_item* selection_item =
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

      if (selection_item)
      {
        namespace PMP = CGAL::Polygon_mesh_processing;
        typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;

        selection_item->polyhedron_item()->setItemIsMulticolor(true);
        selection_item->polyhedron_item()->set_color_vector_read_only(false);

        const Polyhedron& pmesh = *(selection_item->polyhedron());

        boost::property_map<Polyhedron, boost::face_external_index_t>::type fim
          = get(boost::face_external_index, pmesh);
        boost::vector_property_map<int,
          boost::property_map<Polyhedron, boost::face_external_index_t>::type>
          fccmap(fim);

        std::cout << "color CC" << std::endl;

        PMP::connected_components(pmesh
          , fccmap
          , PMP::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
          .face_index_map(fim));

        BOOST_FOREACH(face_descriptor f, faces(pmesh))
          f->set_patch_id(fccmap[f]);

        to_skip.insert(selection_item->polyhedron_item());

        selection_item->changed_with_poly_item();
      }
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}


#include "Join_and_split_polyhedra_plugin.moc"
