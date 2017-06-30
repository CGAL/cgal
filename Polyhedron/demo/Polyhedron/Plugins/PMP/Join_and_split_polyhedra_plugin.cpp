#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#ifdef USE_SURFACE_MESH
#include "Scene_surface_mesh_item.h"
#include <CGAL/Mesh_3/properties_Surface_mesh.h>
#else
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Mesh_3/properties_Polyhedron_3.h>
#endif
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_map.hpp>
#include "Color_map.h"

#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_facegraph_item;
#else
typedef Scene_polyhedron_item Scene_facegraph_item;
#endif
typedef Scene_facegraph_item::Face_graph FaceGraph;
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
      if (qobject_cast<Scene_facegraph_item*>(scene->item(index)))
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
    = scene->selectionIndices().first();
  Scene_facegraph_item* mainSelectionItem
    = qobject_cast<Scene_facegraph_item*>(scene->item(mainSelectionIndex));
  if(!mainSelectionItem)
  {
    std::cerr<<"No selected polyhedron_item"<<std::endl;
    return;
  }
  QList<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices()) {
    if (index == mainSelectionIndex)
      continue;

    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(item)
    {
      indices_to_remove.push_front(index);
      CGAL::copy_face_graph(*item->polyhedron(), *mainSelectionItem->polyhedron());
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
  Polyhedron_appender(std::list<FaceGraph*>& new_polyhedra):
    m_new_polyhedra(new_polyhedra) {}
  void operator()(const FaceGraph& p){
    m_new_polyhedra.push_back( new FaceGraph(p) );
  }
  std::list<FaceGraph*>& m_new_polyhedra;
};

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionSplitPolyhedra_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()) {
    colors_.clear();
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(item)
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      std::list<FaceGraph*> new_polyhedra;
      typedef boost::property_map<FaceGraph,CGAL::face_patch_id_t<int> >::type PatchIDMap;
      PatchIDMap pidmap = get(CGAL::face_patch_id_t<int>(), *item->face_graph());
      int nb_patches = CGAL::Polygon_mesh_processing::connected_components(*item->face_graph(),
                                                          pidmap);


      for(int i=0; i<nb_patches; ++i)
      {
        //std::vector<int> pids;
        //pids.push_back(i);
        CGAL::Face_filtered_graph<FaceGraph> filter_graph(*item->face_graph(), i, pidmap);
        FaceGraph* new_graph = new FaceGraph();
        CGAL::copy_face_graph(filter_graph, *new_graph);
        new_polyhedra.push_back(new_graph);
      }


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
      BOOST_FOREACH(FaceGraph* polyhedron_ptr, new_polyhedra)
      {
        Scene_facegraph_item* new_item=new Scene_facegraph_item(polyhedron_ptr);
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

  std::set<Scene_facegraph_item*> to_skip;

  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(item && to_skip.find(item) == to_skip.end())
    {
      item->setItemIsMulticolor(true);
      typedef boost::property_map<FaceGraph,CGAL::face_patch_id_t<int> >::type PatchIDMap;
      PatchIDMap pidmap = get(CGAL::face_patch_id_t<int>(), *item->face_graph());
      CGAL::Polygon_mesh_processing::connected_components(*item->face_graph(),
                                                          pidmap);

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
        typedef boost::graph_traits<FaceGraph>::face_descriptor   face_descriptor;

        selection_item->polyhedron_item()->setItemIsMulticolor(true);
#ifndef USE_SURFACE_MESH
        selection_item->polyhedron_item()->set_color_vector_read_only(false);
#endif
        FaceGraph& pmesh = *(selection_item->polyhedron());

        boost::property_map<FaceGraph, boost::face_index_t>::type fim
          = get(boost::face_index, pmesh);
        boost::vector_property_map<int,
          boost::property_map<FaceGraph, boost::face_index_t>::type>
          fccmap(fim);
        boost::property_map<FaceGraph, CGAL::face_patch_id_t<int> >::type pid
          = get(CGAL::face_patch_id_t<int>(), pmesh);

        std::cout << "color CC" << std::endl;

        PMP::connected_components(pmesh
          , fccmap
          , PMP::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
          .face_index_map(fim));

        BOOST_FOREACH(face_descriptor f, faces(pmesh))
        {
          put(pid, f, fccmap[f]);
        }

        to_skip.insert(selection_item->polyhedron_item());

        selection_item->changed_with_poly_item();
      }
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}


#include "Join_and_split_polyhedra_plugin.moc"
