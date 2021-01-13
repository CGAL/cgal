#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/function_output_iterator.hpp>
#include <boost/unordered_map.hpp>
#include "Color_map.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;
using namespace CGAL::Three;
class Polyhedron_demo_join_and_split_polyhedra_plugin:
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "join_and_split_polyhedra_plugin.json")
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
    actionSplitPolyhedra->setProperty("subMenuName", "Polygon Mesh Processing");
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
  std::vector<int> indices_to_remove;
  Q_FOREACH(int index, scene->selectionIndices()) {
    if (index == mainSelectionIndex)
      continue;

    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(item)
    {
      indices_to_remove.push_back(index);
      CGAL::copy_face_graph(*item->polyhedron(), *mainSelectionItem->polyhedron());
    }
    else
      std::cerr<<"No selected polyhedron_item"<<std::endl;
  }
  mainSelectionItem->invalidateOpenGLBuffers();
  scene->itemChanged(mainSelectionIndex);

  //remove the other items
  std::sort(indices_to_remove.begin(), indices_to_remove.end(), std::greater<int>());
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

struct Compare{
bool operator()(const FaceGraph& mesh1, const FaceGraph& mesh2)
{
  return num_faces(mesh1) < num_faces(mesh2);
}
};

void Polyhedron_demo_join_and_split_polyhedra_plugin::on_actionSplitPolyhedra_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()) {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(item)
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      std::vector<FaceGraph> new_polyhedra;
      CGAL::Polygon_mesh_processing::split_connected_components(*item->face_graph(),
                                                                new_polyhedra);
      //sort polyhedra by number of faces
      std::sort(new_polyhedra.begin(), new_polyhedra.end(), Compare());


      if (new_polyhedra.size()==1)
      {
        CGAL::Three::Three::information( tr("%1 has only one connected component").arg(item->name()) );
        QApplication::restoreOverrideCursor();
        continue;
      }

      int cc=0;
      std::vector<QColor> color_map;
      if(item->hasPatchIds())
        color_map = item->color_vector();
      else
        compute_color_map(item->color(), new_polyhedra.size(), std::back_inserter(color_map));
      Scene_group_item *group = new Scene_group_item("CC");
       scene->addItem(group);
      for(FaceGraph& poly : new_polyhedra)
      {
        Scene_facegraph_item* new_item=new Scene_facegraph_item(poly);
        new_item->setName(tr("%1 - CC %2").arg(item->name()).arg(cc));
        new_item->setColor(color_map[cc]);
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
      int nb_patch_ids = CGAL::Polygon_mesh_processing::connected_components(*item->face_graph(),
                                                                             pidmap);
      item->computeItemColorVectorAutomatically(true);
      item->invalidateOpenGLBuffers();
      item->setProperty("NbPatchIds", nb_patch_ids);
      scene->itemChanged(item);
    }
    else
    {
      Scene_polyhedron_selection_item* selection_item =
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

      if (selection_item)
      {
        if(selection_item->selected_edges.empty())
        {
          QApplication::restoreOverrideCursor();
          QMessageBox::warning(mw, "Empty Edges", "There are no selected edges. Skipping.");
          continue;
          QApplication::setOverrideCursor(Qt::WaitCursor);
        }
        namespace PMP = CGAL::Polygon_mesh_processing;
        typedef boost::graph_traits<FaceGraph>::face_descriptor   face_descriptor;

        selection_item->polyhedron_item()->setItemIsMulticolor(true);
        selection_item->polyhedron_item()->computeItemColorVectorAutomatically(true);
        FaceGraph& pmesh = *(selection_item->polyhedron());

        boost::property_map<FaceGraph, boost::face_index_t>::type fim
          = get(boost::face_index, pmesh);
        boost::vector_property_map<int,
          boost::property_map<FaceGraph, boost::face_index_t>::type>
          fccmap(static_cast<unsigned>(num_faces(pmesh)),fim);
        boost::property_map<FaceGraph, CGAL::face_patch_id_t<int> >::type pid
          = get(CGAL::face_patch_id_t<int>(), pmesh);

        std::cout << "color CC" << std::endl;

        int nb_patch_ids = PMP::connected_components(pmesh
                                                     , fccmap
                                                     , PMP::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
                                                     .face_index_map(fim));

        for(face_descriptor f : faces(pmesh))
        {
          put(pid, f, fccmap[f]);
        }

        to_skip.insert(selection_item->polyhedron_item());

        selection_item->changed_with_poly_item();
        selection_item->polyhedron_item()->setProperty("NbPatchIds", nb_patch_ids);
      }
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}


#include "Join_and_split_polyhedra_plugin.moc"
