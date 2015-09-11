#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>

template <typename G>
struct Is_border {
  const G& g;
  Is_border(const G& g)
    : g(g)
  {}

 template <typename Descriptor>
  bool operator()(const Descriptor& d) const {
   return is_border(d,g);
  }

  bool operator()(typename boost::graph_traits<G>::vertex_descriptor d) const {
    return is_border(d,g) != boost::none;
  }

};



class Polyhedron_demo_polyhedron_stitching_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionDetectBorders;
  QAction* actionStitchBorders;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectBorders << actionStitchBorders; }
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* /* m */)
  {
    actionDetectBorders= new QAction(tr("Detect polyhedron boundaries"), mainWindow);
    actionStitchBorders= new QAction(tr("Stitch polyhedron duplicated boundaries"), mainWindow);
    actionDetectBorders->setObjectName("actionDetectBorders");
    actionStitchBorders->setObjectName("actionStitchBorders");
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
  void on_actionDetectBorders_triggered();
  void on_actionStitchBorders_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin



struct Polyline_visitor
{
  Scene_polylines_item* new_item;

  Polyline_visitor(Scene_polylines_item* new_item)
    : new_item(new_item)
  {}

  void start_new_polyline()
  {
    new_item->polylines.push_back( Scene_polylines_item::Polyline() );
  }

  void add_node(boost::graph_traits<Polyhedron>::vertex_descriptor vd)
  {
    new_item->polylines.back().push_back(vd->point());
  }
  void end_polyline(){}
};

void Polyhedron_demo_polyhedron_stitching_plugin::on_actionDetectBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Scene_polylines_item* new_item = new Scene_polylines_item();

      Polyhedron* pMesh = item->polyhedron();
      pMesh->normalize_border();

#if 0
      for (Polyhedron::Halfedge_iterator
              it=pMesh->border_halfedges_begin(), it_end=pMesh->halfedges_end();
              it!=it_end; ++it)
      {
        if (!it->is_border()) continue;
        /// \todo build cycles and graph with nodes of valence 2.
        new_item->polylines.push_back( Scene_polylines_item::Polyline() );
        new_item->polylines.back().push_back( it->opposite()->vertex()->point() );
        new_item->polylines.back().push_back( it->vertex()->point() );
      }
#else
      typedef boost::filtered_graph<Polyhedron,Is_border<Polyhedron>, Is_border<Polyhedron> > BorderGraph;
      
      Is_border<Polyhedron> ib(*pMesh);
      BorderGraph bg(*pMesh,ib,ib);
      Polyline_visitor polyline_visitor(new_item); 
      CGAL::split_graph_into_polylines( bg,
                                        polyline_visitor,
                                        CGAL::IsTerminalDefault() );
#endif
      
      if (new_item->polylines.empty())
      {
        delete new_item;
      }
      else
      {
        new_item->setName(tr("Boundary of %1").arg(item->name()));
        new_item->setColor(Qt::red);
        scene->addItem(new_item);
        new_item->invalidate_buffers();
      }
    }
  }
}

void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if(item)
    {
      Polyhedron* pMesh = item->polyhedron();
      CGAL::Polygon_mesh_processing::stitch_borders(*pMesh);
      item->invalidate_buffers();
      scene->itemChanged(item);
    }
  }
}

#include "Polyhedron_demo_polyhedron_stitching_plugin.moc"
