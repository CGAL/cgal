#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/make_constrained_Delaunay_triangulation_3.h>
#include <QAction>
#include <QList>

using namespace CGAL::Three;

using Vertex = T3::Vertex;
using Cell = T3::Cell;
using Wp = Geom_traits::Weighted_point_3;

struct Vertex_converter
{
  template <typename Src> auto operator()(const Src& src) const
  {
    Vertex v;
    v.set_point(Wp{src.point()});
    return v;
  }

  template <typename T1, typename T2> void operator()(const T1&, T2&) const {}
};

struct Cell_converter
{
  template <typename Src> auto operator()(const Src& input_c) const
  {
    Cell result_c;
    result_c.set_subdomain_index(input_c.subdomain_index());
    for(int i = 0; i < 4; ++i) {
      result_c.set_surface_patch_index(i, input_c.surface_patch_index(i));
    }
    return result_c;
  }

  template <typename T1, typename T2> void operator()(const T1&, T2&) const {}
};

class CDT_3_plugin : public QObject, public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

  QAction* actionCDT;
  CGAL::Three::Scene_interface* scene;

  void cdt_3()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    auto* mesh = item->face_graph();
    if(!mesh)
      return;
    using CDT = CGAL::Constrained_Delaunay_triangulation_3<EPICK>;
    CDT cdt = CGAL::make_constrained_Delaunay_triangulation_3(*mesh);
    const auto& cdt_tr = cdt.triangulation();
    auto triangulation_item = std::make_unique<Scene_c3t3_item>();
    auto& item_tr = triangulation_item->triangulation();

    auto inf_v = item_tr.tds().copy_tds(cdt_tr.tds(), cdt_tr.infinite_vertex(), Vertex_converter(), Cell_converter());
    item_tr.set_infinite_vertex(inf_v);
    triangulation_item->triangulation_changed();

    triangulation_item->setParent(item->parent());
    triangulation_item->setName("CDT of " + item->name());
    triangulation_item->show_intersection(false);
    scene->addItem(triangulation_item.release());
  }

public:
  void init(QMainWindow*, CGAL::Three::Scene_interface* scene, Messages_interface*) override
  {
    this->scene = scene;
    actionCDT = new QAction("3D Constrained Delaunay Triangulation", this);
    connect(actionCDT, &QAction::triggered, this, &CDT_3_plugin::cdt_3);
  }

  bool applicable(QAction* action) const override
  {
    if(action != actionCDT)
      return false;
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const override { return QList<QAction*>() << actionCDT; }
};

#include "CDT_3_plugin.moc"
