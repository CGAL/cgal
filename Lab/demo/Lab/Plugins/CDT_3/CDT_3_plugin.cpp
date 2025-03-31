#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Three.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <QAction>
#include <QList>
#include <QMessageBox>

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
    v.set_dimension(src.in_dimension());
    v.set_index(src.index());
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
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "cdt_3_plugin.json")

  QAction* actionCDT;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;

  void cdt_3()
  {
    CGAL::Three::OverrideCursorScopeGuard guard(Qt::WaitCursor);
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_surface_mesh_item* mesh_item = qobject_cast<Scene_surface_mesh_item*>(item);
    Scene_polygon_soup_item* soup_item = qobject_cast<Scene_polygon_soup_item*>(item);
    if(mesh_item)
      item = mesh_item;
    else if(soup_item)
      item = soup_item;
    else
      item = nullptr;

    if(item == nullptr) {
      CGAL::Three::Three::warning(tr("This function is only applicable on PLCs and polygon soups."));
      return;
    }

    auto* const mesh = mesh_item ? mesh_item->face_graph() : nullptr;
    if(!mesh) return;

    using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<EPICK>;
    CDT cdt = std::invoke([&] {
      CDT cdt;
      if(mesh_item) {
        auto patch_id_pmap_opt = mesh->property_map<SMesh::Face_index, int>("f:patch_id");

        if(patch_id_pmap_opt.has_value()) {
          cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(
            *mesh, CGAL::parameters::face_patch_map(*patch_id_pmap_opt));
        }
        else {
          cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(*mesh);
        }
      }
      else if(soup_item) {
        cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(
            soup_item->points(), soup_item->polygons());

      }
      return cdt;
    });
    auto triangulation_item = std::make_unique<Scene_c3t3_item>();
    auto& item_tr = triangulation_item->triangulation();

    const auto cdt_tr = CGAL::convert_to_triangulation_3(std::move(cdt));
    auto inf_v = item_tr.tds().copy_tds(cdt_tr.tds(), cdt_tr.infinite_vertex(), Vertex_converter(), Cell_converter());
    item_tr.set_infinite_vertex(inf_v);
    triangulation_item->triangulation_changed();

    triangulation_item->setParent(item->parent());
    triangulation_item->setName("CDT of " + item->name());
    triangulation_item->show_intersection(false);
    scene->addItem(triangulation_item.release());
    item->setVisible(false);
  }

public:
  void init(QMainWindow* mw, CGAL::Three::Scene_interface* scene, Messages_interface*) override
  {
    this->scene = scene;
    this->mw = mw;
    actionCDT = new QAction("3D Constrained Delaunay Triangulation", this);
    connect(actionCDT, &QAction::triggered, this, &CDT_3_plugin::cdt_3);
  }

  bool applicable(QAction* action) const override
  {
    if(action != actionCDT)
      return false;
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_polygon_soup_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const override { return QList<QAction*>() << actionCDT; }
};

#include "CDT_3_plugin.moc"
