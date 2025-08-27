#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Three.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_data_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <QAction>
#include <QList>
#include <QMessageBox>

#include <functional>

using namespace CGAL::Three;

using Vertex = T3::Vertex;
using Cell = T3::Cell;
using Wp = Geom_traits::Weighted_point_3;

template <typename CDT>
struct Vertex_converter
{
  const CDT* tr;

  Vertex_converter(const CDT* tr) : tr(tr) {}

  template <typename Src> auto operator()(const Src& src) const
  {
    using namespace CGAL;

    Vertex v;
    v.set_point(Wp{src.point()});
      switch(src.ccdt_3_data().vertex_type()) {
      case CDT_3_vertex_type::CORNER:
        v.set_dimension(0);
        v.set_index(0);
        break;
      case CDT_3_vertex_type::STEINER_ON_EDGE:
        v.set_dimension(1);
        v.set_index(static_cast<int>(src.ccdt_3_data().constrained_polyline_id(*tr).index()));
        break;
      case CDT_3_vertex_type::STEINER_IN_FACE:
        v.set_dimension(2);
        v.set_index(src.ccdt_3_data().face_index());
        break;
      case CDT_3_vertex_type::FREE:
        v.set_dimension(3);
        v.set_index(1);
        break;
      default:
        CGAL_error();
        break;
      }
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

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;

    using Vbb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
    using Vb  = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vbb>;

    using Cbb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
    using Cb  = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cbb>;

    using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
    using Tr  = CGAL::Triangulation_3<K, Tds>;
    using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;
    CDT cdt = std::invoke([&] {
      CDT cdt;
      if(mesh_item) {
        auto* const mesh = mesh_item ? mesh_item->face_graph() : nullptr;
        if(!mesh) return cdt;

        auto patch_id_pmap_opt = mesh->property_map<SMesh::Face_index, int>("f:patch_id");

        if(patch_id_pmap_opt.has_value()) {
          cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CDT>(
              *mesh, CGAL::parameters::plc_face_id(*patch_id_pmap_opt));
        } else {
          cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CDT>(*mesh);
        }
      }
      else
        if(soup_item) {
          cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CDT>(
              soup_item->points(), soup_item->polygons());
      }
      return cdt;
    });
    auto triangulation_item = std::make_unique<Scene_c3t3_item>();
    auto& item_tr = triangulation_item->triangulation();

    namespace Tet_remesh = CGAL::Tetrahedral_remeshing;
    const auto cdt_tr = Tet_remesh::get_remeshing_triangulation(std::move(cdt));
    auto inf_v =
        item_tr.tds().copy_tds(cdt_tr.tds(), cdt_tr.infinite_vertex(), Vertex_converter(&cdt), Cell_converter());
    item_tr.set_infinite_vertex(inf_v);
    triangulation_item->c3t3_changed();

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
