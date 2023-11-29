#include "config.h"

#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "ui_Offset_meshing_dialog.h"

#include "C3t3_type.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Timer.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <CGAL/Three/Three.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMenu>
#include <QApplication>
#include <QtPlugin>
#include <QThread>
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"
#include <QInputDialog>
#include <QStringList>
#include <QMessageBox>
#include <QAbstractButton>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>

using namespace CGAL::Three;

namespace CGAL {

template <class TriangleMesh, class GeomTraits>
class Offset_function
{
  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh>;
  using Traits = AABB_traits<GeomTraits, Primitive>;
  using Tree = AABB_tree<Traits>;
  using Side_of = Side_of_triangle_mesh<TriangleMesh, GeomTraits>;

  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

public:
  Offset_function(const TriangleMesh& tm,
                  double offset_distance)
    : m_tree_ptr(std::make_shared<Tree>(std::begin(faces(tm)), std::end(faces(tm)), tm)),
      m_side_of_ptr(std::make_shared<Side_of>(*m_tree_ptr)),
      m_is_inset(offset_distance < 0),
      m_sq_offset_distance(CGAL::square(offset_distance)),
      m_is_closed(is_closed(tm))
  {
    CGAL_assertion(!m_tree_ptr->empty());
  }

  // we only need negative inside, and positive outside, so we can compare square roots
  double operator()(const Point_3& p) const
  {
    const Bounded_side side = m_is_closed ? m_side_of_ptr->operator()(p) : ON_UNBOUNDED_SIDE;

    if(m_is_inset) // also means that the mesh is closed
    {
      // - ON_UNBOUNDED_SIDE is outside the offset since we are insetting
      // - ON_BOUNDARY is outside the offset since we are insetting
      if(side != ON_BOUNDED_SIDE)
        return 1;

      // inside the offset if the distance to the input mesh is greater than the offset distance
      const FT sq_distance = m_tree_ptr->squared_distance(p);
      return (sq_distance > m_sq_offset_distance) ? -1 : 1;
    }
    else // outset
    {
      // - ON_BOUNDED_SIDE can only happen if it's a closed mesh, and in that case, being inside
      // the mesh is being inside the offset
      // - ON_BOUNDARY is in the offset whether the mesh is open or closed
      if(side != ON_UNBOUNDED_SIDE)
        return - 1;

      // inside the offset if the distance to the input mesh is smaller than the offset distance
      const FT sq_distance = m_tree_ptr->squared_distance(p);
      return (sq_distance < m_sq_offset_distance) ? -1 : 1;
    }
  }

private:
  std::shared_ptr<Tree> m_tree_ptr;
  std::shared_ptr<Side_of> m_side_of_ptr;
  const bool m_is_inset;
  const double m_sq_offset_distance;
  const bool m_is_closed;
};

template <typename Points, typename Polygons>
class Polygon_soup_offset_function
{
  using Polygon_iterator = typename Polygons::const_iterator;

  class Polygon_soup_point_property_map
  {
    const Points* points_vector_ptr;

  public:
    using key_type = Polygon_iterator;
    using value_type = EPICK::Point_3;
    using reference = const value_type&;
    using category = boost::readable_property_map_tag;

    Polygon_soup_point_property_map() = default;
    Polygon_soup_point_property_map(const Points* ptr) : points_vector_ptr(ptr) { }

    friend reference get(Polygon_soup_point_property_map map,
                         key_type polygon_it)
    {
      return (*map.points_vector_ptr)[*polygon_it->begin()];
    }
  };

  class Polygon_soup_triangle_property_map
  {
    const Points* points_vector_ptr;

  public:
    using key_type = Polygon_iterator;
    using value_type = EPICK::Triangle_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    Polygon_soup_triangle_property_map() = default;
    Polygon_soup_triangle_property_map(const Points* ptr) : points_vector_ptr(ptr) { }

    friend value_type get(Polygon_soup_triangle_property_map map,
                          key_type polygon_it)
    {
      auto it = polygon_it->begin();
      CGAL_assertion(it != polygon_it->end());
      const auto id0 = *it++;
      CGAL_assertion(it != polygon_it->end());
      const auto id1 = *it++;
      CGAL_assertion(it != polygon_it->end());
      const auto id2 = *it++;
      CGAL_assertion(it == polygon_it->end());

      return value_type((*map.points_vector_ptr)[id0],
                        (*map.points_vector_ptr)[id1],
                        (*map.points_vector_ptr)[id2]);
    }
  };

  struct AABB_polygon_soup_triangle_primitive
    : public CGAL::AABB_primitive<Polygon_iterator,
                                  Polygon_soup_triangle_property_map,
                                  Polygon_soup_point_property_map,
                                  CGAL::Tag_true /*ExternalPropertyMaps*/,
                                  CGAL::Tag_false /*CacheDatum*/>
  {
    using Base = CGAL::AABB_primitive<Polygon_iterator,
                                      Polygon_soup_triangle_property_map,
                                      Polygon_soup_point_property_map,
                                      CGAL::Tag_true,
                                      CGAL::Tag_false>;

    using Id = Polygon_iterator;

    template <typename ObjectPmap, typename PointPmap>
    AABB_polygon_soup_triangle_primitive(Id id,
                                         ObjectPmap&& opmap,
                                         PointPmap&& ppmap)
      : Base(id, std::forward<ObjectPmap>(opmap), std::forward<PointPmap>(ppmap))
    {
    }

    template <typename Iterator, typename ObjectPmap, typename PointPmap>
    AABB_polygon_soup_triangle_primitive(Iterator it,
                                         ObjectPmap&& opmap,
                                         PointPmap&& ppmap)
      : Base(*it, std::forward<ObjectPmap>(opmap), std::forward<PointPmap>(ppmap))
    {
    }
  }; // struct template Polygon_soup_primitive

  using AABB_traits = CGAL::AABB_traits<EPICK, AABB_polygon_soup_triangle_primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;

  std::shared_ptr<AABB_tree> m_tree_ptr;
  double m_sq_offset_distance;

public:
  Polygon_soup_offset_function(const Points& points,
                               const Polygons& polygons,
                               const double offset_distance)
    : m_tree_ptr(std::make_shared<AABB_tree>(std::begin(polygons),
                                             std::end(polygons),
                                             Polygon_soup_triangle_property_map(&points),
                                             Polygon_soup_point_property_map(&points))),
      m_sq_offset_distance(square(offset_distance))
  {
    CGAL_assertion(!m_tree_ptr->empty());
  }

  // we only need negative inside, and positive outside, so we can compare square roots
  double operator()(const EPICK::Point_3& p) const
  {
    // it's a soup so it's open by definition ==> treat inset and outset identically
    const double sq_distance = m_tree_ptr->squared_distance(p);
    return sq_distance - m_sq_offset_distance;
  }

}; // class Polygon_soup_offset_function

} // namespace CGAL

CGAL::Offset_function<SMesh, EPICK>
offset_function(Scene_surface_mesh_item* item, double offset_value)
{
  return { *(item->face_graph()), offset_value };
}

CGAL::Polygon_soup_offset_function<Scene_polygon_soup_item::Points,
                                   Scene_polygon_soup_item::Polygons>
offset_function(Scene_polygon_soup_item* item, double offset_value)
{
  return { item->points(), item->polygons(), offset_value };
}

class MeshGuard
{
  SMesh* mesh;
  bool done;

public:
  MeshGuard(SMesh* mesh) : mesh(mesh), done(false) { }
  void setDone() { done = true; }
  ~MeshGuard()
  {
    if(!done)
      delete mesh;
  }
};

// declare the CGAL function
template<class SourceItem>
SMesh* cgal_off_meshing(QWidget*,
                        SourceItem* source_item,
                        Scene_polylines_item* polylines_item,
                        const double offset_value,
                        const double angle,
                        const double sizing,
                        const double approx,
                        const double edge_size,
                        int tag)
{
  using GT = EPICK;
  using Sphere_3 = GT::Sphere_3;

  using Mesh_domain_base = CGAL::Labeled_mesh_domain_3<GT, int, int>;
  using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain_base>;
  using Tr = C3t3::Triangulation;
  using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

  const CGAL::Bbox_3 bbox = source_item->bbox();

  const GT::Point_3 center((bbox.xmax() + bbox.xmin()) / 2,
                           (bbox.ymax() + bbox.ymin()) / 2,
                           (bbox.zmax() + bbox.zmin()) / 2);
  const double rad = 0.6 * std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                     CGAL::square(bbox.ymax() - bbox.ymin()) +
                                     CGAL::square(bbox.zmax() - bbox.zmin()))
                     + offset_value;
  const double sqrad = CGAL::square(rad);

  CGAL::Timer timer;
  timer.start();

  namespace p = CGAL::parameters;

  Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain
    (p::function = offset_function(source_item, offset_value),
     p::bounding_object = Sphere_3(center, sqrad),
     p::relative_error_bound = 1e-7,
     p::construct_surface_patch_index = [](int i, int j) { return (i * 1000 + j); });

  const CGAL::Mesh_facet_topology topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH;
  auto manifold_option = p::non_manifold();
  if(tag == 1)
    manifold_option = p::manifold_with_boundary();
  if(tag == 2)
    manifold_option = p::manifold();

  Mesh_criteria criteria(p::facet_angle = angle,
                         p::facet_size = sizing,
                         p::facet_distance = approx,
                         p::facet_topology = topology,
                         p::edge_size = edge_size);

  if(polylines_item != nullptr)
  {
    typedef std::vector<Mesh_domain::Surface_patch_index> Surface_patch_ids;
    std::vector<Mesh_domain::Surface_patch_index> surface_patch_ids;

    domain.add_features_and_incidences(polylines_item->polylines.begin(),
                                       polylines_item->polylines.end(),
                                       CGAL::Identity_property_map<Scene_polylines_item::Polyline>(),
                                       CGAL::Constant_property_map<Scene_polylines_item::Polyline, Surface_patch_ids>(surface_patch_ids));
  }

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      p::no_perturb(),
                                      p::no_exude(),
                                      manifold_option);

  timer.stop();
  std::cerr << "done (" << timer.time() << " ms, " << c3t3.triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(c3t3.number_of_facets_in_complex() > 0)
  {
    SMesh* pRemesh = new SMesh();

    // if the thread is interrupted before the mesh is returned, delete it.
    MeshGuard guard(pRemesh);
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, *pRemesh);
    guard.setDone();

    CGAL_postcondition(CGAL::Polygon_mesh_processing::is_outward_oriented(*pRemesh));

    return pRemesh;
  }
  else
  {
    return nullptr;
  }
}

struct Mesher_thread
  : public QThread
{
  Q_OBJECT

private:
  Scene_surface_mesh_item* sm_item;
  Scene_polygon_soup_item* soup_item;
  Scene_polylines_item* polylines_item;

  const double offset_value;
  const double angle;
  const double sizing;
  const double approx;
  const double edge_size;
  int tag_index;

public:
  Mesher_thread(Scene_surface_mesh_item* sm_item,
                Scene_polygon_soup_item* soup_item,
                Scene_polylines_item* polylines_item,
                const double offset_value,
                const double angle,
                const double sizing,
                const double approx,
                const double edge_size,
                int tag)
    : sm_item(sm_item), soup_item(soup_item), polylines_item(polylines_item),
      offset_value(offset_value),
      angle(angle), sizing(sizing), approx(approx), edge_size(edge_size), tag_index(tag)
  {
  }

  void run() override
  {
    SMesh* offset_mesh = nullptr;

    if(soup_item)
    {
      offset_mesh = cgal_off_meshing(Three::mainWindow(),
                                     soup_item, polylines_item,
                                     offset_value,
                                     angle, sizing, approx, edge_size, tag_index);
    }
    else
    {
      offset_mesh = cgal_off_meshing(Three::mainWindow(),
                                     sm_item, polylines_item,
                                     offset_value,
                                     angle, sizing, approx, edge_size, tag_index);
    }

    Three::getMutex()->lock();
    Three::getWaitCondition()->wakeAll();
    Three::getMutex()->unlock();

    Q_EMIT resultReady(offset_mesh);
  }

Q_SIGNALS:
  void resultReady(SMesh *offset_mesh);
};

class Polyhedron_demo_offset_meshing_plugin
  : public QObject,
    protected Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionOffsetMeshing;
  QAction* actionInflateMesh;

  Scene_interface *scene;
  QMainWindow *mw;

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionOffsetMeshing = new QAction(tr("Offset Meshing"), mw);
    actionOffsetMeshing->setProperty("subMenuName", "3D Surface Mesh Generation");
    connect(actionOffsetMeshing, SIGNAL(triggered()),
            this, SLOT(offset_meshing()));

    actionInflateMesh = new QAction(tr("Inflate Mesh"), mw);
    actionInflateMesh->setProperty("subMenuName", "Operations on Polyhedra");
    connect(actionInflateMesh, SIGNAL(triggered()),
            this, SLOT(inflate_mesh()));
  }

  bool applicable(QAction* action) const
  {
    if(action == actionOffsetMeshing)
    {
      if(scene->selectionIndices().size() == 1)
      {
        const int index = scene->mainSelectionIndex();
        return (qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) ||
                qobject_cast<Scene_polygon_soup_item*>(scene->item(index)));
      }

      // Can provide a polyline item for feature protection
      if(scene->selectionIndices().size() != 2)
        return false;

      // One needs to be a surface mesh or polygon soup item, and the other a polyline item
      const int index1 = scene->selectionIndices().at(0);
      const int index2 = scene->selectionIndices().at(1);
      Scene_item* item1 = scene->item(index1);
      Scene_item* item2 = scene->item(index2);

      if((qobject_cast<Scene_surface_mesh_item*>(item1) ||
          qobject_cast<Scene_polygon_soup_item*>(item1)) &&
         qobject_cast<Scene_polylines_item*>(item2))
        return true;

      if((qobject_cast<Scene_surface_mesh_item*>(item2) ||
          qobject_cast<Scene_polygon_soup_item*>(item2)) &&
         qobject_cast<Scene_polylines_item*>(item1))
        return true;
    }
    else if(action == actionInflateMesh)
    {
      if(scene->selectionIndices().size() == 1)
      {
        const int index = scene->mainSelectionIndex();
        return qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      }
    }

    return false;
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionOffsetMeshing
                             << actionInflateMesh;
  }

public Q_SLOTS:
  void offset_meshing();
  void inflate_mesh();
}; // class Polyhedron_demo_offset_meshing_plugin

void
Polyhedron_demo_offset_meshing_plugin::
offset_meshing()
{
  Scene_item* item = nullptr;
  Scene_surface_mesh_item* sm_item = nullptr;
  Scene_polygon_soup_item* soup_item = nullptr;
  Scene_polylines_item* polylines_item = nullptr;

  bool mesh_or_soup_item_found = false;
  for(Scene_interface::Item_id index : scene->selectionIndices())
  {
    if(!mesh_or_soup_item_found)
    {
      sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      if(sm_item == nullptr)
      {
        soup_item = qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
        if(soup_item != nullptr)
        {
          item = scene->item(index);
          mesh_or_soup_item_found = true;
          continue;
        }
      }
      else
      {
        item = scene->item(index);
        mesh_or_soup_item_found = true;
        continue;
      }
    }

    polylines_item = qobject_cast<Scene_polylines_item*>(scene->item(index));
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if(!mesh_or_soup_item_found)
    return;

  if(sm_item)
  {
    if(!is_triangle_mesh(*(sm_item->face_graph())))
    {
      QMessageBox::critical(mw,
                            tr("Offset Meshing"),
                            tr("The selected mesh is not a triangle mesh."));
      return;
    }
  }
  else
  {
    for(const auto& p : soup_item->polygons())
    {
      if(p.size() != 3)
      {
        QMessageBox::critical(mw,
                              tr("Offset Meshing"),
                              tr("The selected polygon soup is not a triangle soup."));
        return;
      }
    }
  }

  double diag;
  if(sm_item)
    diag = sm_item->bboxDiagonal();
  else
    diag = soup_item->bboxDiagonal();

  QApplication::restoreOverrideCursor();

  bool ok = true;
  double offset_value = QInputDialog::getDouble(mw,
                                                QString("Choose Offset Value"),
                                                QString("Offset Value (use a negative number to compute the inset of a closed mesh)"),
                                                0.1 * diag,
                                                - (std::numeric_limits<double>::max)(),
                                                (std::numeric_limits<double>::max)(), 10, &ok);
  if(!ok)
    return;

  if(offset_value < 0 && (!sm_item || !is_closed(*(sm_item->face_graph()))))
  {
    QMessageBox::critical(mw,
                          tr("Offset Meshing"),
                          tr("Insetting is only possible for closed polygon meshes."));
    return;
  }

  QDialog dialog(mw);
  Ui::Offset_meshing_dialog ui;
  ui.setupUi(&dialog);
  ui.angle->setRange(1.0, 30.0);

  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  ui.sizing->setRange(diag * 10e-6, diag);
  ui.sizing->setValue(diag * 0.05); // default value
  ui.approx->setRange(diag * 10e-7, diag);
  ui.approx->setValue(diag * 0.005);

  if(polylines_item != nullptr)
  {
    ui.edge_sizing->setRange(diag * 10e-6, diag);
    ui.edge_sizing->setValue(diag * 0.05); // default value
  }
  else
  {
    ui.edge_sizing->setEnabled(false);
  }

  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;

  const double angle = ui.angle->value();
  const double approx = ui.approx->value();
  const double sizing = ui.sizing->value();
  const double edge_size = (polylines_item != nullptr) ? ui.edge_sizing->value() : 0;
  const int tag_index = ui.tags->currentIndex();

  if(tag_index < 0)
    return;

  QApplication::setOverrideCursor(Qt::BusyCursor);

  std::cerr << "mesh with:"
            << "\n  angle= " << angle
            << "\n  sizing= " << sizing
            << "\n  approx= " << approx
            << "\n  tag= " << tag_index
            << std::boolalpha
            << std::endl;

  Mesher_thread* worker = nullptr;
  if(soup_item)
  {
    worker = new Mesher_thread(nullptr, soup_item, polylines_item,
                               offset_value,
                               angle, sizing, approx, edge_size, tag_index);
  }
  else
  {
    worker = new Mesher_thread(sm_item, nullptr, polylines_item,
                               offset_value,
                               angle, sizing, approx, edge_size, tag_index);
  }

  connect(worker, &QThread::finished,
          worker, &QObject::deleteLater);

  connect(worker, &Mesher_thread::resultReady,
          this, [item, angle, sizing, approx, offset_value/* , index */](SMesh* offset_mesh)
                {
                  if(!offset_mesh)
                  {
                    QApplication::restoreOverrideCursor();

                    Three::getMutex()->lock();
                    Three::isLocked() = false;
                    Three::getMutex()->unlock();

                    return;
                  }

                  Scene_surface_mesh_item* offset_item = new Scene_surface_mesh_item(offset_mesh);
                  offset_item->setName(tr("%1 offset %5 (%2 %3 %4)").arg(item->name())
                                                                   .arg(angle)
                                                                   .arg(sizing)
                                                                   .arg(approx)
                                                                   .arg(offset_value));
                  offset_item->setColor(Qt::magenta);
                  offset_item->setRenderingMode(Wireframe);
                  Three::scene()->addItem(offset_item);

                  QApplication::restoreOverrideCursor();

                  Three::getMutex()->lock();
                  Three::isLocked() = false;
                  Three::getMutex()->unlock();
                });

  QMessageBox* message_box = new QMessageBox(QMessageBox::NoIcon,
                                             "Meshing",
                                             "Offset meshing in progress...",
                                             QMessageBox::Cancel,
                                             mw);
  message_box->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));

  connect(cancelButton, &QAbstractButton::clicked,
          this, [worker](){ worker->terminate(); });
  connect(worker, &Mesher_thread::finished,
          message_box, &QMessageBox::close);

  Three::getMutex()->lock();
  Three::isLocked() = true;
  Three::getMutex()->unlock();

  message_box->open();
  worker->start();
}

void
Polyhedron_demo_offset_meshing_plugin::
inflate_mesh()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_item* item = scene->item(index);
  if(item == nullptr)
    return;

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
  if(sm_item == nullptr)
    return;

  SMesh* sMesh = sm_item->face_graph();
  if(sMesh == nullptr)
    return;

  const double diag = sm_item->bboxDiagonal();

  bool ok = true;
  const double offset_value = QInputDialog::getDouble(mw,
                                                      QString("Choose Inflate Distance"),
                                                      QString("Inflate Distance (use a negative number to deflate)"),
                                                      0.1 * diag,
                                                      -(std::numeric_limits<double>::max)(),
                                                      (std::numeric_limits<double>::max)(),
                                                      10,
                                                      &ok);
  if(!ok)
    return;

  auto vpm = get(CGAL::vertex_point, *sMesh);
  auto vnm = sMesh->property_map<vertex_descriptor, EPICK::Vector_3 >("v:normal").first;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  for(const auto& v : vertices(*sMesh))
  {
    const EPICK::Vector_3& n = get(vnm, v);
    put(vpm, v, get(vpm, v) + offset_value * n);
  }

  sm_item->invalidateOpenGLBuffers();
  sm_item->itemChanged();
  sm_item->itemVisibilityChanged();

  QApplication::restoreOverrideCursor();
}

#include "Offset_meshing_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
