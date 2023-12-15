#ifndef SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_H
#define SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_H

#if defined( scene_aff_transformed_surface_mesh_item_EXPORTS)
#  define SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_EXPORT Q_DECL_IMPORT
#endif

#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_aff_transformed_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <vector>

using namespace CGAL::Three;

struct Scene_aff_transformed_surface_mesh_item_priv
{
  using Point = Kernel::Point_3;

  Scene_surface_mesh_item* sm_item;
  mutable std::vector<float> positions_lines;
  mutable std::size_t nb_lines;

  CGAL::qglviewer::Vec center_;

public:
  Scene_aff_transformed_surface_mesh_item_priv(Scene_surface_mesh_item* sm_item,
                                               const CGAL::qglviewer::Vec& pos)
    : sm_item(sm_item),
      center_(pos)
  {
    nb_lines = 0;
  }

  ~Scene_aff_transformed_surface_mesh_item_priv() { }

  void compute_elements() const
  {
    if(!sm_item)
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    const auto& sm_ptr = sm_item->face_graph();
    auto vpm = get(CGAL::vertex_point, *sm_ptr);

    positions_lines.resize(0);
    for(auto e : edges(*sm_ptr))
    {
      const Point a = get(vpm, target(halfedge(e, *sm_ptr), *sm_ptr));
      const Point b = get(vpm, target(opposite(halfedge(e, *sm_ptr), *sm_ptr), *sm_ptr));

      positions_lines.push_back(a.x() - center_.x);
      positions_lines.push_back(a.y() - center_.y);
      positions_lines.push_back(a.z() - center_.z);

      positions_lines.push_back(b.x() - center_.x);
      positions_lines.push_back(b.y() - center_.y);
      positions_lines.push_back(b.z() - center_.z);
    }

    QApplication::restoreOverrideCursor();
  }
};

class SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_EXPORT Scene_aff_transformed_surface_mesh_item
  : public Scene_aff_transformed_item
{
  Q_OBJECT

protected:
  friend Scene_aff_transformed_surface_mesh_item_priv;
  Scene_aff_transformed_surface_mesh_item_priv* d;

public:
  Scene_aff_transformed_surface_mesh_item(Scene_surface_mesh_item* item,
                                          const CGAL::qglviewer::Vec& pos);

  ~Scene_aff_transformed_surface_mesh_item();

  Scene_surface_mesh_item* item() { return d->sm_item; }
  const CGAL::qglviewer::Vec& center() const override { return d->center_; }

  CGAL::Three::Scene_item* clone() const override { return nullptr; }
  QString name() const override { return tr("%1_transformed").arg(d->sm_item->name()); }
  QString toolTip() const override;

  bool isEmpty() const override { return (d->nb_lines == 0); }

  void updateCache();

  virtual bool supportsRenderingMode(RenderingMode m) const override { return m == Wireframe ; }
  virtual void invalidateOpenGLBuffers() override;
  void initializeBuffers(CGAL::Three::Viewer_interface*) const override;

  void compute_bbox() const override;
  void computeElements() const override;
  void drawEdges(CGAL::Three::Viewer_interface*) const override;
};

#endif // SCENE_AFF_TRANSFORMED_SURFACE_MESH_ITEM_H
