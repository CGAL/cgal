#ifndef SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_H
#define SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_H

#if defined( scene_aff_transformed_point_set_item_EXPORTS)
#  define SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_EXPORT Q_DECL_IMPORT
#endif

#include "Kernel_type.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_aff_transformed_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Three.h>

#include <vector>

using namespace CGAL::Three;

struct Scene_aff_transformed_point_set_item_priv
{
  using Point = Kernel::Point_3;

  Scene_points_with_normal_item* pts_item;
  mutable std::vector<float> points;
  mutable std::size_t nb_points;

  CGAL::qglviewer::Vec center_;

  Scene_aff_transformed_point_set_item_priv(Scene_points_with_normal_item *pts_item,
                                            const CGAL::qglviewer::Vec& pos)
    : pts_item(pts_item),
      center_(pos)
  {
    nb_points = 0;
  }

  ~Scene_aff_transformed_point_set_item_priv() { }

  void compute_elements() const
  {
    const Point_set& ps = *(pts_item->point_set());
    points.reserve(3*ps.size());
    for(Point_set::const_iterator it = ps.begin(); it != ps.first_selected(); it++)
    {
      const Point& p = ps.point(*it);
      points.push_back(p.x() - center_.x);
      points.push_back(p.y() - center_.y);
      points.push_back(p.z() - center_.z);
    }
  }
};

class SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_EXPORT Scene_aff_transformed_point_set_item
  : public Scene_aff_transformed_item
{
  Q_OBJECT

  using Point_set = Point_set_3<Kernel>;

protected:
  friend Scene_aff_transformed_point_set_item_priv;
  Scene_aff_transformed_point_set_item_priv* d;

public:
  Scene_aff_transformed_point_set_item(Scene_points_with_normal_item *pts_item,
                                       const CGAL::qglviewer::Vec& pos);

  ~Scene_aff_transformed_point_set_item();

  Scene_points_with_normal_item* item() { return d->pts_item; }
  const CGAL::qglviewer::Vec& center() const override { return d->center_; }

  CGAL::Three::Scene_item* clone() const override { return nullptr; }
  QString name() const override { return tr("%1_transformed").arg(d->pts_item->name()); }
  QString toolTip() const override;

  bool isEmpty() const override { return (d->nb_points == 0); }

  void updateCache();

  virtual bool supportsRenderingMode(RenderingMode m) const override { return m == Points ; }

  virtual void invalidateOpenGLBuffers() override;
  void initializeBuffers(CGAL::Three::Viewer_interface* v) const override;

  void compute_bbox() const override;
  void computeElements() const override;
  void drawPoints(CGAL::Three::Viewer_interface *viewer) const override;
};

#endif // SCENE_AFF_TRANSFORMED_POINT_SET_ITEM_H
