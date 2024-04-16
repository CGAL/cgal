#ifndef SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_H
#define SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_H

#if defined( scene_aff_transformed_polygon_soup_item_EXPORTS)
#  define SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_EXPORT Q_DECL_IMPORT
#endif

#include "Kernel_type.h"
#include "Scene_polygon_soup_item.h"
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

struct Scene_aff_transformed_polygon_soup_item_priv
{
  using Point = Kernel::Point_3;

  Scene_polygon_soup_item* ps_item;
  mutable std::vector<float> positions_lines;
  mutable std::size_t nb_lines;

  CGAL::qglviewer::Vec center_;

public:
  Scene_aff_transformed_polygon_soup_item_priv(Scene_polygon_soup_item* ps_item,
                                               const CGAL::qglviewer::Vec& pos)
    : ps_item(ps_item),
      center_(pos)
  {
    nb_lines = 0;
  }

  ~Scene_aff_transformed_polygon_soup_item_priv() { }

  void compute_elements() const
  {
    if(!ps_item)
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    const auto& points = ps_item->points();
    const auto& polygons = ps_item->polygons();

    for(std::size_t i=0; i<polygons.size(); ++i)
    {
      const std::size_t size = polygons[i].size();
      for(std::size_t j=0; j<size; ++j)
      {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][(j+1 < size) ? j+1: 0];
        const Point& a = points[i0];
        const Point& b = points[i1];

        positions_lines.push_back(a.x() - center_.x);
        positions_lines.push_back(a.y() - center_.y);
        positions_lines.push_back(a.z() - center_.z);

        positions_lines.push_back(b.x() - center_.x);
        positions_lines.push_back(b.y() - center_.y);
        positions_lines.push_back(b.z() - center_.z);
      }
    }

    QApplication::restoreOverrideCursor();
  }
};

class SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_EXPORT Scene_aff_transformed_polygon_soup_item
  : public Scene_aff_transformed_item
{
  Q_OBJECT

protected:
  friend Scene_aff_transformed_polygon_soup_item_priv;
  Scene_aff_transformed_polygon_soup_item_priv* d;

public:
  Scene_aff_transformed_polygon_soup_item(Scene_polygon_soup_item* item,
                                          const CGAL::qglviewer::Vec& pos);

  ~Scene_aff_transformed_polygon_soup_item();

  Scene_polygon_soup_item* item() { return d->ps_item; }
  const CGAL::qglviewer::Vec& center() const override { return d->center_; }

  CGAL::Three::Scene_item* clone() const override { return nullptr; }
  QString name() const override { return tr("%1_transformed").arg(d->ps_item->name()); }
  QString toolTip() const override;

  void updateCache();

  bool isEmpty() const override { return (d->nb_lines == 0); }

  virtual bool supportsRenderingMode(RenderingMode m) const override { return m == Wireframe ; }
  virtual void invalidateOpenGLBuffers() override;
  void initializeBuffers(CGAL::Three::Viewer_interface*) const override;

  void compute_bbox() const override;
  void computeElements() const override;
  void drawEdges(CGAL::Three::Viewer_interface*) const override;
};

#endif // SCENE_AFF_TRANSFORMED_POLYGON_SOUP_ITEM_H
