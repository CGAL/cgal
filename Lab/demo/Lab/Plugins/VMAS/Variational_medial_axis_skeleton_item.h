#ifndef VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_H
#define VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_H

#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Simple_cartesian.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_spheres_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polygon_soup_item.h"
#include "create_sphere.h"
#include "Variational_medial_axis_skeleton_item_config.h"

struct Variational_medial_axis_skeleton_item_priv;

class VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_EXPORT Variational_medial_axis_skeleton_item:
    public CGAL::Three::Scene_group_item
{
  Scene_polygon_soup_item* faces=nullptr;
  Scene_spheres_item* spheres=nullptr;
  Scene_polylines_item* edges=nullptr;

  Variational_medial_axis_skeleton_item_priv* d;

    Q_OBJECT
  public:
    Variational_medial_axis_skeleton_item(CGAL::Three::Scene_interface* scene_interface,
                                          const Scene_surface_mesh_item* sm_item,
                                          std::size_t nb_spheres);
    ~Variational_medial_axis_skeleton_item();
    void fill_subitems();

    // bool eventFilter(QObject *, QEvent *);
};
#endif // VARIATIONAL_MEDIAL_AXIS_SKELETON
