#ifndef SCENE_MCF_ITEM_H
#define SCENE_MCF_ITEM_H



#if defined( scene_mcf_item_EXPORTS )
#  define SCENE_MCF_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_MCF_ITEM_EXPORT Q_DECL_IMPORT
#endif


#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include "SMesh_type.h"
typedef SMesh Face_graph;
typedef CGAL::Mean_curvature_flow_skeletonization<Face_graph>      Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton Skeleton;

struct Scene_mcf_item_priv;
// This class represents a polyhedron in the OpenGL scene
class SCENE_MCF_ITEM_EXPORT Scene_mcf_item
        : public CGAL::Three::Scene_group_item {
    Q_OBJECT
public:

  Scene_mcf_item(Face_graph* graph,
                    Scene_interface::Item_id id,
                    QString name);

  ~Scene_mcf_item();

public:
  Mean_curvature_skeleton* mcs;
  Face_graph* meso_skeleton; // a copy of the meso_skeleton that is displayed
  Face_graph* input_triangle_mesh;

  int fixedPointsItemIndex;
  int nonFixedPointsItemIndex;
  int poleLinesItemIndex;
  int skeletonItemIndex;
  int contractedItemIndex;
  int InputMeshItemIndex;

  Skeleton skeleton_curve;
}; // end class Scene_mcf_item
#endif // SCENE_MCF_ITEM_H
