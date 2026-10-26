#include "Scene_mcf_item.h"

Scene_mcf_item::Scene_mcf_item(Face_graph* graph,
                  Scene_interface::Item_id id,
                  QString name)
{
  mcs = NULL;
  meso_skeleton = NULL;
  input_triangle_mesh = graph;
  fixedPointsItemIndex = -1;
  skeletonItemIndex = -1;
  nonFixedPointsItemIndex = -1;
  poleLinesItemIndex = -1;
  contractedItemIndex = -1;
  InputMeshItemIndex = id;
  meso_skeleton = NULL;
  this->setName(name);
}

Scene_mcf_item::~Scene_mcf_item()
{
  delete mcs;
  delete meso_skeleton;
}
