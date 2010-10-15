#include "Scene_item.h"

const QColor Scene_item::defaultColor = QColor(100, 100, 255);

Scene_item::~Scene_item() {}

void Scene_item::itemAboutToBeDestroyed(Scene_item* item) {
  if(this == item)
    emit aboutToBeDestroyed();
}

// Rendering mode list as a human readable string
QString Scene_item::renderingModeNameList()
{
    return tr("Rendering mode (points/points+normals/splatting/wireframe/flat/flat+edges/Gouraud)");
}

// Rendering mode as a human readable string
QString Scene_item::renderingModeName() const
{
  switch(renderingMode()) 
  {
    case Points:
      return tr("points");
    case PointsPlusNormals:
      return tr("pts+normals");
    case Splatting:
      return tr("splats");
    case Wireframe:
      return tr("wire");
    case Flat:
      return tr("flat");
    case FlatPlusEdges:
      return tr("flat+edges");
    case Gouraud:
      return tr("Gouraud");
    default:
      Q_ASSERT(false);
      return tr("unknown");
  }
} 

#include "Scene_item.moc"

