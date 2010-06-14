#include <CGAL_demo/Scene_item.h>

const QColor Scene_item::defaultColor = QColor(100, 100, 255);

Scene_item::~Scene_item() {}

void Scene_item::itemAboutToBeDestroyed(Scene_item* item) {
  if(this == item)
    emit aboutToBeDestroyed();
}

// Rendering mode as a human readable string
QString Scene_item::renderingModeName() const
{
  switch(renderingMode()) 
  {
    case Points:
      return tr("points");
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
