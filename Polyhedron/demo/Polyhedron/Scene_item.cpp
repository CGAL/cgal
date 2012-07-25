#include "Scene_item.h"
#include "Scene_interface.h"
#include <QMenu>

const QColor Scene_item::defaultColor = QColor(100, 100, 255);

Scene_item::~Scene_item() {
  delete defaultContextMenu;
}

void Scene_item::itemAboutToBeDestroyed(Scene_item* item) {
  if(this == item)
    emit aboutToBeDestroyed();
}


QString modeName(RenderingMode mode) {
  switch(mode) 
  {
    case Points:
      return QObject::tr("points");
    case Wireframe:
      return QObject::tr("wire");
    case Flat:
      return QObject::tr("flat");
    case FlatPlusEdges:
      return QObject::tr("flat+edges");
    case Gouraud:
      return QObject::tr("Gouraud");
    case PointsPlusNormals:
      return QObject::tr("pts+normals");
    default:
      Q_ASSERT(false);
      return QObject::tr("unknown");
  }
}

const char* slotName(RenderingMode mode) {
  switch(mode) 
  {
    case Points:
      return SLOT(setPointsMode());
    case Wireframe:
      return SLOT(setWireframeMode());
    case Flat:
      return SLOT(setFlatMode());
    case FlatPlusEdges:
      return SLOT(setFlatPlusEdgesMode());
    case Gouraud:
      return SLOT(setGouraudMode());
    case PointsPlusNormals:
      return SLOT(setPointsPlusNormalsMode());
    default:
      Q_ASSERT(false);
      return "";
  }
}

// Rendering mode as a human readable string
QString Scene_item::renderingModeName() const
{
  return modeName(renderingMode());
} 

QMenu* Scene_item::contextMenu()
{
  if(defaultContextMenu) {
    defaultContextMenu->setTitle(name());
    return defaultContextMenu;
  }

  defaultContextMenu = new QMenu(name());
  // defaultContextMenu->addAction(name());
  // defaultContextMenu->addSeparator();
  // QMenu* modeMenu = new QMenu(QObject::tr("Rendering mode"),
  //                             defaultContextMenu);
  for(unsigned int mode = 0; mode < NumberOfRenderingMode;
      ++mode) 
  {
    if(!supportsRenderingMode(RenderingMode(mode))) continue;
    QString mName = modeName(RenderingMode(mode));
    QAction* action = 
      defaultContextMenu->addAction(tr("Set %1 mode")
                                    .arg(mName),
                                    this,
                                    slotName(RenderingMode(mode)));
    QObject::connect(action, SIGNAL(triggered()),
                     this, SIGNAL(itemChanged()));
  }
  // defaultContextMenu->addAction(modeMenu->menuAction());
  return defaultContextMenu;
}

void Scene_item::changed() {
  // emit itemChanged();
}

void Scene_item::select(double /*orig_x*/,
                        double /*orig_y*/,
                        double /*orig_z*/,
                        double /*dir_x*/,
                        double /*dir_y*/,
                        double /*dir_z*/)
{
}

#include "Scene_item.moc"

