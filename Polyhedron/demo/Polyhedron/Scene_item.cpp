#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_group_item.h>
#include <QMenu>
#include <QWidgetAction>
#include <QApplication>
#include <QMutex>
#include <iostream>
#include <QDebug>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>

struct D{
  D(Scene_item* item):
    has_group(0),
    name_("unnamed"),
    color_(QColor(100, 100, 255)),
    visible_(true),
    parent_group(0),
    is_selected(false),
    was_selected(false),
    rendering_mode(FlatPlusEdges),
    defaultContextMenu(NULL),
    is_locked(false),
    is_reading(0),
    item(item)
  {
  }

  int has_group;

  //The name of the item.
  QString name_;
  //The color of the item.
  QColor color_;
  //The visibility of the item.
  bool visible_;
  //The parent group, or 0 if the item is not in a group.
  Scene_group_item* parent_group;
  //Specifies if the item is currently selected.
  bool is_selected;
  //Last selection state
  bool was_selected;
  RenderingMode rendering_mode;
  //The default context menu.
  QMenu* defaultContextMenu;
  //when this is true, all the operations and options of this item are disabled.
  //Used in multithreading context.
  bool is_locked;
  int is_reading;
  QMutex mutex;
  int cur_id;
  Scene_item* item;

}; //end D

const QColor Scene_item::defaultColor()
{
  return QColor(100, 100, 255);
}

QColor Scene_item::selectionColor()
{
  QColor c = color();
  return (QColor::fromHsv(
            c.hue()<270 && c.hue() > 90 ? 0
                                        : 180,
            255, 255));
}

namespace CT = CGAL::Three;
CGAL::Three::Scene_item::Scene_item()
  : d(new D(this)){}

CGAL::Three::Scene_item::~Scene_item() {
  if(d && d->defaultContextMenu != NULL)
    d->defaultContextMenu->deleteLater();
}


void CGAL::Three::Scene_item::itemAboutToBeDestroyed(CGAL::Three::Scene_item* item) {
    if(this == item)
    {
      Q_EMIT aboutToBeDestroyed();
    }
}


QString modeName(RenderingMode mode) {
    switch(mode)
    {
    case Points:
        return QObject::tr("points");
    case ShadedPoints:
        return QObject::tr("shaded points");
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
    case ShadedPoints:
      return SLOT(setShadedPointsMode());
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
QString CGAL::Three::Scene_item::renderingModeName() const
{
    return modeName(renderingMode());
}
QMenu* CGAL::Three::Scene_item::contextMenu()
{
    if(d->defaultContextMenu) {
        d->defaultContextMenu->setTitle(name());
        return d->defaultContextMenu;
    }

    this->moveToThread(QApplication::instance()->thread());
    d->defaultContextMenu = new QMenu(name());
    for(unsigned int mode = 0; mode < NumberOfRenderingMode;
        ++mode)
    {
        if(!supportsRenderingMode(RenderingMode(mode))) continue;
        QString mName = modeName(RenderingMode(mode));
        d->defaultContextMenu->addAction(tr("Set %1 Mode")
                                      .arg(mName),
                                      this,
                                      slotName(RenderingMode(mode)));
    }
    return d->defaultContextMenu;
}

CGAL::Three::Scene_group_item* CGAL::Three::Scene_item::parentGroup() const {
  return d->has_group > 0 ? d->parent_group : NULL;
}

void CGAL::Three::Scene_item::
moveToGroup(CGAL::Three::Scene_group_item* group) {
  d->parent_group = group;
  if(group)
    d->has_group = group->hasGroup() + 1;
  else
    d->has_group = 0;
}

void CGAL::Three::Scene_item::invalidate(Gl_data_names) {}

void CGAL::Three::Scene_item::selection_changed(bool b)
{
  if(b && !d->was_selected)
  {
    setSelected(true);
    d->was_selected = true;
    redraw();
    QTimer::singleShot(500, this, [this]{setSelected(false); redraw();});
  }
  else{
    if(!b)
      setSelected(false);
    d->was_selected = b;
  }

}
void CGAL::Three::Scene_item::setVisible(bool b)
{
  d->visible_ = b;
  if(b)
    itemVisibilityChanged();
}


void CGAL::Three::Scene_item::select(double /*orig_x*/,
                        double /*orig_y*/,
                        double /*orig_z*/,
                        double /*dir_x*/,
                        double /*dir_y*/,
                        double /*dir_z*/)
{
}


CGAL::Three::Scene_item::Header_data CGAL::Three::Scene_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  return data;
}

QString CGAL::Three::Scene_item::computeStats(int )
{
  return QString();
}

#include <CGAL/double.h>




float Scene_item::alpha() const
{
  return 1.0f;
}


void Scene_item::writing(){
  d->mutex.lock();
  d->is_locked = true;
  d->mutex.unlock();
  itemChanged();
}
bool Scene_item::isWriting()const{
  d->mutex.lock();
  bool res = d->is_locked;
  d->mutex.unlock();
  return res;
}
void Scene_item::doneWriting() {
  d->mutex.lock();
  d->is_locked = false;
  d->mutex.unlock();
  itemChanged();
}
int Scene_item::isReading()const{
  d->mutex.lock();
  bool res = d->is_reading;
  d->mutex.unlock();
  return res;
}


void Scene_item::reading(){
  d->mutex.lock();
  ++d->is_reading;
  d->mutex.unlock();
  itemChanged();
}


void Scene_item::doneReading() {
  d->mutex.lock();
  --d->is_reading;
  d->mutex.unlock();
  itemChanged();
}

void Scene_item::setRenderingMode(RenderingMode m) {
  if (supportsRenderingMode(m))
    d->rendering_mode = m;
  Q_EMIT redraw();
}

int            Scene_item::hasGroup() const      { return d->has_group; }
 QColor        Scene_item::color() const         { return d->color_;         }
 QString       Scene_item::name() const          { return d->name_;          }
 bool          Scene_item::visible() const       { return d->visible_;       }
 RenderingMode Scene_item::renderingMode() const { return d->rendering_mode; }
bool           Scene_item::isSelected() const    { return d->is_selected; }
int            Scene_item::getId() const         { return d->cur_id; }
void           Scene_item::setSelected(bool b)   { d->is_selected = b; }
void           Scene_item::setColor(QColor c)    { d->color_ = c;}
void           Scene_item::setName(QString n)    { d->name_ = n; }
void           Scene_item::setHasGroup(int hg)   { d->has_group = hg; }
void           Scene_item::setId(int id)         {d->cur_id = id; }
