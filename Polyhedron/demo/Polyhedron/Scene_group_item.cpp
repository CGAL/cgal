#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <QDebug>

using namespace CGAL::Three;
Scene_group_item::Scene_group_item(QString name, int nb_vbos, int nb_vaos )
    :  Scene_item(nb_vbos, nb_vaos)
    , scene(NULL)
{
    this->name_ = name;
    expanded = true;
    already_drawn = false;
}

bool Scene_group_item::isFinite() const
{
    Q_FOREACH(Scene_item *item, children)
        if(!item->isFinite()){
            return false;
        }
    return true;
}

bool Scene_group_item::isEmpty() const {
    Q_FOREACH(Scene_item *item, children)
        if(!item->isEmpty()){
            return true;
        }
    return true;
}

Scene_group_item::Bbox Scene_group_item::bbox() const
{
    return Bbox(0, 0, 0, 0, 0,0);
}


bool Scene_group_item::supportsRenderingMode(RenderingMode m) const {

    Q_FOREACH(Scene_item* item, children)
        if(!item->supportsRenderingMode(m))
            return false;
    return true;

}

QString Scene_group_item::toolTip() const {
    QString str =
            QObject::tr( "<p>Number of children: %1<br />").arg(children.size());
    str+="</p>";
    return str;
}

void Scene_group_item::addChild(Scene_item* new_item)
{  
    if(!children.contains(new_item))
    {
        children.append(new_item);
        update_group_number(new_item, has_group+1);
    }

}

void Scene_group_item::update_group_number(Scene_item * new_item, int n)
{

    Scene_group_item* group =
            qobject_cast<Scene_group_item*>(new_item);
    if(group)
      Q_FOREACH(Scene_item* child, group->getChildren())
          update_group_number(child,n+1);
    new_item->has_group = n;
}
void Scene_group_item::setColor(QColor c)
{
  Scene_item::setColor(c);
  Q_FOREACH(Scene_item* child, children)
  {
    child->setColor(c);
  }
}

void Scene_group_item::setRenderingMode(RenderingMode m)
{
  Scene_item::setRenderingMode(m);
  Q_FOREACH(Scene_item* child, children)
  {
    if(child->supportsRenderingMode(m))
      child->setRenderingMode(m);
  }
}

void Scene_group_item::setVisible(bool b)
{
  Scene_item::setVisible(b);
  Q_FOREACH(Scene_item* child, children)
  {
    child->setVisible(b);
    child->itemChanged();
  }
    Q_EMIT itemChanged();
}

bool Scene_group_item::isExpanded() const
{
  return expanded;
}

void Scene_group_item::setExpanded(bool b)
{
    expanded = b;
}

void Scene_group_item::moveDown(int i)
{
    children.move(i, i+1);
}

void Scene_group_item::moveUp(int i)
{
    children.move(i, i-1);
}

void Scene_group_item::draw(CGAL::Three::Viewer_interface* viewer) const  {
  if(viewer->inDrawWithNames() || already_drawn ) return;
  Q_FOREACH(Scene_item* child, children) {
    if(!child->visible()) continue;
    switch(child->renderingMode()) {
    case Flat:
    case FlatPlusEdges:
    case Gouraud:
      child->draw(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case FlatPlusEdges:
    case Wireframe:
    case PointsPlusNormals:
      child->drawEdges(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case Points:
    case ShadedPoints:
    case PointsPlusNormals:
      child->drawPoints(viewer); break;
    default: break;
    }
  }
  already_drawn = true;
}

void Scene_group_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(viewer->inDrawWithNames() || already_drawn ) return;
  Q_FOREACH(Scene_item* child, children) {
    if(!child->visible()) continue;
    switch(child->renderingMode()) {
    case FlatPlusEdges:
    case Wireframe:
    case PointsPlusNormals:
      child->drawEdges(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case Flat:
    case FlatPlusEdges:
    case Gouraud:
      child->draw(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case Points:
    case ShadedPoints:
    case PointsPlusNormals:
      child->drawPoints(viewer); break;
    default: break;
    }
  }
  already_drawn = true;
}

void Scene_group_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
  if(viewer->inDrawWithNames() || already_drawn ) return;
  Q_FOREACH(Scene_item* child, children) {
    if(!child->visible()) continue;
    switch(child->renderingMode()) {
    case Points:
    case ShadedPoints:
    case PointsPlusNormals:
      child->drawPoints(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case Flat:
    case FlatPlusEdges:
    case Gouraud:
      child->draw(viewer); break;
    default: break;
    }
    switch(child->renderingMode()) {
    case FlatPlusEdges:
    case Wireframe:
    case PointsPlusNormals:
      child->drawEdges(viewer); break;
    default: break;
    }
  }
  already_drawn = true;
}

void Scene_group_item::lockChild(Scene_item *child)
{
  if(!children.contains(child))
    return;
  child->setProperty("lock", true);
}
void Scene_group_item::unlockChild(Scene_item *child)
{
  if(!children.contains(child))
    return;
  child->setProperty("lock", false);
}
bool Scene_group_item::isChildLocked(Scene_item *child)
{
  if(!children.contains(child)
     || (!child->property("lock").toBool()) )
    return false;
  return true;
}
