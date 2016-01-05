
#include <CGAL/Three/Scene_group_item.h>
#include <QDebug>

using namespace CGAL::Three;
Scene_group_item::Scene_group_item(QString name)
    :  Scene_item(0,0)
{
    this->name_ = name;
    expanded = true;
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
    return !children.isEmpty();

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
        add_group_number(new_item);
    }

}

void Scene_group_item::add_group_number(Scene_item * new_item)
{

    Scene_group_item* group =
            qobject_cast<Scene_group_item*>(new_item);
    if(group)
      Q_FOREACH(Scene_item* child, group->getChildren())
          add_group_number(child);
    new_item->has_group++;
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
    child->setRenderingMode(m);
  }
}

void Scene_group_item::setVisible(bool b)
{
  Scene_item::setVisible(b);
  Q_FOREACH(Scene_item* child, children)
  {
    child->setVisible(b);
  }
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
