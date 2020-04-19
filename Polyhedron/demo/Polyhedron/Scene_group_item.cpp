#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <QDebug>

using namespace CGAL::Three;
Scene_group_item::Scene_group_item(QString name)
{
    this->name_ = name;
    expanded = true;
    already_drawn = false;
    scene = Three::scene();
}

bool Scene_group_item::isFinite() const
{
  Q_FOREACH(Scene_interface::Item_id id, children)
    if(!getChild(id)->isFinite()){      return false;
    }
  return true;
}

bool Scene_group_item::isEmpty() const {
  Q_FOREACH(Scene_interface::Item_id id, children)
    if(!getChild(id)->isEmpty()){
      return false;
    }
  return true;
}

Scene_group_item::Bbox Scene_group_item::bbox() const
{
 Scene_item* first_non_empty = nullptr;
 Q_FOREACH(Scene_interface::Item_id id, children)
   if(!getChild(id)->isEmpty())
   {
     first_non_empty = getChild(id);
   }

 if(first_non_empty)
 {
   Bbox b =first_non_empty->bbox();
   Q_FOREACH(Scene_interface::Item_id id, children)
     b+=getChild(id)->bbox();
   return b;
 }
 return Bbox(0,0,0,0,0,0);
}


bool Scene_group_item::supportsRenderingMode(RenderingMode m) const {

  Q_FOREACH(Scene_interface::Item_id id, children)
    if(!getChild(id)->supportsRenderingMode(m))
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
    if(!children.contains(scene->item_id(new_item)))
    {
        children.append(scene->item_id(new_item));
        update_group_number(new_item, has_group+1);
    }

}

void Scene_group_item::addChild(Scene_interface::Item_id new_id)
{
  if(!children.contains(new_id))
  {
    children.append(new_id);
    update_group_number(getChild(new_id), has_group+1);
  }
}

void Scene_group_item::update_group_number(Scene_item * new_item, int n)
{

    Scene_group_item* group =
            qobject_cast<Scene_group_item*>(new_item);
    if(group)
      Q_FOREACH(Scene_interface::Item_id id, group->getChildren()){

        update_group_number(getChild(id),n+1);
      }
    new_item->has_group = n;
}
void Scene_group_item::setColor(QColor c)
{
  Scene_item::setColor(c);
  Q_FOREACH(Scene_interface::Item_id id, children)
  {
    getChild(id)->setColor(c);
  }
}

void Scene_group_item::setRenderingMode(RenderingMode m)
{
  Scene_item::setRenderingMode(m);
  Q_FOREACH(Scene_interface::Item_id id, children)
  {
    Scene_item* child = getChild(id);
    if(child->supportsRenderingMode(m))
      child->setRenderingMode(m);
  }
}

void Scene_group_item::setVisible(bool b)
{
  Scene_item::setVisible(b);
  Q_FOREACH(Scene_interface::Item_id id, children)
  {
    Scene_item* child = getChild(id);
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

void Scene_group_item::draw(CGAL::Three::Viewer_interface* ) const  {

}

void Scene_group_item::drawEdges(CGAL::Three::Viewer_interface* ) const
{

}

void Scene_group_item::drawPoints(CGAL::Three::Viewer_interface* ) const
{

}

void Scene_group_item::renderChildren(Viewer_interface *viewer,
                            QMap<float, int>& picked_item_IDs,
                            const QPoint& picked_pixel,
                            bool with_names)
{

  Q_FOREACH(Scene_interface::Item_id id, children){
    if(with_names) {
      viewer->glClearDepthf(1.0f);
      viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    if(id == scene->mainSelectionIndex()|| scene->selectionIndices().contains(id))
    {
      getChild(id)->selection_changed(true);
    }
    else
    {
      getChild(id)->selection_changed(false);
    }
    if(getChild(id)->visible() &&
       (getChild(id)->renderingMode() == Flat ||
        getChild(id)->renderingMode() == FlatPlusEdges ||
        getChild(id)->renderingMode() == Gouraud ||
        getChild(id)->renderingMode() == GouraudPlusEdges))
    {
      getChild(id)->draw(viewer);
    }

    if(getChild(id)->visible() &&
       (getChild(id)->renderingMode() == FlatPlusEdges
        || getChild(id)->renderingMode() == Wireframe
        || getChild(id)->renderingMode() == PointsPlusNormals
        || getChild(id)->renderingMode() == GouraudPlusEdges))
    {
      getChild(id)->drawEdges(viewer);
    }

    if(getChild(id)->visible() &&
       (getChild(id)->renderingMode() == Points  ||
       (getChild(id)->renderingMode() == PointsPlusNormals)  ||
       (getChild(id)->renderingMode() == ShadedPoints)))
    {
      getChild(id)->drawPoints(viewer);
    }
    if(with_names) {
      //    read depth buffer at pick location;
      float depth = 1.0;
      viewer->glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
      if (depth != 1.0)
      {
        //add object to list of picked objects;
        picked_item_IDs[depth] = id;
      }
    }
    CGAL::Three::Scene_group_item* group =
        qobject_cast<CGAL::Three::Scene_group_item*>(getChild(id));
    if(group)
      group->renderChildren(viewer, picked_item_IDs, picked_pixel, with_names);
  }
}

void Scene_group_item::lockChild(Scene_item *child)
{
  lockChild(scene->item_id(child));
}

void Scene_group_item::lockChild(Scene_interface::Item_id id)
{
  if(!children.contains(id))
    return;
  getChild(id)->setProperty("lock", true);
}

void Scene_group_item::unlockChild(Scene_interface::Item_id id)
{
  if(!children.contains(id))
       return;
  getChild(id)->setProperty("lock", false);
}
void Scene_group_item::unlockChild(Scene_item *child)
{
  unlockChild(scene->item_id(child));
}
bool Scene_group_item::isChildLocked(Scene_interface::Item_id id)
{
  if(!children.contains(id)
     || (!getChild(id)->property("lock").toBool()) )
     return false;
   return true;
 }
bool Scene_group_item::isChildLocked(Scene_item *child)
{
  return isChildLocked(scene->item_id(child));
}

void Scene_group_item::setAlpha(int )
{

  Q_FOREACH(Scene_interface::Item_id id, children)
  {
    scene->item(id)->setAlpha(static_cast<int>(alpha()*255));
  }
}
