#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <QDebug>
#include <QSlider>

using namespace CGAL::Three;
Scene_group_item::Scene_group_item(QString name, Scene_interface *scene)
  : scene(scene)
{
  setName(name);
  expanded = true;
}

bool Scene_group_item::isFinite() const
{
  Q_FOREACH(Scene_interface::Item_id id, children)
    if(!getChild(id)->isFinite()){
      return false;
    }
  return true;
}

bool Scene_group_item::isEmpty() const {
  Q_FOREACH(Scene_interface::Item_id id, children)
    if(!getChild(id)->isEmpty()){
      return true;
    }
  return true;
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

void Scene_group_item::addChild(Scene_interface::Item_id new_id)
{
  if(!children.contains(new_id))
  {
    children.append(new_id);
    update_group_number(getChild(new_id), hasGroup()+1);
  }
}

void Scene_group_item::addChild(Scene_item* new_item)
{  
  addChild(scene->item_id(new_item));
}

void Scene_group_item::update_group_number(Scene_item * new_item, int n)
{

  Scene_group_item* group =
      qobject_cast<Scene_group_item*>(new_item);
  if(group)
    Q_FOREACH(Scene_interface::Item_id id, group->getChildrenIds())
    {
      update_group_number(getChild(id),n+1);
    }
  new_item->setHasGroup(n);
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
    if(getChild(id)->supportsRenderingMode(m))
      getChild(id)->setRenderingMode(m);
  }
}

void Scene_group_item::setVisible(bool b)
{
  Scene_item::setVisible(b);
  Q_FOREACH(Scene_interface::Item_id id, children)
  {
    getChild(id)->setVisible(b);
    getChild(id)->itemChanged();
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

void Scene_group_item::draw(CGAL::Three::Viewer_interface* viewer, int pass , bool is_writing, QOpenGLFramebufferObject *fbo) const
{
  if(!isInit())
    initGL();
  Q_FOREACH(Scene_interface::Item_id id, children){
    if(getChild(id)->visible() &&
       (getChild(id)->renderingMode() == Flat ||
        getChild(id)->renderingMode() == FlatPlusEdges ||
        getChild(id)->renderingMode() == Gouraud))
    {
      getChild(id)->draw(viewer, pass, is_writing, fbo);
    }
  }
}

void Scene_group_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit())
    initGL();
  Q_FOREACH(Scene_interface::Item_id id, children){
    if(getChild(id)->visible() &&
       (getChild(id)->renderingMode() == FlatPlusEdges
        || getChild(id)->renderingMode() == Wireframe
        || getChild(id)->renderingMode() == PointsPlusNormals))
    {
      getChild(id)->drawEdges(viewer);
    }
  }
}

void Scene_group_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit())
    initGL();
  Q_FOREACH(Scene_interface::Item_id id, children){
    if(getChild(id)->visible())
    {

      if(getChild(id)->renderingMode() == Points  ||
         (getChild(id)->renderingMode() == PointsPlusNormals)  ||
         (getChild(id)->renderingMode() == ShadedPoints))
      {
        getChild(id)->drawPoints(viewer);
      }
    }
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

void Scene_group_item::initGL()const
{
  Scene_item_rendering_helper::initGL();
  Scene_group_item* ncthis = const_cast<Scene_group_item*>(this);
  connect(ncthis->alphaSlider(), &QSlider::valueChanged,
          ncthis, &Scene_group_item::setAlpha);
}

void Scene_group_item::compute_bbox() const
{
  setBbox(Bbox(0, 0, 0, 0, 0, 0));
}
