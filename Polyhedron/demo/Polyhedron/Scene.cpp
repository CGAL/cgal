#include "GlSplat/GlSplat.h"



#include "config.h"
#include "Scene.h"
#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_print_interface_item.h>

#include <QObject>
#include <QMetaObject>
#include <QString>
#include <QGLWidget>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QColorDialog>
#include <QApplication>
#include <QPointer>
#include <QList>
#include <QAbstractProxyModel>
#include <QMimeData>

GlSplat::SplatRenderer* Scene::ms_splatting = 0;
int Scene::ms_splattingCounter = 0;
GlSplat::SplatRenderer* Scene::splatting()
{
    assert(ms_splatting!=0 && "A Scene object must be created before requesting the splatting object");
    return ms_splatting;
}

Scene::Scene(QObject* parent)
    : QStandardItemModel(parent),
      selected_item(-1),
      item_A(-1),
      item_B(-1)
{

    connect(this, SIGNAL(selectionRay(double, double, double,
                                      double, double, double)),
            this, SLOT(setSelectionRay(double, double, double,
                                       double, double, double)));
    connect(this, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
            this, SLOT(s_itemAboutToBeDestroyed(CGAL::Three::Scene_item*)));
    if(ms_splatting==0)
        ms_splatting  = new GlSplat::SplatRenderer();
    ms_splattingCounter++;
    picked = false;
    gl_init = false;

}
Scene::Item_id
Scene::addItem(CGAL::Three::Scene_item* item)
{
    Bbox bbox_before = bbox();
    m_entries.push_back(item);
    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    connect(item, SIGNAL(redraw()),
            this, SLOT(callDraw()));
    if(item->isFinite()
            && !item->isEmpty()
            && bbox_before + item->bbox() != bbox_before
            )
    {
        Q_EMIT updated_bbox();
    }
    QList<QStandardItem*> list;
    for(int i=0; i<5; i++)
    {
        list<<new QStandardItem();
        list.at(i)->setEditable(false);
    }
    invisibleRootItem()->appendRow(list);
    for(int i=0; i<5; i++){
        index_map[list.at(i)->index()] = m_entries.size() -1;
    }
    Q_EMIT updated();
    Item_id id = m_entries.size() - 1;
    Q_EMIT newItem(id);
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item);
    if(group)
        addGroup(group);
    return id;
}

CGAL::Three::Scene_item*
Scene::replaceItem(Scene::Item_id index, CGAL::Three::Scene_item* item, bool emit_item_about_to_be_destroyed)
{
    if(index < 0 || index >= m_entries.size())
        return 0;

    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(m_entries[index]);
    if(group)
    {
      QList<int> children;
      Q_FOREACH(CGAL::Three::Scene_item* child, group->getChildren())
      {
        group->unlockChild(child);
        children << item_id(child);
      }
      erase(children);
    }
    std::swap(m_entries[index], item);
    if ( item->isFinite() && !item->isEmpty() &&
         m_entries[index]->isFinite() && !m_entries[index]->isEmpty() &&
         item->bbox()!=m_entries[index]->bbox() )
    {
    Q_EMIT updated_bbox();
    }

    if(emit_item_about_to_be_destroyed) {
      Q_EMIT itemAboutToBeDestroyed(item);
    }

    Q_EMIT updated();
    group =
            qobject_cast<CGAL::Three::Scene_group_item*>(m_entries[index]);
    if(group)
    {
        addGroup(group);
    }
    itemChanged(index);
    Q_EMIT restoreCollapsedState();
    redraw_model();
    Q_EMIT selectionChanged(index);
    return item;
}

Scene::Item_id
Scene::erase(Scene::Item_id index)
{
  CGAL::Three::Scene_item* item = m_entries[index];
  if(item->parentGroup()
     && item->parentGroup()->isChildLocked(item))
    return -1;
  //clears the Scene_view
    clear();
    index_map.clear();
    if(index < 0 || index >= m_entries.size())
        return -1;
  if(item->parentGroup())
    item->parentGroup()->removeChild(item);

  //removes the item from all groups that contain it
  Q_EMIT itemAboutToBeDestroyed(item);
    item->deleteLater();
    m_entries.removeAll(item);
    selected_item = -1;
    //re-creates the Scene_view
    Q_FOREACH(Scene_item* item, m_entries)
    {
        organize_items(item, invisibleRootItem(), 0);
    }
    QStandardItemModel::beginResetModel();
    Q_EMIT updated();
    QStandardItemModel::endResetModel();
    Q_EMIT restoreCollapsedState();
    if(--index >= 0)
      return index;
    if(!m_entries.isEmpty())
      return 0;
    return -1;

}

int
Scene::erase(QList<int> indices)
{
  if(indices.empty())
    return -1;
  QList<CGAL::Three::Scene_item*> to_be_removed;
  int max_index = -1;
  Q_FOREACH(int index, indices) {
    if(index < 0 || index >= m_entries.size())
      continue;

    max_index = (std::max)(max_index, index);
    CGAL::Three::Scene_item* item = m_entries[index];
    if(item->parentGroup()
       && item->parentGroup()->isChildLocked(item))
      if(!indices.contains(item_id(item->parentGroup())))
        continue;
    if(!to_be_removed.contains(item))
      to_be_removed.push_back(item);
  }

  Q_FOREACH(Scene_item* item, to_be_removed) {
      if(item->parentGroup())
        item->parentGroup()->removeChild(item);
    Q_EMIT itemAboutToBeDestroyed(item);
    item->deleteLater();
    m_entries.removeAll(item);
  }
  clear();
  index_map.clear();
  selected_item = -1;
  Q_FOREACH(Scene_item* item, m_entries)
  {
    organize_items(item, invisibleRootItem(), 0);
  }
  QStandardItemModel::beginResetModel();
  Q_EMIT updated();
  QStandardItemModel::endResetModel();
  Q_EMIT restoreCollapsedState();

  int index = max_index + 1 - indices.size();
  if(index >= m_entries.size()) {
    index = m_entries.size() - 1;
  }
  if(index >= 0)
    return index;
  if(!m_entries.isEmpty())
    return 0;
  return -1;

}

void Scene::remove_item_from_groups(Scene_item* item)
{
    CGAL::Three::Scene_group_item* group = item->parentGroup();
    if(group)
    {
        group->removeChild(item);
    }
}
Scene::~Scene()
{
    Q_FOREACH(CGAL::Three::Scene_item* item_ptr, m_entries)
    {
         item_ptr->deleteLater();
    }
    m_entries.clear();

    if((--ms_splattingCounter)==0)
        delete ms_splatting;
}

CGAL::Three::Scene_item*
Scene::item(Item_id index) const
{
    return m_entries.value(index); // QList::value checks bounds
}

Scene::Item_id 
Scene::item_id(CGAL::Three::Scene_item* scene_item) const
{
    return m_entries.indexOf(scene_item);
}

int
Scene::numberOfEntries() const
{
    return m_entries.size();
}

// Duplicate a scene item.
// Return the ID of the new item (-1 on error).
Scene::Item_id
Scene::duplicate(Item_id index)
{
    if(index < 0 || index >= m_entries.size())
        return -1;

    const CGAL::Three::Scene_item* item = m_entries[index];
    CGAL::Three::Scene_item* new_item = item->clone();
    if(new_item) {
        new_item->setName(tr("%1 (copy)").arg(item->name()));
        new_item->setColor(item->color());
        new_item->setVisible(item->visible());
        addItem(new_item);
        return m_entries.size() - 1;
    }
    else
        return -1;
}

void Scene::initializeGL()
{
    ms_splatting->init();

    //Setting the light options

    // Create light components
    GLfloat ambientLight[] = { 0.4f, 0.4f, 0.4f, 1.0f };
    GLfloat diffuseLight[] = { 1.0f, 1.0f, 1.0, 1.0f };
    GLfloat specularLight[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat position[] = { 0.0f, 0.0f, 1.0f, 1.0f };

    // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    gl_init = true;
}

void Scene::s_itemAboutToBeDestroyed(CGAL::Three::Scene_item *rmv_itm)
{
 Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
 {
   if(item == rmv_itm)
     item->itemAboutToBeDestroyed(item);
 }
}
bool
Scene::keyPressEvent(QKeyEvent* e){
    bool res=false;
    for (QList<int>::iterator it=selected_items_list.begin(),endit=selected_items_list.end();
         it!=endit;++it)
    {
        CGAL::Three::Scene_item* item=m_entries[*it];
        res |= item->keyPressEvent(e);
    }
    return res;
}

void
Scene::draw(CGAL::Three::Viewer_interface* viewer)
{
    draw_aux(false, viewer);
}
void
Scene::drawWithNames(CGAL::Three::Viewer_interface* viewer)
{
    draw_aux(true, viewer);
}

bool item_should_be_skipped_in_draw(Scene_item* item) {
  if(!item->visible()) return true;
  if(item->has_group == 0) return false;
  Scene_group_item* group = item->parentGroup();
  while(group != 0) {
    if(!group->visible()) return false;
    group = group->parentGroup();
  }
  return true;
}

void 
Scene::draw_aux(bool with_names, CGAL::Three::Viewer_interface* viewer)
{
    QMap<float, int> picked_item_IDs;
    if(with_names)
      glEnable(GL_DEPTH_TEST);
    if(!ms_splatting->viewer_is_set)
        ms_splatting->setViewer(viewer);
    if(!gl_init)
        initializeGL();
    // Flat/Gouraud OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        CGAL::Three::Scene_item& item = *m_entries[index];
        if(index == selected_item || selected_items_list.contains(index))
        {
            item.selection_changed(true);
        }
        else

        {
            item.selection_changed(false);
        }
        if(!with_names && item_should_be_skipped_in_draw(&item)) continue;
        if(item.visible())
        {
            if(item.renderingMode() == Flat || item.renderingMode() == FlatPlusEdges || item.renderingMode() == Gouraud)
            {
                if(with_names) {
                    glClearDepth(1.0);
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                }
                viewer->glEnable(GL_LIGHTING);
                viewer->glPointSize(2.f);
                viewer->glLineWidth(1.0f);
                if(item.renderingMode() == Gouraud)
                    viewer->glShadeModel(GL_SMOOTH);
                else
                    viewer->glShadeModel(GL_FLAT);
                if(viewer)
                    item.draw(viewer);
                else
                    item.draw();

                if(with_names) {

                    //    read depth buffer at pick location;
                    float depth = 1.0;
                    glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
                    if (depth != 1.0)
                    {
                        //add object to list of picked objects;
                        picked_item_IDs[depth] = index;
                    }
                }
            }
        }
    }
    glDepthFunc(GL_LEQUAL);
    // Wireframe OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        CGAL::Three::Scene_item& item = *m_entries[index];
        if(index == selected_item || selected_items_list.contains(index))
        {
            item.selection_changed(true);
        }
        else
        {
            item.selection_changed(false);
        }

        if(!with_names && item_should_be_skipped_in_draw(&item)) continue;
        if(item.visible())
        {
            if((item.renderingMode() == Wireframe || item.renderingMode() == PointsPlusNormals )
                    && with_names)
            {
                glClearDepth(1.0);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            }
            if((!with_names && item.renderingMode() == FlatPlusEdges )
                    || item.renderingMode() == Wireframe)
            {
                viewer->glDisable(GL_LIGHTING);
                viewer->glPointSize(2.f);
                viewer->glLineWidth(1.0f);

                if(viewer)
                    item.drawEdges(viewer);
                else
                    item.drawEdges();
            }
            else{
                if( item.renderingMode() == PointsPlusNormals ){
                    viewer->glDisable(GL_LIGHTING);
                    viewer->glPointSize(2.f);
                    viewer->glLineWidth(1.0f);
                    if(index == selected_item || selected_items_list.contains(index))
                    {

                        item.selection_changed(true);
                    }
                    else
                    {

                        item.selection_changed(false);
                    }
                    if(viewer)
                        item.drawEdges(viewer);
                    else
                        item.drawEdges();
                }
            }
            if((item.renderingMode() == Wireframe || item.renderingMode() == PointsPlusNormals )
                    && with_names)
            {

                //    read depth buffer at pick location;
                float depth = 1.0;
                glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
                if (depth != 1.0)
                {
                    //add object to list of picked objects;
                    picked_item_IDs[depth] = index;
                }
            }
        }
    }
    // Points OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        CGAL::Three::Scene_item& item = *m_entries[index];
        if(!with_names && item_should_be_skipped_in_draw(&item)) continue;
        if(item.visible())
        {
            if(item.renderingMode() == Points && with_names) {
                glClearDepth(1.0);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            }
            if(item.renderingMode() == Points  ||
                    (!with_names && item.renderingMode() == PointsPlusNormals)  ||
                 (!with_names && item.renderingMode() == ShadedPoints))
            {
                viewer->glDisable(GL_LIGHTING);
                viewer->glPointSize(2.0f);
                viewer->glLineWidth(1.0f);

                if(viewer)
                    item.drawPoints(viewer);
                else
                    item.drawPoints();
            }
            if(item.renderingMode() == Points && with_names) {
                //    read depth buffer at pick location;
                float depth = 1.0;
                glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
                if (depth != 1.0)
                {
                    //add object to list of picked objects;
                    picked_item_IDs[depth] = index;
                }
            }

            if(!with_names)
            {
                glDepthFunc(GL_LESS);
                // Splatting
                if(!with_names && ms_splatting->isSupported())
                {
                    ms_splatting->beginVisibilityPass();
                    for(int index = 0; index < m_entries.size(); ++index)
                    {
                        CGAL::Three::Scene_item& item = *m_entries[index];
                        if(!with_names && item_should_be_skipped_in_draw(&item)) continue;
                        if(item.visible() && item.renderingMode() == Splatting)
                        {

                          if(viewer)
                          {
                             item.drawSplats(viewer);
                          }
                          else
                              item.drawSplats();
                        }

                    }
                    ms_splatting->beginAttributePass();
                    for(int index = 0; index < m_entries.size(); ++index)
                    {  CGAL::Three::Scene_item& item = *m_entries[index];
                        if(item.visible() && item.renderingMode() == Splatting)
                        {
                            viewer->glColor4d(item.color().redF(), item.color().greenF(), item.color().blueF(), item.color().alphaF());
                            if(viewer)
                                item.drawSplats(viewer);
                            else
                                item.drawSplats();
                        }
                    }
                    ms_splatting->finalize();
                }
                else
                    item.drawSplats();
            }
        }
    }
    if(with_names)
    {
        QList<float> depths = picked_item_IDs.keys();
        if(!depths.isEmpty())
        {
            qSort(depths);
            int id = picked_item_IDs[depths.first()];
            setSelectedItemIndex(id);
            viewer->setSelectedName(id);

        }
    }
    if(with_names)
        picked = true;
    else
        picked = false;
    //scrolls the sceneView to the selected item's line.
    if(picked)
    {
        Q_EMIT(itemPicked(index_map.key(mainSelectionIndex())));
    }
    Q_EMIT drawFinished();
}

// workaround for Qt-4.2 (see above)
#undef lighter
QVariant
Scene::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
    {
        return QVariant();
    }

    int id = index_map[index];
    if(id < 0 || id >= m_entries.size())
        return QVariant();
    if(role == ::Qt::ToolTipRole)
    {
        return m_entries[id]->toolTip();
    }
    switch(index.column())
    {
    case ColorColumn:
        if(role == ::Qt::DecorationRole)
            return m_entries.value(id)->color();
        break;
    case NameColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(id)->name();
        if(role == ::Qt::FontRole)
            return m_entries.value(id)->font();
        break;
    case RenderingModeColumn:
        if(role == ::Qt::DisplayRole) {
            return m_entries.value(id)->renderingModeName();
        }
        else if(role == ::Qt::EditRole) {
            return static_cast<int>(m_entries.value(id)->renderingMode());
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignCenter;
        }
        break;
    case ABColumn:
        if(role == ::Qt::DisplayRole) {
            if(id == item_A)
                return "A";
            if(id == item_B)
                return "B";
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignLeft;
        }
        break;
    case VisibleColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(id)->visible();
        break;
    default:
        return QVariant();
    }
    return QVariant();
}

QVariant
Scene::headerData ( int section, ::Qt::Orientation orientation, int role ) const
{
    if(orientation == ::Qt::Horizontal)  {
        if (role == ::Qt::DisplayRole)
        {
            switch(section)
            {
            case NameColumn:
                return tr("Name");
                break;
            case ColorColumn:
                return tr("#");
                break;
            case RenderingModeColumn:
                return tr("Mode");
            case ABColumn:
                return tr("A/B");
                break;
            case VisibleColumn:
                return tr("View");
                break;
            default:
                return QVariant();
            }
        }
        else if(role == ::Qt::ToolTipRole) {
            if(section == RenderingModeColumn) {
                return tr("Rendering mode (points/wireframe/flat/flat+edges/Gouraud)");
            }
            else if(section == ABColumn) {
                return tr("Selection A/Selection B");
            }
        }
    }
    return QStandardItemModel::headerData(section, orientation, role);
}

Qt::ItemFlags
Scene::flags ( const QModelIndex & index ) const
{
    if (index.isValid() && index.column() == NameColumn) {
        return QStandardItemModel::flags(index) | ::Qt::ItemIsEditable;
    }
    else {
        return QStandardItemModel::flags(index);
    }
}

bool
Scene::setData(const QModelIndex &index,
               const QVariant &value,
               int role)
{

    if( role != ::Qt::EditRole || !index.isValid() )
        return false;

    int id = index_map[index];
    if(id < 0 || id >= m_entries.size()){
        return false;
    }

    CGAL::Three::Scene_item* item = m_entries[id];

    if(!item) return false;
    switch(index.column())
    {
    case NameColumn:
        item->setName(value.toString());
    Q_EMIT dataChanged(index, index);
        return true;
        break;
    case ColorColumn:
        item->setColor(value.value<QColor>());
    Q_EMIT dataChanged(index, index);
        return true;
        break;
    case RenderingModeColumn:
    {
        RenderingMode rendering_mode = static_cast<RenderingMode>(value.toInt());
        // Find next supported rendering mode
        while ( ! item->supportsRenderingMode(rendering_mode)
      //          || (rendering_mode==Splatting && !Scene::splatting()->isSupported())
                )
        {
            rendering_mode = static_cast<RenderingMode>( (rendering_mode+1) % NumberOfRenderingMode );
        }
        item->setRenderingMode(rendering_mode);
        QModelIndex nindex = createIndex(m_entries.size()-1,RenderingModeColumn+1);
    Q_EMIT dataChanged(index, nindex);
        return true;
        break;
    }
    case VisibleColumn:
        item->setVisible(value.toBool());
    Q_EMIT dataChanged(index, createIndex(m_entries.size()-1,VisibleColumn+1));
        return true;
    default:
        return false;
    }
    return false;
}

bool Scene::dropMimeData(const QMimeData * /*data*/,
                         Qt::DropAction /*action*/,
                         int /*row*/,
                         int /*column*/,
                         const QModelIndex &parent)
{
    //gets the moving items
    QList<Scene_item*> items;
    QList<int> groups_children;

    //get IDs of all children of selected groups
    Q_FOREACH(int i, selected_items_list)
    {
      CGAL::Three::Scene_group_item* group =
          qobject_cast<CGAL::Three::Scene_group_item*>(item(i));
      if(group)
        Q_FOREACH(Scene_item* child, group->getChildren())
          groups_children << item_id(child);
    }
    // Insure that children of selected groups will not be added twice
    Q_FOREACH(int i, selected_items_list)
    {
      if(!groups_children.contains(i))
      {
        items << item(i);
      }
    }
    //Gets the group at the drop position
    CGAL::Three::Scene_group_item* group = NULL;
    if(parent.isValid())
        group = qobject_cast<CGAL::Three::Scene_group_item*>(this->item(index_map[parent]));
    bool one_contained = false;
    if(group)
    {
      Q_FOREACH(int id, selected_items_list)
        if(group->getChildren().contains(item(id)))
        {
          one_contained = true;
          break;

        }
    }
    //if the drop item is not a group_item or if it already contains the item, then the drop action must be ignored
    if(!group ||one_contained)
    {
      //unless the drop zone is empty, which means the item should be removed from all groups.
      if(!parent.isValid())
      {
        Q_FOREACH(Scene_item* item, items)
        {
          if(item->parentGroup())
          {
            item->parentGroup()->removeChild(item);
          }
        }
        redraw_model();
        return true;
      }
      return false;
    }
    Q_FOREACH(Scene_item* item, items)
      changeGroup(item, group);
    redraw_model();
    return true;
}

void Scene::moveRowUp()
{
    Scene_item* selected_item = item(mainSelectionIndex());
    if(index_map.key(mainSelectionIndex()).row() > 0)
    {
        if(item(mainSelectionIndex())->has_group >0)
        {
            Scene_group_item* group = selected_item->parentGroup();
            if(group)
            {
                int id = group->getChildren().indexOf(selected_item);
                group->moveUp(id);
            }
        }
        else
        {
            //if not in group
            QModelIndex baseId = index_map.key(mainSelectionIndex());
            int newId = index_map.value(index(baseId.row()-1, baseId.column(),baseId.parent())) ;
            m_entries.move(mainSelectionIndex(), newId);
        }
        redraw_model();
        setSelectedItem(m_entries.indexOf(selected_item));
    }
}
void Scene::moveRowDown()
{
    Scene_item* selected_item = item(mainSelectionIndex());
    if(index_map.key(mainSelectionIndex()).row() < rowCount(index_map.key(mainSelectionIndex()).parent())-1)
    {
        if(item(mainSelectionIndex())->has_group >0)
        {
            Scene_group_item* group = selected_item->parentGroup();
            if(group)
            {
                int id = group->getChildren().indexOf(selected_item);
                group->moveDown(id);
            }
        }
        else
        {
            //if not in group
            QModelIndex baseId = index_map.key(mainSelectionIndex());
            int newId = index_map.value(index(baseId.row()+1, baseId.column(),baseId.parent())) ;
            m_entries.move(mainSelectionIndex(), newId);
        }
        redraw_model();
        setSelectedItem(m_entries.indexOf(selected_item));
    }
}
Scene::Item_id Scene::mainSelectionIndex() const {
    return selected_item;
}

QList<int> Scene::selectionIndices() const {
    return selected_items_list;
}

int Scene::selectionAindex() const {
    return item_A;
}

int Scene::selectionBindex() const {
    return item_B;
}

QItemSelection Scene::createSelection(int i)
{

    return QItemSelection(index_map.keys(i).at(0),
                          index_map.keys(i).at(4));
}

QItemSelection Scene::createSelectionAll()
{
    return QItemSelection(index_map.keys(0).at(0),
                          index_map.keys(m_entries.size() - 1).at(4));
}

void Scene::itemChanged()
{
    CGAL::Three::Scene_item* item = qobject_cast<CGAL::Three::Scene_item*>(sender());
    if(item)
        itemChanged(item);
}

void Scene::itemChanged(Item_id i)
{
    if(i < 0 || i >= m_entries.size())
        return;

  Q_EMIT dataChanged(this->createIndex(i, 0),
                     this->createIndex(i, LastColumn));
  //  Q_EMIT restoreCollapsedState();
}

void Scene::itemChanged(CGAL::Three::Scene_item* /* item */)
{
  Q_EMIT dataChanged(this->createIndex(0, 0),
                     this->createIndex(m_entries.size() - 1, LastColumn));
}

bool SceneDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
                                const QStyleOptionViewItem &option,
                                const QModelIndex &index)
{
    QAbstractProxyModel* proxyModel = dynamic_cast<QAbstractProxyModel*>(model);
    Q_ASSERT(proxyModel);
    Scene *scene = dynamic_cast<Scene*>(proxyModel->sourceModel());
    Q_ASSERT(scene);
    int id = scene->index_map[proxyModel->mapToSource(index)];
    switch(index.column()) {
    case Scene::VisibleColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                int x = mouseEvent->pos().x() - option.rect.x();
                if(x >= (option.rect.width() - size)/2 &&
                        x <= (option.rect.width() + size)/2) {
                    model->setData(index, !model->data(index).toBool());
                }
            }
            return false; //so that the selection can change
        }
        return true;
        break;
    case Scene::ColorColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                QColor color =
                        QColorDialog::getColor(model->data(index).value<QColor>(),
                                               0/*,
                                               tr("Select color"),
                                               QColorDialog::ShowAlphaChannel*/);
                if (color.isValid()) {
                    model->setData(index, color );
                }
            }
        }
        else if(event->type() == QEvent::MouseButtonDblClick) {
            return true; // block double-click
        }
        return false;
        break;
    case Scene::RenderingModeColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                // Switch rendering mode
                /*RenderingMode*/int rendering_mode = model->data(index, ::Qt::EditRole).toInt();
                rendering_mode = (rendering_mode+1) % NumberOfRenderingMode;
                model->setData(index, rendering_mode);
            }
        }
        else if(event->type() == QEvent::MouseButtonDblClick) {
            return true; // block double-click
        }
        return false;
        break;
    case Scene::ABColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            if(id == scene->item_B) {
                scene->item_A = id;
                scene->item_B = -1;
            }
            else if(id == scene->item_A) {
                scene->item_B = id;
                scene->item_A = -1;
            }
            else if(scene->item_A == -1) {
                scene->item_A = id;
            }
            else {
                scene->item_B = id;
            }
            scene->dataChanged(scene->createIndex(0, Scene::ABColumn),
                               scene->createIndex(scene->rowCount() - 1, Scene::ABColumn));
        }
        return false;
        break;
    default:
        return QItemDelegate::editorEvent(event, model, option, index);
    }
}

void SceneDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const
{
    QModelIndex test = proxy->mapToSource(index);
    if (index.column() != Scene::VisibleColumn) {
        QItemDelegate::paint(painter, option, index);
    } else {
        const QAbstractItemModel *model = index.model();

        QPalette::ColorGroup cg = (option.state & QStyle::State_Enabled) ?
                    (option.state & QStyle::State_Active) ? QPalette::Normal : QPalette::Inactive : QPalette::Disabled;

        if (option.state & QStyle::State_Selected)
            painter->fillRect(option.rect, option.palette.color(cg, QPalette::Highlight));
        bool checked = model->data(index, ::Qt::DisplayRole).toBool();
        int width = option.rect.width();
        int height = option.rect.height();
        size = (std::min)(width, height);
        int x = option.rect.x() + (option.rect.width() / 2) - (size / 2);;
        int y = option.rect.y() + (option.rect.height() / 2) - (size / 2);
        if(test.row()>=0 && test.row()<scene->m_entries.size()){

            if(checked) {
                painter->drawPixmap(x, y, checkOnPixmap.scaled(QSize(size, size),
                                                               ::Qt::KeepAspectRatio,
                                                               ::Qt::SmoothTransformation));
            }
            else {
                painter->drawPixmap(x, y, checkOffPixmap.scaled(QSize(size, size),
                                                                ::Qt::KeepAspectRatio,
                                                                ::Qt::SmoothTransformation));
            }
        }
        drawFocus(painter, option, option.rect); // since we draw the grid ourselves
    }
}

void Scene::setItemVisible(int index, bool b)
{
    if( index < 0 || index >= m_entries.size() )
        return;
    m_entries[index]->setVisible(b);
  Q_EMIT dataChanged(this->createIndex(index, VisibleColumn),
                     this->createIndex(index, VisibleColumn));
}

void Scene::setSelectionRay(double orig_x,
                            double orig_y,
                            double orig_z,
                            double dir_x,
                            double dir_y,
                            double dir_z)
{
    CGAL::Three::Scene_item* item = this->item(selected_item);
    if(item) item->select(orig_x,
                          orig_y,
                          orig_z,
                          dir_x,
                          dir_y,
                          dir_z);
}

void Scene::setItemA(int i)
{
    item_A = i;
    if(item_A == item_B)
    {
        item_B = -1;
    }
  Q_EMIT dataChanged(this->createIndex(0, ABColumn),
                     this->createIndex(m_entries.size()-1, ABColumn));
}

void Scene::setItemB(int i)
{
    item_B = i;
    if(item_A == item_B)
    {
        item_A = -1;
    }
  Q_EMIT updated();
  Q_EMIT dataChanged(this->createIndex(0, ABColumn),
                     this->createIndex(m_entries.size()-1, ABColumn));
}

Scene::Bbox Scene::bbox() const
{
    if(m_entries.empty())
        return Bbox(0,0,0,0,0,0);

    bool bbox_initialized = false;
    Bbox bbox = Bbox(0,0,0,0,0,0);
    Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
    {
        if(item->isFinite() && !item->isEmpty() && item->visible()) {
            if(bbox_initialized) {

                bbox = bbox + item->bbox();
            }
            else {
                bbox = item->bbox();
                bbox_initialized = true;

            }
        }

    }
    return bbox;
}

QList<Scene_item*> Scene::item_entries() const
{
    return m_entries;
}
void Scene::redraw_model()
{
    //makes the hierarchy in the tree
    //clears the model
    clear();
    index_map.clear();
    //fills the model
    Q_FOREACH(Scene_item* item, m_entries)
    {
        organize_items(item, invisibleRootItem(), 0);
    }
    Q_EMIT restoreCollapsedState();
}
void Scene::changeGroup(Scene_item *item, CGAL::Three::Scene_group_item *target_group)
{
    //remove item from the containing group if any
    if(item->parentGroup())
    {
      if(item->parentGroup()->isChildLocked(item))
        return;
      item->parentGroup()->removeChild(item);
    }
    //add the item to the target group
    target_group->addChild(item);
    item->moveToGroup(target_group);
    redraw_model();
    Q_EMIT updated();
}

float Scene::get_bbox_length() const
{
    return bbox().ymax()-bbox().ymin();
}

void Scene::printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface* viewer)
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printPrimitiveId if the item is a Scene_print_interface_item
    Scene_print_interface_item* item= dynamic_cast<Scene_print_interface_item*>(it);
    if(item)
      item->printPrimitiveId(point, viewer);
  }
}
void Scene::printPrimitiveIds(CGAL::Three::Viewer_interface* viewer)
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printPrimitiveIds if the item is a Scene_print_interface_item
    Scene_print_interface_item* item= dynamic_cast<Scene_print_interface_item*>(it);
    if(item)
      item->printPrimitiveIds(viewer);
  }
}
void Scene::updatePrimitiveIds(CGAL::Three::Viewer_interface* viewer, CGAL::Three::Scene_item* it)
{
  if(it)
  {
    //Only call printPrimitiveIds if the item is a Scene_print_interface_item
    Scene_print_interface_item* item= dynamic_cast<Scene_print_interface_item*>(it);
    if(item)
    {
      //As this function works as a toggle, the first call hides the ids and the second one  shows them again,
      //thereby triggering their re-computation.
      item->printPrimitiveIds(viewer);
      item->printPrimitiveIds(viewer);
    }
  }
}
bool Scene::testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer)
{
    CGAL::Three::Scene_item *i = item(mainSelectionIndex());
    if(i && i->visible())
    {
        bool res = i->testDisplayId(x,y,z, viewer);
        return res;
    }
    else
      return false;
}
#include "Scene_find_items.h"

void Scene::organize_items(Scene_item* item, QStandardItem* root, int loop)
{
    if(item->has_group <= loop)
    {
        QList<QStandardItem*> list;
        for(int i=0; i<5; i++)
        {
            list<<new QStandardItem();
            list.at(i)->setEditable(false);

        }
        root->appendRow(list);
        for(int i=0; i<5; i++){
            index_map[list.at(i)->index()] = m_entries.indexOf(item);
        }
        CGAL::Three::Scene_group_item* group =
                qobject_cast<CGAL::Three::Scene_group_item*>(item);
        if(group)
        {
            Q_FOREACH(Scene_item*child, group->getChildren())
            {
                organize_items(child, list.first(), loop+1);
            }
        }
    }
}

void Scene::setExpanded(QModelIndex id)
{
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item(getIdFromModelIndex(id)));
    if(group)
    {
        group->setExpanded(true);
    }
}
void Scene::setCollapsed(QModelIndex id)
{
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item(getIdFromModelIndex(id)));
    if(group)
    {
        group->setExpanded(false);
    }
}

int Scene::getIdFromModelIndex(QModelIndex modelId)const
{
    return index_map.value(modelId);
}

QList<QModelIndex> Scene::getModelIndexFromId(int id) const
{
    return index_map.keys(id);
}

void Scene::addGroup(Scene_group_item* group)
{
    connect(this, SIGNAL(drawFinished()), group, SLOT(resetDraw()));
    group->setScene(this);
}

namespace scene { namespace details {

Q_DECL_EXPORT
CGAL::Three::Scene_item*
findItem(const CGAL::Three::Scene_interface* scene_interface,
         const QMetaObject& metaobj,
         QString name, Scene_item_name_fn_ptr fn) {
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    if(!scene) return 0;
    Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
       CGAL::Three::Scene_item* ptr = qobject_cast<CGAL::Three::Scene_item*>(metaobj.cast(item));
        if(ptr && ((ptr->*fn)() == name)) return ptr;
    }
    return 0;
}

Q_DECL_EXPORT
QList<CGAL::Three::Scene_item*>
findItems(const CGAL::Three::Scene_interface* scene_interface,

          const QMetaObject&,
          QString name, Scene_item_name_fn_ptr fn)
{
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    QList<CGAL::Three::Scene_item*> list;
    if(!scene) return list;

    Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
        CGAL::Three::Scene_item* ptr = qobject_cast<CGAL::Three::Scene_item*>(item);
        if(ptr && ((ptr->*fn)() == name)) {
            list << ptr;
        }
    }
    return list;
}

} // end namespace details
                } // end namespace scene
