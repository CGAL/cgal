#include "GlSplat/GlSplat.h"



#include "config.h"
#include "Scene.h"

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

    if(ms_splatting==0)
        ms_splatting  = new GlSplat::SplatRenderer();
    ms_splattingCounter++;
    viewItem = invisibleRootItem();

    m_group_entries.append(new Scene_group_item("new group"));

}
Scene::Item_id
Scene::addItem(Scene_item* item)
{
    Bbox bbox_before = bbox();
    m_entries.push_back(item);
    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    connect(item, SIGNAL(renderingModeChanged()),
            this, SLOT(callDraw()));
    if(bbox_before + item->bbox() != bbox_before)
{ Q_EMIT updated_bbox(); }
    QStandardItemModel::beginResetModel();
    QList<QStandardItem*> list;
    for(int i=0; i<5; i++)
    {
        list<<new QStandardItem();
        list.at(i)->setEditable(false);
    }
    viewItem->appendRow(list);
    for(int i=0; i<5; i++){
        index_map[list.at(i)->index()] = m_entries.size() -1;
    }

    Q_EMIT updated();
    QStandardItemModel::endResetModel();
    Item_id id = m_entries.size() - 1;
  Q_EMIT newItem(id);
    return id;
}

Scene_item*
Scene::replaceItem(Scene::Item_id index, Scene_item* item, bool emit_item_about_to_be_destroyed)
{
    if(index < 0 || index >= m_entries.size())
        return 0;

    if(emit_item_about_to_be_destroyed) {
    Q_EMIT itemAboutToBeDestroyed(m_entries[index]);
    }

    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    std::swap(m_entries[index], item);
    if ( item->isFinite() && !item->isEmpty() &&
         m_entries[index]->isFinite() && !m_entries[index]->isEmpty() &&
         item->bbox()!=m_entries[index]->bbox() )
    {
    Q_EMIT updated_bbox();
    }
  Q_EMIT updated();
    itemChanged(index);
    // QStandardItemModel::reset();
    return item;
}

int
Scene::erase(int index)
{
    if(index < 0 || index >= m_entries.size())
        return -1;

    Scene_item* item = m_entries[index];
    Q_FOREACH(Scene_group_item* group, m_group_entries)
    {
        if(group->getChildren().contains(item))
            group->removeChild(item);
        if (group->getChildren().isEmpty())
        {
            m_group_entries.removeOne(group);
            m_entries.removeOne(group);
            check_first_group();
        }
    }

  Q_EMIT itemAboutToBeDestroyed(item);
    delete item;
    m_entries.removeAt(index);
  selected_item = -1;

  QStandardItemModel::beginResetModel();
  Q_EMIT updated();
  QStandardItemModel::endResetModel();
    if(--index >= 0)
        return index;
    if(!m_entries.isEmpty())
        return 0;
    return -1;

}
void Scene::check_empty_group(Scene_item* item)
{
    Q_FOREACH(Scene_group_item* group, m_group_entries)
    {
        if(group->getChildren().contains(item))
        {
            group->removeChild(item);
            if (group->getChildren().isEmpty())
            {
                check_empty_group(group);
                m_group_entries.removeOne(group);
                m_entries.removeOne(group);
            }
        }
    }
       check_first_group();
}

int
Scene::erase(QList<int> indices)
{
    clear();
    index_map.clear();
    QList<Scene_item*> to_be_removed;
    viewItem = invisibleRootItem();
    int max_index = -1;
    Q_FOREACH(int index, indices) {
        if(index < 0 || index >= m_entries.size())
            continue;

        max_index = (std::max)(max_index, index);

        Scene_item* item = m_entries[index];
        to_be_removed.push_back(item);
    }

    Q_FOREACH(Scene_item* item, to_be_removed) {
      Scene_group_item* group =
              qobject_cast<Scene_group_item*>(item);
    if(group)
    {
        QList<int> to_erase;
        Q_FOREACH(Scene_item* child, group->getChildren())
            to_erase<<m_entries.indexOf(child);
        erase(to_erase);
    }
    check_empty_group(item);

    Q_EMIT itemAboutToBeDestroyed(item);
    delete item;
    m_entries.removeAll(item);
  }
  selected_item = -1;
  Q_FOREACH(Scene_item* item, m_entries)
  {
    organize_items(item, viewItem, 0);
  }
  QStandardItemModel::beginResetModel();
  Q_EMIT updated();
  QStandardItemModel::endResetModel();

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

Scene::~Scene()
{
    Q_FOREACH(Scene_item* item_ptr, m_entries)
    {
        delete item_ptr;
    }
    m_entries.clear();

    if((--ms_splattingCounter)==0)
        delete ms_splatting;
}

Scene_item*
Scene::item(Item_id index) const
{
    return m_entries.value(index); // QList::value checks bounds
}

Scene::Item_id
Scene::item_id(Scene_item* scene_item) const
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

    const Scene_item* item = m_entries[index];
    Scene_item* new_item = item->clone();
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

}

bool
Scene::keyPressEvent(QKeyEvent* e){
    bool res=false;
    for (QList<int>::iterator it=selected_items_list.begin(),endit=selected_items_list.end();
         it!=endit;++it)
    {
        Scene_item* item=m_entries[*it];
        res |= item->keyPressEvent(e);
    }
    return res;
}

void
Scene::draw()
{
    draw_aux(false, 0);
}
void
Scene::draw(Viewer_interface* viewer)
{
    draw_aux(false, viewer);
}
void
Scene::drawWithNames()
{
    draw_aux(true, 0);
}
void
Scene::drawWithNames(Viewer_interface* viewer)
{
    draw_aux(true, viewer);
}

void
Scene::draw_aux(bool with_names, Viewer_interface* viewer)
{
    if(!ms_splatting->viewer_is_set)
        ms_splatting->setViewer(viewer);
    // Flat/Gouraud OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        if(with_names) {
            viewer->glPushName(index);
        }
        Scene_item& item = *m_entries[index];
        if(item.visible())
        {
            if(item.renderingMode() == Flat || item.renderingMode() == FlatPlusEdges || item.renderingMode() == Gouraud)
            {
                viewer->glEnable(GL_LIGHTING);
                viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
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

                if(item.renderingMode() == Gouraud)
                    viewer->glShadeModel(GL_SMOOTH);
                else
                    viewer->glShadeModel(GL_FLAT);

                item.contextual_changed();
                if(viewer)
                    item.draw(viewer);
                else
                    item.draw();
            }
        }
        if(with_names) {
            viewer->glPopName();
        }
    }
glDepthFunc(GL_LEQUAL);
    // Wireframe OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        if(with_names) {
            viewer->glPushName(index);
        }
        Scene_item& item = *m_entries[index];
        if(item.visible())
        {
            if(item.renderingMode() == FlatPlusEdges || item.renderingMode() == Wireframe)
            {
                viewer->glDisable(GL_LIGHTING);
                viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
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
                    item.draw_edges(viewer);
                else
                    item.draw_edges();
            }
            else{
                if( item.renderingMode() == PointsPlusNormals ){
                    viewer->glDisable(GL_LIGHTING);
                    viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
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
                        item.draw_edges(viewer);
                    else
                        item.draw_edges();
                }
            }
        }
         item.contextual_changed();
        if(with_names) {
            viewer->glPopName();
        }
    }


    // Points OpenGL drawing
    for(int index = 0; index < m_entries.size(); ++index)
    {
        if(with_names) {
            viewer->glPushName(index);
        }
        Scene_item& item = *m_entries[index];
        if(item.visible())
        {
            if(item.renderingMode() == Points  || item.renderingMode() == PointsPlusNormals)
            {
                viewer->glDisable(GL_LIGHTING);
                viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
                viewer->glPointSize(2.f);
                viewer->glLineWidth(1.0f);

                if(viewer)
                    item.draw_points(viewer);
                else
                    item.draw_points();
            }
        }
         item.contextual_changed();
        if(with_names) {
            viewer->glPopName();
        }
    }
    glDepthFunc(GL_LESS);
    // Splatting
    if(!with_names && ms_splatting->isSupported())
    {

        ms_splatting->beginVisibilityPass();
        for(int index = 0; index < m_entries.size(); ++index)
        {
            Scene_item& item = *m_entries[index];
            if(item.visible() && item.renderingMode() == Splatting)
            {

                if(viewer)
                {
                    item.draw_splats(viewer);
                }
                else
                    item.draw_splats();
            }

        }
       ms_splatting->beginAttributePass();
         for(int index = 0; index < m_entries.size(); ++index)
        {  Scene_item& item = *m_entries[index];
            if(item.visible() && item.renderingMode() == Splatting)
            {
                viewer->glColor4d(item.color().redF(), item.color().greenF(), item.color().blueF(), item.color().alphaF());
                if(viewer)
                    item.draw_splats(viewer);
                else
                    item.draw_splats();
            }
        }
        ms_splatting->finalize();

    }
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

    int id = index_map[index];\
    if(id < 0 || id >= m_entries.size())
        return QVariant();
    if(role == ::Qt::ToolTipRole)
    {
        return m_entries[id]->toolTip();
    }
    switch(index.column())
    {
    case ColorColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(id)->color();
        else if(role == ::Qt::DecorationRole)
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
            return ::Qt::AlignCenter;
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
                return tr("Color");
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

    Scene_item* item = m_entries[id];
    if(!item) return false;
    switch(index.column())
    {
    case NameColumn:
        item->setName(value.toString());
        check_first_group();
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

bool Scene::dropMimeData(const QMimeData *data,
                         Qt::DropAction action,
                         int row,
                         int column,
                         const QModelIndex &parent)
{
    //gets the moving items
    QList<Scene_item*> items;
    Q_FOREACH(int i, selected_items_list)
        items << item(i);
    //Gets the group at the drop position
    Scene_group_item* group =
            qobject_cast<Scene_group_item*>(this->item(index_map[parent]));
    //if the drop item is not a group_item, then the drop action must be ignored
    if(!group)
    {
        //unless the drop zone is empty, which means the item should be removed from all groups.
        if(!parent.isValid())
        {
          Q_FOREACH(Scene_item* item, items)
            while(item->has_group!=0)
            {
              Q_FOREACH(Scene_group_item* group_item, m_group_entries)
                if(group_item->getChildren().contains(item))
                {
                  group_item->removeChild(item);
                  break;
                }
            }
        group_added();
        return true;
        }
        return false;
    }
      Q_FOREACH(Scene_item* item, items)
    changeGroup(item, group);
    //group->addChild(item(mainSelectionIndex()));
    group_added();
    return true;


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
    Scene_item* item = qobject_cast<Scene_item*>(sender());
    if(item)
        itemChanged(item);
}

void Scene::itemChanged(Item_id i)
{
    if(i < 0 || i >= m_entries.size())
        return;

  Q_EMIT dataChanged(this->createIndex(i, 0),
                     this->createIndex(i, LastColumn));
}

void Scene::itemChanged(Scene_item* /* item */)
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
    Scene_item* item = this->item(selected_item);
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
        return Bbox();

    bool bbox_initialized = false;
    Bbox bbox;
    Q_FOREACH(Scene_item* item, m_entries)
    {
        if(item->isFinite() && !item->isEmpty()) {
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
QList<Scene_group_item*> Scene::group_entries() const
{
    return m_group_entries;
}

QList<Scene_item*> Scene::item_entries() const
{
    return m_entries;
}

void Scene::check_first_group()
{


    if(m_group_entries.isEmpty() || m_group_entries.first()->name() != "new group")
    {
      m_group_entries.prepend(new Scene_group_item("new group"));
    }
}
void Scene::group_added()
{
    Q_FOREACH(Scene_group_item *group, m_group_entries)
    {
        if (group->getChildren().isEmpty())
        {
            m_group_entries.removeOne(group);
            m_entries.removeOne(group);
        }
    }
    check_first_group();

//makes the hierarchy in the tree
    //clears the model
    clear();
    index_map.clear();
    viewItem = invisibleRootItem();
    //fills the model
    Q_FOREACH(Scene_item* item, m_entries)
    {
    organize_items(item, viewItem, 0);
    }
}
//TO DO gerer les has_group
void Scene::changeGroup(Scene_item *item, Scene_group_item *target_group)
{
    //remove item from the containing group if any
 if(item->has_group!=0)
  Q_FOREACH(Scene_group_item* group, m_group_entries)
  {
    if(group->getChildren().contains(item))
    {
      group->removeChild(item);
      break;
    }
  }
 //add the item to the target group
 target_group->addChild(item);
 item->has_group = target_group->has_group +1;
}
void Scene::setGroupName(QString name)
{
 m_group_entries[0]->setName(name);
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
        Scene_group_item* group =
                qobject_cast<Scene_group_item*>(item);
        if(group)
        {
            Q_FOREACH(Scene_item*child, group->getChildren())
            {
                organize_items(child, list.first(), loop+1);
            }
        }
    }
}
namespace scene { namespace details {

Q_DECL_EXPORT
Scene_item*
findItem(const Scene_interface* scene_interface,
         const QMetaObject& metaobj,
         QString name, Scene_item_name_fn_ptr fn) {
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    if(!scene) return 0;
    Q_FOREACH(Scene_item* item, scene->entries()) {
        Scene_item* ptr = qobject_cast<Scene_item*>(metaobj.cast(item));
        if(ptr && ((ptr->*fn)() == name)) return ptr;
    }
    return 0;
}

Q_DECL_EXPORT
QList<Scene_item*>
findItems(const Scene_interface* scene_interface,
          const QMetaObject&,
          QString name, Scene_item_name_fn_ptr fn)
{
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    QList<Scene_item*> list;
    if(!scene) return list;

    Q_FOREACH(Scene_item* item, scene->entries()) {
        Scene_item* ptr = qobject_cast<Scene_item*>(item);
        if(ptr && ((ptr->*fn)() == name)) {
            list << ptr;
        }
    }
    return list;
}

} // end namespace details
                } // end namespace scene
