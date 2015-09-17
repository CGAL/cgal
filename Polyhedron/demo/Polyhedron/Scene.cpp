#include "GlSplat/GlSplat.h"



#include "config.h"
#include "Scene.h"
#include "Scene_item.h"

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




GlSplat::SplatRenderer* Scene::ms_splatting = 0;
int Scene::ms_splattingCounter = 0;
GlSplat::SplatRenderer* Scene::splatting()
{
    assert(ms_splatting!=0 && "A Scene object must be created before requesting the splatting object");
    return ms_splatting;
}

Scene::Scene(QObject* parent)
    : QAbstractListModel(parent),
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
    QAbstractListModel::beginResetModel();
    Q_EMIT updated();
    QAbstractListModel::endResetModel();
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
    // QAbstractListModel::reset();
    return item;
}

int
Scene::erase(int index)
{
    if(index < 0 || index >= m_entries.size())
        return -1;

    Scene_item* item = m_entries[index];
  Q_EMIT itemAboutToBeDestroyed(item);
    delete item;
    m_entries.removeAt(index);

  selected_item = -1;

  QAbstractListModel::beginResetModel();
  Q_EMIT updated();
  QAbstractListModel::endResetModel();

    if(--index >= 0)
        return index;
    if(!m_entries.isEmpty())
        return 0;
    return -1;
}

int
Scene::erase(QList<int> indices)
{
    QList<Scene_item*> to_be_removed;

    int max_index = -1;
    Q_FOREACH(int index, indices) {
        if(index < 0 || index >= m_entries.size())
            continue;
        max_index = (std::max)(max_index, index);
        Scene_item* item = m_entries[index];
        to_be_removed.push_back(item);
    }



  Q_FOREACH(Scene_item* item, to_be_removed) {
    Q_EMIT itemAboutToBeDestroyed(item);
    delete item;
    m_entries.removeAll(item);
  }

  selected_item = -1;
  QAbstractListModel::beginResetModel();
  Q_EMIT updated();
  QAbstractListModel::endResetModel();
 
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
                if(index == selected_item)
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
                if(index == selected_item)
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
                    if(index == selected_item)
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

int 
Scene::rowCount(const QModelIndex & parent) const
{
    if (parent.isValid())
        return 0;
    else
        return m_entries.size();
}

int 
Scene::columnCount(const QModelIndex & parent) const
{
    if (parent.isValid())
        return 0;
    else
        return NumberOfColumns;
}

QVariant 
Scene::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if(index.row() < 0 || index.row() >= m_entries.size())
        return QVariant();

    if(role == ::Qt::ToolTipRole)
    {
        return m_entries[index.row()]->toolTip();
    }
    switch(index.column())
    {
    case ColorColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(index.row())->color();
        else if(role == ::Qt::DecorationRole)
            return m_entries.value(index.row())->color();
        break;
    case NameColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(index.row())->name();
        if(role == ::Qt::FontRole)
            return m_entries.value(index.row())->font();
        break;
    case RenderingModeColumn:
        if(role == ::Qt::DisplayRole) {
            return m_entries.value(index.row())->renderingModeName();
        }
        else if(role == ::Qt::EditRole) {
            return static_cast<int>(m_entries.value(index.row())->renderingMode());
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignCenter;
        }
        break;
    case ABColumn:
        if(role == ::Qt::DisplayRole) {
            if(index.row() == item_A)
                return "A";
            if(index.row() == item_B)
                return "B";
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignCenter;
        }
        break;
    case VisibleColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(index.row())->visible();
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
    return QAbstractListModel::headerData(section, orientation, role);
}

Qt::ItemFlags 
Scene::flags ( const QModelIndex & index ) const
{
    if (index.isValid() && index.column() == NameColumn) {
        return QAbstractListModel::flags(index) | ::Qt::ItemIsEditable;
    }
    else {
        return QAbstractListModel::flags(index);
    }
}

bool 
Scene::setData(const QModelIndex &index, 
               const QVariant &value,
               int role)
{
    if( role != ::Qt::EditRole || !index.isValid() )
        return false;

    if(index.row() < 0 || index.row() >= m_entries.size())
        return false;

    Scene_item* item = m_entries[index.row()];
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
    Q_EMIT dataChanged(index, index);
        return true;
        break;
    }
    case VisibleColumn:
        item->setVisible(value.toBool());
    Q_EMIT dataChanged(index, index);
        return true;
    default:
        return false;
    }
    return false;
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
    return QItemSelection(this->createIndex(i, 0),
                          this->createIndex(i, LastColumn));
}

QItemSelection Scene::createSelectionAll()
{
    return QItemSelection(this->createIndex(0, 0),
                          this->createIndex(m_entries.size() - 1 , LastColumn));
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
    switch(index.column()) {
    case Scene::VisibleColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                int x = mouseEvent->pos().x() - option.rect.x();
                if(x >= (option.rect.width() - size)/2 &&
                        x <= (option.rect.width() + size)/2) {
                    model->setData(index, ! model->data(index).toBool() );
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
            if(index.row() == scene->item_B) {
                scene->item_A = index.row();
                scene->item_B = -1;
            }
            else if(index.row() == scene->item_A) {
                scene->item_B = index.row();
                scene->item_A = -1;
            }
            else if(scene->item_A == -1) {
                scene->item_A = index.row();
            }
            else {
                scene->item_B = index.row();
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

#include "Scene_find_items.h"

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
