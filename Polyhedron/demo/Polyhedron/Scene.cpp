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

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
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
}

Scene::Item_id
Scene::addItem(Scene_item* item)
{
  m_entries.push_back(item);

  connect(item, SIGNAL(itemChanged()),
          this, SLOT(itemChanged()));
  emit updated_bbox();
  emit updated();
  QAbstractListModel::reset();
  Item_id id = m_entries.size() - 1;
  emit newItem(id);
  return id;
}

Scene_item*
Scene::replaceItem(Scene::Item_id index, Scene_item* item)
{
  if(index < 0 || index >= m_entries.size())
    return 0;

  connect(item, SIGNAL(itemChanged()),
          this, SLOT(itemChanged()));
  std::swap(m_entries[index], item);

  emit updated_bbox();
  emit updated();
  itemChanged(index);
  QAbstractListModel::reset();
  return item;
}

int
Scene::erase(int index)
{
  if(index < 0 || index >= m_entries.size())
    return -1;

  Scene_item* item = m_entries[index];
  emit itemAboutToBeDestroyed(item);
  delete item;
  m_entries.removeAt(index);

  selected_item = -1;
  emit updated();
  QAbstractListModel::reset();

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
    emit itemAboutToBeDestroyed(item);
    delete item;
  }

  Q_FOREACH(Scene_item* item, to_be_removed) {
    m_entries.removeAll(item);
  }

  selected_item = -1;
  emit updated();
  QAbstractListModel::reset();

  int index = max_index + 1 - indices.size();
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
}

Scene_item*
Scene::item(Item_id index) const
{
  return m_entries.value(index); // QList::value checks bounds
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
}

// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define lighter light
#endif

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
  draw_aux(false);
}
void 
Scene::drawWithNames()
{
  draw_aux(true);
}

void 
Scene::draw_aux(bool with_names)
{
  // Flat/Gouraud OpenGL drawing
  for(int index = 0; index < m_entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *m_entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == Flat || item.renderingMode() == FlatPlusEdges || item.renderingMode() == Gouraud)
      {
	::glEnable(GL_LIGHTING);
	::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        ::glPointSize(2.f);
        ::glLineWidth(1.0f);
	if(index == selected_item)
	  CGALglcolor(item.color().lighter(120));
	else
	  CGALglcolor(item.color());
	if(item.renderingMode() == Gouraud)
	  ::glShadeModel(GL_SMOOTH);
	else
	  ::glShadeModel(GL_FLAT);

        item.draw();
      }
    }
    if(with_names) {
      ::glPopName();
    }
  }

  // Wireframe OpenGL drawing
  for(int index = 0; index < m_entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *m_entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == FlatPlusEdges || item.renderingMode() == Wireframe)
      {
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        ::glPointSize(2.f);
        ::glLineWidth(1.0f);
        if(index == selected_item)
          CGALglcolor(Qt::black);
        else
          CGALglcolor(item.color().lighter(50));
        
        item.draw_edges();
      }
      else{
        if( item.renderingMode() == PointsPlusNormals ){
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        ::glPointSize(2.f);
        ::glLineWidth(1.0f);
        if(index == selected_item)
          CGALglcolor(item.color().lighter(120));
        else
          CGALglcolor(item.color());
        item.draw_edges();
        }
      }
      if(with_names) {
        ::glPopName();
      }
    }
  }

  // Points OpenGL drawing
  for(int index = 0; index < m_entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *m_entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == Points  || item.renderingMode() == PointsPlusNormals)
      {
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
        ::glPointSize(2.f);
        ::glLineWidth(1.0f);
        if(index == selected_item)
          CGALglcolor(Qt::black);
        else
          CGALglcolor(item.color().lighter(50));
        
        item.draw_points();
      }
      if(with_names) {
        ::glPopName();
      }
    }
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
    item->changed();
    emit dataChanged(index, index);
    return true;
    break;
  case ColorColumn:
    item->setColor(value.value<QColor>());
    item->changed();
    emit dataChanged(index, index);
    return true;
    break;
  case RenderingModeColumn:
  {
    RenderingMode rendering_mode = static_cast<RenderingMode>(value.toInt());
    // Find next supported rendering mode
    while ( ! item->supportsRenderingMode(rendering_mode) ) {
      rendering_mode = static_cast<RenderingMode>( (rendering_mode+1) % NumberOfRenderingMode );
    }
    item->setRenderingMode(rendering_mode);
    item->changed();
    emit dataChanged(index, index);
    return true;
    break;
  }
  case VisibleColumn:
    item->setVisible(value.toBool());
    item->changed();
    emit dataChanged(index, index);
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
  return QItemSelection(QAbstractItemModel::createIndex(i, 0),
    QAbstractItemModel::createIndex(i, LastColumn));
}

QItemSelection Scene::createSelectionAll()
{
  return QItemSelection(QAbstractItemModel::createIndex(0, 0),
    QAbstractItemModel::createIndex(m_entries.size() - 1 , LastColumn));
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

  m_entries[i]->changed();
  emit dataChanged(QAbstractItemModel::createIndex(i, 0),
    QAbstractItemModel::createIndex(i, LastColumn));
}

void Scene::itemChanged(Scene_item* item)
{
  item->changed();
  emit dataChanged(QAbstractItemModel::createIndex(0, 0),
    QAbstractItemModel::createIndex(m_entries.size() - 1, LastColumn));
}

bool SceneDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
				const QStyleOptionViewItem &option,
				const QModelIndex &index)
{
  Scene *scene = static_cast<Scene*>(model);
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
      scene->dataChanged(scene->createIndex(Scene::ABColumn, 0),
	scene->createIndex(Scene::ABColumn, scene->rowCount()));
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
  emit dataChanged(QAbstractItemModel::createIndex(index, VisibleColumn),
    QAbstractItemModel::createIndex(index, VisibleColumn));
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
  emit dataChanged(QAbstractItemModel::createIndex(0, ABColumn),
    QAbstractItemModel::createIndex(m_entries.size()-1, ABColumn));
}

void Scene::setItemB(int i)
{
  item_B = i;
  if(item_A == item_B)
  {
    item_A = -1;
  }
  emit updated();
  emit dataChanged(QAbstractItemModel::createIndex(0, ABColumn),
    QAbstractItemModel::createIndex(m_entries.size()-1, ABColumn));
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

#include "Scene.moc"
