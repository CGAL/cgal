#include "config.h"
#include "Scene.h"
#include "Scene_item.h"

#include <QString>
#include <QGLWidget>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QColorDialog>
#include <QApplication>
#include <QPointer>

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

Scene::Scene(QObject* parent)
  : QAbstractListModel(parent),
    selected_item(-1),
    viewEdges(true)
{
}

Scene::Item_id
Scene::addItem(Scene_item* item)
{
  entries.push_back(item);

  emit updated_bbox();
  emit updated();
  QAbstractListModel::reset();
  return entries.size() - 1;
}

int
Scene::erase(int index)
{
  if(index < 0 || index >= entries.size())
    return -1;

  Scene_item* item = entries[index];
  emit itemAboutToBeDestroyed(item);
  delete item;
  entries.removeAt(index);

  selected_item = -1;
  emit updated();
  QAbstractListModel::reset();

  if(--index >= 0)
    return index;
  if(!entries.isEmpty())
    return 0;
  return -1;
}

Scene::~Scene()
{
  Q_FOREACH(Scene_item* item_ptr, entries)
  {
    delete item_ptr;
  }
  entries.clear();
}

Scene_item*
Scene::item(Item_id index) const
{
  return entries.value(index); // QList::value checks bounds
}

size_t
Scene::numberOfEntries() const
{
  return entries.size();
}

Scene::Item_id
Scene::duplicate(Item_id index)
{
  if(index < 0 || index >= entries.size())
    return -1;

  const Scene_item* item = entries[index];
  Scene_item* new_item = item->clone();
  if(new_item) {
    new_item->setName(tr("%1 (copy)").arg(item->name()));
    new_item->setColor(item->color());
    new_item->setVisible(item->visible());
    addItem(new_item);
    return entries.size() - 1;
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
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == Fill)
      {
	::glEnable(GL_LIGHTING);
	::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	if(index == selected_item)
	  CGALglcolor(item.color().lighter(120));
	else
	  CGALglcolor(item.color());

        item.draw();
      }
    }
    if(with_names) {
      ::glPopName();
    }
  }
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
    if(item.visible())
    {
      if(viewEdges || item.renderingMode() == Wireframe)
      {
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        if(index == selected_item)
          CGALglcolor(Qt::black);
        else
          CGALglcolor(item.color().lighter(50));
        
        item.draw_edges();
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
    return entries.size();
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

  if(index.row() < 0 || index.row() >= entries.size())
    return QVariant();

  if(role == ::Qt::ToolTipRole)
  {
    return entries[index.row()]->toolTip();
  }
  switch(index.column())
  {
  case ColorColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return entries.value(index.row())->color();
    else if(role == ::Qt::DecorationRole)
      return entries.value(index.row())->color();
    break;
  case NameColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return entries.value(index.row())->name();
    if(role == ::Qt::FontRole)
      return entries.value(index.row())->font();
    break;
  case RenderingModeColumn:
    if(role == ::Qt::DisplayRole) {
      if(entries.value(index.row())->renderingMode() == Wireframe)
	return tr("wire");
      else return tr("fill");
    }
    else if(role == ::Qt::EditRole) {
      return static_cast<int>(entries.value(index.row())->renderingMode());
    }
    else if(role == ::Qt::TextAlignmentRole) {
      return ::Qt::AlignCenter;
    }
    break;
  case VisibleColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return entries.value(index.row())->visible();
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
      case VisibleColumn:
	return tr("View");
	break;
      default:
	return QVariant();
      }
    }
    else if(role == ::Qt::ToolTipRole) {
      if(section == RenderingModeColumn) {
	return tr("Rendering mode (fill/fireframe)");
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

  if(index.row() < 0 || index.row() >= entries.size())
    return false;

  Scene_item* item = entries[index.row()];
  if(!item) return false;
  switch(index.column())
  {
  case NameColumn:
    item->setName(value.toString());
    emit dataChanged(index, index);
    return true;
    break;
  case ColorColumn:
    item->setColor(value.value<QColor>());
    emit dataChanged(index, index);
    return true;
    break;
  case RenderingModeColumn:
    item->setRenderingMode(static_cast<RenderingMode>(value.toInt()));
    emit dataChanged(index, index);
    return true;
    break;
  case VisibleColumn:
    item->setVisible(value.toBool());
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

QItemSelection Scene::createSelection(int i)
{
  return QItemSelection(QAbstractItemModel::createIndex(i, 0),
    QAbstractItemModel::createIndex(i, LastColumn));
}

void Scene::itemChanged(Item_id i)
{
  if(i < 0 || i >= entries.size())
    return;

  entries[i]->changed();
  emit dataChanged(QAbstractItemModel::createIndex(i, 0),
    QAbstractItemModel::createIndex(i, LastColumn));
}

void Scene::itemChanged(Scene_item* item)
{
  item->changed();
  emit dataChanged(QAbstractItemModel::createIndex(0, 0),
    QAbstractItemModel::createIndex(entries.size() - 1, LastColumn));
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
      QColor color = 
        QColorDialog::getColor(model->data(index).value<QColor>(),
                               0/*,
                               tr("Select color"),
                               QColorDialog::ShowAlphaChannel*/);
      if (color.isValid()) {
	model->setData(index, color );
      }
    }
    else if(event->type() == QEvent::MouseButtonDblClick) {
      return true; // block double-click
    }
    return false;
    break;
  case Scene::RenderingModeColumn:
    if (event->type() == QEvent::MouseButtonPress) {
      RenderingMode rendering_mode = 
	static_cast<RenderingMode>(model->data(index, ::Qt::EditRole).toInt());
      if(rendering_mode == Wireframe)
	model->setData(index, static_cast<int>(Fill));
      else 
	model->setData(index, static_cast<int>(Wireframe));
    }
    else if(event->type() == QEvent::MouseButtonDblClick) {
      return true; // block double-click
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
  if( index < 0 || index >= entries.size() )
    return;
  entries[index]->setVisible(b);
  emit dataChanged(QAbstractItemModel::createIndex(index, VisibleColumn),
    QAbstractItemModel::createIndex(index, VisibleColumn));
}

Scene::Bbox Scene::bbox() const
{
  if(entries.empty())
    return Bbox();

  bool bbox_initialized = false;
  Bbox bbox;
  Q_FOREACH(Scene_item* item, entries) 
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
