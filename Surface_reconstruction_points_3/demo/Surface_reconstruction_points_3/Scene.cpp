
#ifdef CGAL_GLEW_ENABLED
# include "GlSplat/GlSplat.h"
#endif

#include "config.h"
#include "Scene.h"
#include "Scene_item.h"
#include "Scene_polyhedron_item.h"
#include "Point_set_scene_item.h"

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

#ifdef CGAL_GLEW_ENABLED
GlSplat::SplatRenderer* Scene::ms_splatting = 0;
int Scene::ms_splattingCounter = 0;
GlSplat::SplatRenderer* Scene::splatting()
{
  assert(ms_splatting!=0 && "A Scene object must be created before requesting the splatting object");
  return ms_splatting;
}
#endif

Scene::Scene(QObject* parent)
  : QAbstractListModel(parent),
    selected_item(-1)
{
#ifdef CGAL_GLEW_ENABLED
  if(ms_splatting==0)
    ms_splatting  = new GlSplat::SplatRenderer();
  ms_splattingCounter++;
#endif
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

// Erases a scene item.
// Returns the index of the polyhedra just before the one that is erased,
//  or just after. Returns -1 if the list is empty.
Scene::Item_id Scene::erase(Item_id index)
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

#ifdef CGAL_GLEW_ENABLED
  if((--ms_splattingCounter)==0)
    delete ms_splatting;
#endif
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

// Duplicates a scene item.
// Returns the ID of the new item (-1 on error).
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
    Item_id new_index = addItem(new_item);
    return new_index;
  }
  else
    return -1;
}

// Converts a polyhedron to a point set.
// Returns the ID of the new item (-1 on error).
Scene::Item_id
Scene::convertToPointSet(Item_id index)
{
  // Check index
  if(index < 0 || index >= entries.size())
    return -1;

  // Check if scene item is a polyhedron
  Scene_item* item = entries[index];
  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(item);
  if(poly_item == NULL || poly_item->polyhedron() == NULL)
    return -1;

  // Converts polyhedron to a point set
  Point_set_scene_item* new_item = new Point_set_scene_item(*poly_item->polyhedron());
  if(new_item) {
    new_item->setName(tr("%1 (point set)").arg(item->name()));
    new_item->setColor(item->color());
    new_item->setVisible(item->visible());
    Item_id new_index = addItem(new_item);

    // Hide polyhedron
    poly_item->setVisible(false);
    itemChanged(index);

    return new_index;
  }
  else
    return -1;
}

// Delete selection in a scene item
void Scene::deleteSelection(Item_id index)
{
  // Check index
  if(index < 0 || index >= entries.size())
    return;

  Scene_item* item = entries[index];
  if (item->isSelectionEmpty())
    return;

  item->deleteSelection();
  itemChanged(index);
}

// Reset selection mark in a scene item
void Scene::resetSelection(Item_id index)
{
  if(index < 0 || index >= entries.size())
    return;

  Scene_item* item = entries[index];
  item->resetSelection();
  itemChanged(index);
}

void Scene::initializeGL()
{
#ifdef CGAL_GLEW_ENABLED
  ms_splatting->init();
#endif
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
  // Flat/Gouraud OpenGL drawing
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
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
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
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
      if(with_names) {
        ::glPopName();
      }
    }
  }

  // Points OpenGL drawing
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == Points || item.renderingMode() == PointsPlusNormals)
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

#ifdef CGAL_GLEW_ENABLED
  // Splatting
  if(ms_splatting->isSupported())
  {
    ms_splatting->beginVisibilityPass();
    for(int index = 0; index < entries.size(); ++index)
    {
      Scene_item& item = *entries[index];
      if(item.visible() && item.renderingMode() == Splatting)
      {
        item.draw_splats();
      }
    }
    ms_splatting->beginAttributePass();
    for(int index = 0; index < entries.size(); ++index)
    {
      Scene_item& item = *entries[index];
      if(item.visible() && item.renderingMode() == Splatting)
      {
        if(index == selected_item)
          CGALglcolor(item.color().lighter(120));
        else
          CGALglcolor(item.color());
        item.draw_splats();
      }
    }
    ms_splatting->finalize();
  }
#endif

  // Normals OpenGL drawing
  for(int index = 0; index < entries.size(); ++index)
  {
    if(with_names) {
      ::glPushName(index);
    }
    Scene_item& item = *entries[index];
    if(item.visible())
    {
      if(item.renderingMode() == PointsPlusNormals)
      {
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        ::glPointSize(2.f);
        ::glLineWidth(1.0f);
        if(index == selected_item)
          CGALglcolor(item.color().lighter(120));
        else
          CGALglcolor(item.color());

        item.draw_normals();
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
      return entries.value(index.row())->renderingModeName();
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
	return Scene_item::renderingModeNameList();
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
    while ( !item->supportsRenderingMode(rendering_mode)
#ifdef CGAL_GLEW_ENABLED
         || (rendering_mode==Splatting && !Scene::splatting()->isSupported())
#endif
    )
    {
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
//  Scene *scene = static_cast<Scene*>(model);
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
      // Switch rendering mode
      /*RenderingMode*/int rendering_mode = model->data(index, ::Qt::EditRole).toInt();
      rendering_mode = (rendering_mode+1) % NumberOfRenderingMode;
      model->setData(index, rendering_mode);
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
