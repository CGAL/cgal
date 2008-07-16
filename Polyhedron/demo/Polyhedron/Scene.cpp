#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QGLWidget>
#include <QMessageBox>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QColorDialog>
#include <QApplication>

#include "Scene_rendering.h"

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

const QColor Scene::defaultColor = QColor(100, 100, 255);

Scene::Scene(QObject* parent)
  : QAbstractListModel(parent),
    selected_item(-1),
    item_A(-1),
    item_B(-1)
{
}

Scene::~Scene()
{
  for(Polyhedra::iterator 
        poly_it = polyhedra.begin(),
        poly_end = polyhedra.end();
      poly_it != poly_end; ++poly_it) {
    this->destroy(poly_it->polyhedron_ptr);
  }
  polyhedra.clear();
}

int
Scene::open(QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"...").arg(filename);

  QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

  QFileInfo fileinfo(filename);
  std::ifstream in(filename.toUtf8());

  if(!in || !fileinfo.isFile() || ! fileinfo.isReadable()) {
    QMessageBox::critical(qobject_cast<QWidget*>(QObject::parent()),
                          tr("Cannot open file"),
                          tr("File %1 is not a readable file.").arg(filename));
    QApplication::restoreOverrideCursor();
		cerr << QString("\n");
    return -1;
  }

  // allocate new polyhedron
  Polyhedron* poly = this->new_polyhedron();
  this->load_polyhedron(poly, in);
  if(!in)
  {
    QMessageBox::critical(qobject_cast<QWidget*>(QObject::parent()),
                          tr("Cannot read file"),
                          tr("File %1 is not a valid OFF file.").arg(filename));
    QApplication::restoreOverrideCursor();
		cerr << QString("\n");
		destroy(poly);
    return -1;
  }

  addPolyhedron(poly, fileinfo.baseName());
  QApplication::restoreOverrideCursor();

  return polyhedra.size() - 1;
}

bool Scene::save(int index,
		 QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Saving file \"%1\"...").arg(filename);

  Polyhedron_entry entry = polyhedra[index];
  Polyhedron* poly = entry.polyhedron_ptr;

  QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

  QFileInfo fileinfo(filename);
  std::ofstream out(filename.toUtf8());

  if(!out || !fileinfo.isFile() || ! fileinfo.isWritable())
  {
    QMessageBox::critical(qobject_cast<QWidget*>(QObject::parent()),
                          tr("Cannot open file"),
                          tr("File %1 is not a writable file.").arg(filename));
    QApplication::restoreOverrideCursor();
	  cerr << QString("\n");
    return false;
  }

  this->save_polyhedron(poly, out);
  cerr << QString("ok\n");

  QApplication::restoreOverrideCursor();

  return true;
}

void Scene::addPolyhedron(Polyhedron* p,
                          QString name,
                          QColor color,
                          bool activated,
                          RenderingMode mode)
{
  Polyhedron_entry entry;
  entry.polyhedron_ptr = p;
  entry.name = name;
  entry.color = color;
  entry.activated = activated;
  entry.rendering_mode = mode;
  polyhedra.push_back(entry);

  selected_item = -1;
  emit updated_bbox();
  emit updated();
  QAbstractListModel::reset();
}

int
Scene::erase(int polyhedron_index)
{
  this->destroy(polyhedra[polyhedron_index].polyhedron_ptr);
  polyhedra.removeAt(polyhedron_index--);

  selected_item = -1;
  emit updated();
  QAbstractListModel::reset();

  if(polyhedron_index >= 0)
    return polyhedron_index;
  if(!polyhedra.isEmpty())
    return 0;
  return -1;
}

int
Scene::duplicate(int polyhedron_index)
{
  const Polyhedron_entry& entry = polyhedra[polyhedron_index];
  Polyhedron* poly = this->copy_polyhedron(entry.polyhedron_ptr);

  addPolyhedron(poly,
                tr("%1 (copy)").arg(entry.name),
                entry.color,
                entry.activated);

  return polyhedra.size() - 1;
}

void 
Scene::draw()
{
  for(int index = 0; index < polyhedra.size(); ++index)
  {
    Polyhedron_entry& entry = polyhedra[index];
    if(entry.activated)
		{
      Polyhedron* poly = entry.polyhedron_ptr;
      if(entry.rendering_mode == Fill)
			{
        ::glEnable(GL_LIGHTING);
        if(index == selected_item)
          CGALglcolor(entry.color.lighter(120));
        else
          CGALglcolor(entry.color);
        
	gl_render_facets(poly);
      }

      if(index == selected_item)
        CGALglcolor(entry.color.lighter(70));
      else
        CGALglcolor(entry.color.lighter(50));

			// superimpose edges
      ::glDisable(GL_LIGHTING);
      gl_render_edges(poly);
    }
  }
}

int 
Scene::rowCount(const QModelIndex & parent) const
{
  if (parent.isValid())
    return 0;
  else
    return polyhedra.size();
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
  
  switch(index.column())
  {
  case ColorColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return polyhedra.value(index.row()).color;
    else if(role == ::Qt::DecorationRole)
      return polyhedra.value(index.row()).color;
    break;
  case NameColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return polyhedra.value(index.row()).name;
    break;
  case RenderingModeColumn:
    if(role == ::Qt::DisplayRole) {
      if(polyhedra.value(index.row()).rendering_mode == Scene::Wireframe)
        return tr("wire");
      else return tr("fill");
    }
    else if(role == ::Qt::EditRole) {
      return static_cast<int>(polyhedra.value(index.row()).rendering_mode);
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
  case ActivatedColumn:
    if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
      return polyhedra.value(index.row()).activated;
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
      case ActivatedColumn:
        return tr("Activated");
        break;
      default:
        return QVariant();
      }
    }
    else if(role == ::Qt::ToolTipRole) {
      if(section == RenderingModeColumn) {
        return tr("Rendering mode (fill/fireframe)");
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

  Polyhedron_entry& entry = polyhedra[index.row()];
  switch(index.column())
  {
  case NameColumn:
    entry.name = value.toString();
    emit dataChanged(index, index);
    return true;
    break;
  case ColorColumn:
    entry.color = value.value<QColor>();
    emit dataChanged(index, index);
    return true;
    break;
  case RenderingModeColumn:
    entry.rendering_mode = static_cast<RenderingMode>(value.toInt());
    emit dataChanged(index, index);
    return true;
    break;
  case ActivatedColumn:
    entry.activated = value.toBool();
    emit dataChanged(index, index);
    return true;
  default:
    return false;
  }
  return false;
}

Polyhedron* Scene::polyhedron(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return 0;
  else 
    return polyhedra[index].polyhedron_ptr;
}

QString Scene::polyhedronName(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return QString();
  else 
    return polyhedra[index].name;
}

QColor Scene::polyhedronColor(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return QColor();
  else 
    return polyhedra[index].color;
}

bool Scene::isPolyhedronActivated(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return false;
  else 
    return polyhedra[index].activated;
}

Scene::RenderingMode Scene::polyhedronRenderingMode(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return RenderingMode();
  else 
    return polyhedra[index].rendering_mode;
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

void Scene::polyhedronChanged(int i)
{
  emit dataChanged(QAbstractItemModel::createIndex(i, 0),
                   QAbstractItemModel::createIndex(i, LastColumn));
}

void Scene::polyhedronChanged(Polyhedron*)
{
  emit dataChanged(QAbstractItemModel::createIndex(0, 0),
                   QAbstractItemModel::createIndex(polyhedra.size() - 1, LastColumn));
}

bool SceneDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
                                const QStyleOptionViewItem &option,
                                const QModelIndex &index)
{
  Scene *scene = static_cast<Scene*>(model);
  Q_ASSERT(scene);
  switch(index.column()) {
  case Scene::ActivatedColumn:
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
      QColor color = QColorDialog::getColor(::Qt::green, 0);
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
      Scene::RenderingMode rendering_mode = 
        static_cast<Scene::RenderingMode>(model->data(index, ::Qt::EditRole).toInt());
      std::cerr << "render mode = " << rendering_mode << "\n";
      if(rendering_mode == Scene::Wireframe)
        model->setData(index, static_cast<int>(Scene::Fill));
      else 
        model->setData(index, static_cast<int>(Scene::Wireframe));
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
  if (index.column() != Scene::ActivatedColumn) {
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
