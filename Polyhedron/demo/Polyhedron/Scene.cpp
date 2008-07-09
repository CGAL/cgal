#include "Scene.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QGLWidget>
#include <QMessageBox>

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

Scene::Scene(QObject* parent)
  : QAbstractListModel(parent),
    selected_item(-1)
{
}

Scene::~Scene()
{
  for(Polyhedra::iterator 
        poly_it = polyhedra.begin(),
        poly_end = polyhedra.end();
      poly_it != poly_end; ++poly_it) {
    delete poly_it->polyhedron_ptr;
  }
  polyhedra.clear();
}

bool
Scene::open(QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"...").arg(filename);

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  QFileInfo fileinfo(filename);
  std::ifstream in(filename.toUtf8());
  if(!in) return false;

  Polyhedron* poly = new Polyhedron;
  in >> *poly;
  if(!in)
  {
    QMessageBox::critical(qobject_cast<QWidget*>(QObject::parent()),
                          tr("Cannot open file"),
                          tr("File %1 is not a valid OFF file.").arg(filename));
    QApplication::restoreOverrideCursor();
    return false;
  }
  poly->compute_normals();

  Polyhedron_entry entry;
  entry.polyhedron_ptr = poly;
  entry.name = fileinfo.baseName();
  entry.color=QColor(100, 100, 255);
  entry.activated = true;
  polyhedra.push_back(entry);
  cerr << " done.\n";
  QApplication::restoreOverrideCursor();

  emit updated_bbox();
  emit updated();
  QAbstractListModel::reset();

  return true;
}

void 
Scene::erase(int polyhedron_index)
{
  delete polyhedra[polyhedron_index].polyhedron_ptr;
  polyhedra.removeAt(polyhedron_index);

  emit updated();
  QAbstractListModel::reset();
}

void 
Scene::duplicate(int polyhedron_index)
{
  const Polyhedron_entry& entry = polyhedra[polyhedron_index];
  Polyhedron* poly = new Polyhedron(*entry.polyhedron_ptr);

  poly->compute_normals();

  Polyhedron_entry entry2;
  entry2.polyhedron_ptr = poly;
  entry2.name = QString("%1 (copy)").arg(entry.name);
  entry2.color=entry.color;
  entry2.activated = entry.activated;
  polyhedra.push_back(entry2);

  emit updated();
  QAbstractListModel::reset();
}

CGAL::Bbox_3 
Scene::bbox()
{
  if(polyhedra.empty()) {
    return CGAL::Bbox_3(0, 0, 0, 1, 1, 1);
  }
  else
  {
    Point p = polyhedra.begin()->polyhedron_ptr->vertices_begin()->point();
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedra::iterator 
          poly_it = polyhedra.begin(),
          poly_end = polyhedra.end();
        poly_it != poly_end; ++poly_it) {
      for(Polyhedron::Vertex_iterator
            v = poly_it->polyhedron_ptr->vertices_begin(),
            v_end = poly_it->polyhedron_ptr->vertices_end();
          v != v_end; ++v)
      {
        bbox = bbox + v->point().bbox();
      }
    }
    return bbox;
  }
}

void 
Scene::draw()
{
  for(int index = 0; index < polyhedra.size(); ++index)
  {
    Polyhedron_entry& entry = polyhedra[index];
    if(entry.activated) {
      std::cerr << "Drawing " << qPrintable(entry.name)
                << std::endl;
      Polyhedron* poly = entry.polyhedron_ptr;
      ::glEnable(GL_LIGHTING);
      if(index == selected_item) {
        CGALglcolor(entry.color.lighter(120));
      }
      else {
        CGALglcolor(entry.color);
      }
      poly->gl_draw_direct_triangles(false,
                                     true);
      if(index == selected_item) {
        CGALglcolor(entry.color.lighter(70));
      }
      else {
        CGALglcolor(entry.color.lighter(50));
      }
      ::glDisable(GL_LIGHTING);
      poly->superimpose_edges(true,false);
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
    return AfterLastColumn;
}

QVariant 
Scene::data(const QModelIndex &index, int role) const
{
  if (!index.isValid())
    return QVariant();
  
  switch(index.column())
  {
  case ColorColumn:
    if(role == Qt::DisplayRole || role == Qt::EditRole)
      return polyhedra.value(index.row()).color;
    else if(role == Qt::DecorationRole)
      return polyhedra.value(index.row()).color;
    break;
  case NameColumn:
    if(role == Qt::DisplayRole || role == Qt::EditRole)
      return polyhedra.value(index.row()).name;
    break;
  case ActivatedColumn:
    if(role == Qt::DisplayRole || role == Qt::EditRole)
      return polyhedra.value(index.row()).activated;
  default:
    return QVariant();
  }
  return QVariant();
}

QVariant 
Scene::headerData ( int section, Qt::Orientation orientation, int role ) const
{
  if (role == Qt::DisplayRole && orientation == Qt::Horizontal)
  {
    switch(section)
    {
    case NameColumn:
      return tr("Name");
      break;
    case ColorColumn:
      return tr("Color");
      break;
    case ActivatedColumn:
      return tr("Activated");
      break;
    default:
      return QVariant();
    }
  }
  return QAbstractListModel::headerData(section, orientation, role);
}

Qt::ItemFlags 
Scene::flags ( const QModelIndex & index ) const
{
  if (index.isValid()) {
    return QAbstractListModel::flags(index) | Qt::ItemIsEditable;
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
  if( role != Qt::EditRole || !index.isValid() )
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
  case ActivatedColumn:
    entry.activated = value.toBool();
    emit dataChanged(index, index);
    return true;
  default:
    return false;
  }
  return false;
}
