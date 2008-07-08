#ifndef SCENE_H
#define SCENE_H

#include <QAbstractListModel>
#include <QString>
#include <QColor>
#include <QList>

// CGAL
#include <CGAL/basic.h>

// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// surface
#include <CGAL/Polyhedron_3.h>
#include <CGAL/enriched_polyhedron.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle_3;

class Scene  :
  public QAbstractListModel
{
  Q_OBJECT
public:
  typedef Enriched_polyhedron<Kernel,
                              Enriched_items> Polyhedron;
  typedef Kernel::Point_3 Point_3;

private:
  struct Polyhedron_entry {
    Polyhedron* polyhedron_ptr;
    QString name;
    QColor color;
    bool activated;
  };

public:
  Scene();
  ~Scene();

  bool open(QString);
  void erase(int);
  void duplicate(int);

  void draw();
  CGAL::Bbox_3 bbox();

  enum Columns { NameColumn = 0, ColorColumn, ActivatedColumn,
                 AfterLastColumn = ActivatedColumn + 1};

  // QAbstractItemModel functions
  int rowCount ( const QModelIndex & parent = QModelIndex() ) const;
  int columnCount ( const QModelIndex & parent = QModelIndex() ) const;
  QVariant data ( const QModelIndex & index, int role = Qt::DisplayRole ) const;
  QVariant headerData ( int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const;
  Qt::ItemFlags flags ( const QModelIndex & index ) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role);

signals:
  void updated_bbox();
  void updated();

private:
  typedef QList<Polyhedron_entry> Polyhedra;
  Polyhedra polyhedra;
}; // end class Scene

#endif // SCENE_H
