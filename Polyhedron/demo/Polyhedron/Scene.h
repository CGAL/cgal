#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <QAbstractListModel>
#include <QString>
#include <QColor>
#include <QList>
#include <QItemDelegate>
#include <QPixmap>
#include <QItemSelection>

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

struct Polyhedron : public Enriched_polyhedron<Kernel,
                                               Enriched_items> {};

class QEvent;
class QMouseEvent;

class Scene  :
  public QAbstractListModel
{
  Q_OBJECT

private:
  struct Polyhedron_entry {
    Polyhedron* polyhedron_ptr;
    QString name;
    QColor color;
    bool activated;
    RenderingMode rendering_mode;
  };

public:
  Scene(QObject*  parent);
  ~Scene();

  int open(QString);  // Returns the index of the new polyhedra (-1 if
                      // error)

  int erase(int);     // Returns the index of the polyhedra just before the
                      // one that is erased, or just after. Returns -1 if
                      // the list is empty.

  int duplicate(int); // Returns the index of the new polyhedra

  inline QItemSelection createSelection(int i);

  // TODO: move that elsewhere (in MainWindow)
  void convex_hull(int);
  void simplify(int);

  inline Polyhedron* getPolyhedron(int);

  void draw();
  CGAL::Bbox_3 bbox();

  enum RenderingMode { Fill = 0, Wireframe, LastRenderingMode = Wireframe,
                       NumberOfRenderingMode = Wireframe+1};
  enum Columns { NameColumn = 0, ColorColumn, ActivatedColumn,
                 LastColumn = ActivatedColumn,
                 NumberOfColumns = LastColumn + 1};

  // QAbstractItemModel functions
  int rowCount ( const QModelIndex & parent = QModelIndex() ) const;
  int columnCount ( const QModelIndex & parent = QModelIndex() ) const;
  QVariant data ( const QModelIndex & index, int role = Qt::DisplayRole ) const;
  QVariant headerData ( int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const;
  Qt::ItemFlags flags ( const QModelIndex & index ) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role);

public slots:
  void setSelectedItem(int i )
  {
    selected_item = i;
  };

signals:
  void updated_bbox();
  void updated();

private:
  typedef QList<Polyhedron_entry> Polyhedra;
  Polyhedra polyhedra;
  int selected_item;
}; // end class Scene

class SceneDelegate : public QItemDelegate
{
public:
  SceneDelegate(QObject * parent = 0)
    : QItemDelegate(parent),
      checkOnPixmap(":/cgal/icons/check-on.png"),
      checkOffPixmap(":/cgal/icons/check-off.png")
  {
  }

  bool editorEvent(QEvent *event, QAbstractItemModel *model,
                   const QStyleOptionViewItem &option,
                   const QModelIndex &index);
  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const;

private:
  QPixmap checkOnPixmap;
  QPixmap checkOffPixmap;
  mutable int size;
}; // end class SceneDelegate

/*
 * Inlined functions
 */ 

Polyhedron* Scene::getPolyhedron(int index)
{
  if( index < 0 || index >= polyhedra.size() )
    return 0;
  else 
    return polyhedra[index].polyhedron_ptr;
}

QItemSelection Scene::createSelection(int i)
{
  return QItemSelection(QAbstractItemModel::createIndex(i, 0),
                        QAbstractItemModel::createIndex(i, LastColumn));
}

#endif // SCENE_H
