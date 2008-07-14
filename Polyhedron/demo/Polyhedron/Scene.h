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
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;

// Boolean operations work only with exact kernel
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_polyhedron;

struct Polyhedron : public CGAL::Polyhedron_3<Kernel> {};

class QEvent;
class QMouseEvent;

class Scene  :
  public QAbstractListModel
{
  Q_OBJECT

  friend class SceneDelegate;
public:
  enum RenderingMode { Fill = 0, 
                       Wireframe, 
                       LastRenderingMode = Wireframe,
                       NumberOfRenderingMode = Wireframe+1};

  enum Columns { NameColumn = 0, 
                 ColorColumn, 
                 RenderingModeColumn, 
                 ActivatedColumn,
                 ABColumn,
                 LastColumn = ABColumn,
                 NumberOfColumns = LastColumn + 1};


  Scene(QObject*  parent);
  ~Scene();

  void addPolyhedron(Polyhedron* p,
                     QString name,
                     QColor color = defaultColor,
                     bool activated = true,
                     RenderingMode mode = Fill);

  int open(QString);  // Returns the index of the new polyhedra (-1 if
                      // error)

  int erase(int);     // Returns the index of the polyhedra just before the
                      // one that is erased, or just after. Returns -1 if
                      // the list is empty.

  int duplicate(int); // Returns the index of the new polyhedra

  // accessors
  Polyhedron* polyhedron(int i);
  QColor polyhedronColor(int);
  QString polyhedronName(int);
  bool isPolyhedronActivated(int);
  RenderingMode polyhedronRenderingMode(int);

  // for backward compatibility
  Polyhedron* getPolyhedron(int i) { return polyhedron(i); }

  // draw() is called by Viewer::draw()
  void draw();
  CGAL::Bbox_3 bbox();

  // QAbstractItemModel functions
  int rowCount ( const QModelIndex & parent = QModelIndex() ) const;
  int columnCount ( const QModelIndex & parent = QModelIndex() ) const;
  QVariant data ( const QModelIndex & index, int role = ::Qt::DisplayRole ) const;
  QVariant headerData ( int section, ::Qt::Orientation orientation, int role = ::Qt::DisplayRole ) const;
  ::Qt::ItemFlags flags ( const QModelIndex & index ) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role);

  // auxiliary public function for QMainWindow
  QItemSelection createSelection(int i);
public slots:
  void polyhedronChanged(int i);
  void polyhedronChanged(Polyhedron*);
  void setSelectedItem(int i )
  {
    selected_item = i;
  };

signals:
  void updated_bbox();
  void updated();

private:
  static const QColor defaultColor; // defined in Scene.cpp

  struct Polyhedron_entry {
    Polyhedron_entry() : rendering_mode(Fill) {};

    Polyhedron* polyhedron_ptr;
    QString name;
    QColor color;
    bool activated;
    RenderingMode rendering_mode;
  };

  typedef QList<Polyhedron_entry> Polyhedra;
  Polyhedra polyhedra;
  int selected_item;
  int item_A;
  int item_B;
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

#endif // SCENE_H
