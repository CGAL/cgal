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

#include "Polyhedron_type_fwd.h"

#include <iostream>

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
  bool save(int,QString); // Returns true upon successful save

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
  int selectionAindex() const;
  int selectionBindex() const;

  // for backward compatibility
  Polyhedron* getPolyhedron(int i) { return polyhedron(i); }

  // draw() is called by Viewer::draw()
  void draw();

  struct Bbox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
  };
  // defined in Scene_polyhedron_operations.cpp
  Bbox bbox();

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
  // functions that need to know the type Polyhedron
  // defined in Scene_polyhedron_operations.cpp
  Polyhedron* new_polyhedron();
  Polyhedron* copy_polyhedron(Polyhedron* poly);
  void destroy(Polyhedron*);
  bool load_polyhedron(Polyhedron* poly, std::istream& in); // return true
                                                            // iif the
                                                            // loading is ok.
  bool save_polyhedron(Polyhedron* poly, std::ostream& out); // return true
                                                             // iif the
                                                             // save is ok.

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
