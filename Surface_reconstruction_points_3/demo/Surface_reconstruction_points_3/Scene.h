#ifndef SCENE_H
#define SCENE_H
#include "config.h"

#include "Scene_interface.h"
#include "Scene_draw_interface.h"

#include <QtOpenGL/qgl.h>
#include <QAbstractListModel>
#include <QString>
#include <QColor>
#include <QList>
#include <QItemDelegate>
#include <QPixmap>
#include <QItemSelection>

#include <iostream>
#include <cmath>
#include <boost/variant.hpp>

class QEvent;
class QMouseEvent;
namespace GlSplat { class SplatRenderer; }

class Scene  :
  public QAbstractListModel, public Scene_interface, public Scene_draw_interface
{
  Q_OBJECT

  friend class SceneDelegate;

public:
  enum Columns { NameColumn = 0,
                 ColorColumn,
                 RenderingModeColumn,
                 VisibleColumn,
                 LastColumn = VisibleColumn,
                 NumberOfColumns = LastColumn + 1};

  Scene(QObject*  parent);
  ~Scene();

  Item_id addItem(Scene_item* item);

  // Erases a scene item.
  // Returns the index of the polyhedra just before the one that is erased,
  //  or just after. Returns -1 if the list is empty.
  Item_id erase(Item_id index);

  // Duplicates a scene item. Returns the ID of the new item (-1 on error).
  Item_id duplicate(Item_id index);
  // Converts a polyhedron to a point set.
  // Returns the ID of the new item (-1 on error).
  Item_id convertToPointSet(Item_id index);

  // Delete selection in a scene item
  void deleteSelection(Item_id index);
  // Reset selection mark in a scene item.
  void resetSelection(Item_id index);

  // Accessors (getters)
  size_t numberOfEntries() const;
  Scene_item* item(Item_id) const ;
  Item_id mainSelectionIndex() const;

  // initializeGL() is called by Viewer::initializeGL()
  void initializeGL();
  // draw() is called by Viewer::draw()
  void draw();
  void drawWithNames();

  // Gets scene bounding box
  Bbox bbox() const;
  double len_diagonal() const
  {
    Bbox box = bbox();
    double dx = box.xmax - box.xmin;
    double dy = box.ymax - box.ymin;
    double dz = box.zmax - box.zmin;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

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
  // Notify the scene that an item was modified
  void itemChanged(Item_id i);
  void itemChanged(Scene_item*);

  void setSelectedItem(Item_id i )
  {
    selected_item = i;
  };

  // Accessors (setters)
  void setItemVisible(int, bool b);

signals:
  void updated_bbox();
  void updated();
  void itemAboutToBeDestroyed(Scene_item*);

private:
  void draw_aux(bool with_names);
  typedef QList<Scene_item*> Entries;
  Entries entries;
  int selected_item;
#ifdef CGAL_GLEW_ENABLED
  static GlSplat::SplatRenderer* ms_splatting;
  static int ms_splattingCounter;
public:
  static GlSplat::SplatRenderer* splatting();
#endif

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
