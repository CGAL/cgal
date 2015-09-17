#ifndef SCENE_H
#define SCENE_H
#include "config.h"
#include "Scene_config.h"

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
#include <QGLViewer/qglviewer.h>

#include <iostream>
#include <cmath>
#include <boost/variant.hpp>

class QEvent;
class QMouseEvent;
namespace GlSplat { class SplatRenderer; }

class Viewer_interface;

class SCENE_EXPORT Scene  :
  public QAbstractListModel, public Scene_interface, public Scene_draw_interface
{
  Q_OBJECT
  Q_PROPERTY(int numberOfEntries READ numberOfEntries)

  friend class SceneDelegate;

public:
  enum Columns { NameColumn = 0, 
                 ColorColumn, 
                 RenderingModeColumn, 
                 VisibleColumn,
                 ABColumn,
                 LastColumn = ABColumn,
                 NumberOfColumns = LastColumn + 1};

  Scene(QObject*  parent);
  ~Scene();

  int addItem(Scene_item* item);
  Scene_item* replaceItem(int index, Scene_item* item, bool emit_item_about_to_be_destroyed = false);

  Q_INVOKABLE int erase(int);  
  int erase(QList<int>);  
  // Returns the index of the polyhedra just before the
  // one that is erased, or just after. Returns -1 if
  // the list is empty.

  // Duplicate a scene item. Return the ID of the new item (-1 on error).
  int duplicate(int index); 

  // Accessors (getters)
  int numberOfEntries() const;
  const QList<Scene_item*>& entries() const { return m_entries; }
  Q_INVOKABLE Scene_item* item(int) const ;
  Item_id item_id(Scene_item*) const;
  
  //! \todo Replace Index based selection functionality with those
  //! functions.
  ///@{
  Scene_item* selectedItem() const;
  QList<Scene_item*> selectedItems() const;
  QList<Scene_item*> selectionA() const;
  QList<Scene_item*> selectionB() const;
  ///@}

  int mainSelectionIndex() const;
  QList<int> selectionIndices() const;
  int selectionAindex() const;
  int selectionBindex() const;

  // initializeGL() is called by Viewer::initializeGL()
  void initializeGL();
  // draw() is called by Viewer::draw()
  void draw();
  void drawWithNames();
  void draw(Viewer_interface*);
  void drawWithNames(Viewer_interface*);

  bool keyPressEvent(QKeyEvent* e);

  // Get scene bounding box
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
  QItemSelection createSelectionAll();

public Q_SLOTS:
  // Notify the scene that an item was modified
  void itemChanged(); // slots called by items themself
  void itemChanged(int i); 
  void itemChanged(Scene_item*);

  void setSelectedItemIndex(int i)
  {
    selected_item = i;
  }

  void setSelectedItem(int i )
  {
    selected_item = i;
    Q_EMIT selectionChanged(i);
  };

  void setSelectedItem(Scene_item* item_to_select)
  {
    int i=0;
    Q_FOREACH(Scene_item* item, m_entries)
    {
      if (item==item_to_select)
      {
        Q_EMIT setSelectedItem(i);
        break;
      }
      ++i;
    }
  };

  void setSelectedItemsList(QList<int> l )
  {
    selected_items_list = l;
  };

  // Accessors (setters)
  void setItemVisible(int, bool b);
  void setItemA(int i);
  void setItemB(int i);

Q_SIGNALS:
  void newItem(int);
  void updated_bbox();
  void updated();
  void itemAboutToBeDestroyed(Scene_item*);
  void selectionRay(double, double, double, double, double, double);
  void selectionChanged(int i);

private Q_SLOTS:
  void setSelectionRay(double, double, double, double, double, double);
  void callDraw(){  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin(); viewer->update();}

private:
  void draw_aux(bool with_names, Viewer_interface*);
  typedef QList<Scene_item*> Entries;
  Entries m_entries;
  int selected_item;
  QList<int> selected_items_list;
  int item_A;
  int item_B;
  static GlSplat::SplatRenderer* ms_splatting;
  static int ms_splattingCounter;
public:
  static GlSplat::SplatRenderer* splatting();

}; // end class Scene

class SCENE_EXPORT SceneDelegate : public QItemDelegate
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
