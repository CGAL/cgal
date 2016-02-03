//! \file Scene.h
#ifndef SCENE_H
#define SCENE_H
#include "config.h"
#include "Scene_config.h"
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_draw_interface.h>
#include <CGAL/Three/Viewer_config.h>

#include <QtOpenGL/qgl.h>
#include <QStandardItemModel>
#include <QString>
#include <QColor>
#include <QList>
#include <QMap>
#include <QItemDelegate>
#include <QPixmap>
#include <QItemSelection>
#include <QGLViewer/qglviewer.h>
#include <QDebug>
#include <iostream>
#include <cmath>
#include <boost/variant.hpp>
#include <CGAL/Three/Scene_group_item.h>
class QEvent;
class QMouseEvent;
namespace GlSplat { class SplatRenderer; }
namespace CGAL { namespace Three{ class Viewer_interface;}}


class SCENE_EXPORT Scene  :
  public QStandardItemModel, public CGAL::Three::Scene_interface, public CGAL::Three::Scene_draw_interface

{
  Q_OBJECT
  Q_PROPERTY(int numberOfEntries READ numberOfEntries)

  friend class SceneDelegate;

public:
  QList<QModelIndex> getModelIndexFromId(int id) const;
  int getIdFromModelIndex(QModelIndex modelId) const;
  enum Columns { NameColumn = 0, 
                 ColorColumn, 
                 RenderingModeColumn, 
                 VisibleColumn,
                 ABColumn,
                 LastColumn = ABColumn,
                 NumberOfColumns = LastColumn + 1};
  Scene(QObject*  parent);
  ~Scene();

  //!Adds item to the items list, gives it an ID and
  //!updates the bounding box if needed.
  int addItem(CGAL::Three::Scene_item* item);

  void changeGroup(CGAL::Three::Scene_item* item, CGAL::Three::Scene_group_item* target_group);
  //!Sets item as the item at index and calls @ref Scene_item#changed().
  //!If emit_item_about_to_be_destroyed is set to true, emits
  //!an itemAboutToBeDestroyed signal.
  CGAL::Three::Scene_item* replaceItem(int index, CGAL::Three::Scene_item* item, bool emit_item_about_to_be_destroyed = false);
  /*! Deletes the item with the target index.
   * @returns  the index of the polyhedra just before the
   * one that is erased, or just after. -1 if
   * the list is empty.
   */
  Q_INVOKABLE int erase(int);  

  /*! Deletes the items with the target indexes.
   * @returns the index of the polyhedra just before the
   * one that is erased, or just after. Returns -1 if
   * the list is empty.
   */
  int erase(QList<int>);

  /*! Duplicate a scene item.
   * @returns the ID of the new item (-1 on error).
   */
  int duplicate(int index); 

  // Accessors (getters)
  //! @returns the number of items.
  int numberOfEntries() const;
  //! @returns the list of items.
  const QList<CGAL::Three::Scene_item*>& entries() const { return m_entries; }
  //! @returns the item at the target index.
  Q_INVOKABLE CGAL::Three::Scene_item* item(int) const ;
  //! @returns the id of the target item.
  Item_id item_id(CGAL::Three::Scene_item*) const;
  
  //! \todo Replace Index based selection functionality with those
  //! functions.
  ///@{
  CGAL::Three::Scene_item* selectedItem() const;
  QList<CGAL::Three::Scene_item*> selectedItems() const;
  QList<CGAL::Three::Scene_item*> selectionA() const;
  QList<CGAL::Three::Scene_item*> selectionB() const;
  ///@}

  //!@returns the currently selected item's index.
  int mainSelectionIndex() const;
  //!@returns the list of currently selected items indices.
  QList<int> selectionIndices() const;
  //!@returns the index of the Item_A
  int selectionAindex() const;
  //!@returns the index of the Item_B
  int selectionBindex() const;

  /*! Is called by Viewer::initializeGL(). Allows all the initialization
   * of OpenGL code that needs a context.
   */
  void initializeGL();
  /*! Is called by Viewer::draw(). Is deprecated and does nothing.*/
  void draw();
  /*! Is deprecated and does nothing.*/
  void drawWithNames();
  /*! Is called by Viewer::draw(Viewer_interface*). Calls draw_aux(false, viewer).
   * @see draw_aux(bool with_names, Viewer_interface).*/
  void draw(CGAL::Three::Viewer_interface*);
  /*! Is called by Viewer::drawWithNames(Viewer_interface*). Calls draw_aux(true, viewer).
   * @see draw_aux(bool with_names, Viewer_interface).*/
  void drawWithNames(CGAL::Three::Viewer_interface*);
  /*! Manages the key events.
   * @returns true if the keyEvent executed well.
   */
  bool keyPressEvent(QKeyEvent* e);

  //!@returns the scene bounding box
  Bbox bbox() const;
  float get_bbox_length() const;
  //!@returns the length of the bounding box's diagonal.
  double len_diagonal() const
  {
    Bbox box = bbox();
    double dx = box.xmax - box.xmin;
    double dy = box.ymax - box.ymin;
    double dz = box.zmax - box.zmin;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }


  // QStandardItemModel functions
  void moveRowUp();
  void moveRowDown();
  bool dropMimeData(const QMimeData *, Qt::DropAction, int, int, const QModelIndex &parent);
  QVariant data ( const QModelIndex & index, int role = ::Qt::DisplayRole ) const;
  //!@returns the type of data correspondind to the role.
  QVariant headerData ( int section, ::Qt::Orientation orientation, int role = ::Qt::DisplayRole ) const;
  //!@returns the flags for the item at the target index.
  ::Qt::ItemFlags flags ( const QModelIndex & index ) const;
  /*! Sets the column data for the target index. Returns false if index is not valid and
   * if role is not EditRole.*/
  bool setData(const QModelIndex &index, const QVariant &value, int role);
  QList<CGAL::Three::Scene_group_item*> group_entries() const ;
  QList<CGAL::Three::Scene_item*> item_entries() const ;
  // auxiliary public function for QMainWindow
  //!Selects the row at index i in the sceneView.
  QItemSelection createSelection(int i);
  //!Selects all the rows in the sceneView.
  QItemSelection createSelectionAll();

public Q_SLOTS:
  //!Specifies a group as Expanded for the view
  void setExpanded(QModelIndex);
  //!Specifies a group as Collapsed for the view
  void setCollapsed(QModelIndex);
  /*! This is an overloaded function.
   * Notifies the scene that the sender item was modified.
   * Called by the items. Calls @ref Scene_item#changed().
   * This function is called by the items.*/
  void itemChanged();
  /*! Notifies the scene that the item at index i was modified.
   * Called by the items. Calls @ref Scene_item#changed().
   * This function is called by the items.*/
  void itemChanged(int i); 
  /*! Notifies the scene that the item was modified.
   *  Calls @ref Scene_item#changed().
   * This function is called by the items.*/
  void itemChanged(CGAL::Three::Scene_item*);
  //!Removes item from all the groups of the scene.
  void remove_item_from_groups(CGAL::Three::Scene_item* item);

  void add_group(Scene_group_item* group);
  //!Re-organizes the sceneView.
  void group_added();
  //! Sets the selected item to the target index.
  void setSelectedItemIndex(int i)
  {
    selected_item = i;
  }
  //! Sets the selected item to the target index and emits selectionChanged(i).
  void setSelectedItem(int i )
  {
    selected_item = i;
    Q_EMIT selectionChanged(i);
  }

  //! Sets the target item as selected and emits setSelectedItem for its index.
  void setSelectedItem(CGAL::Three::Scene_item* item_to_select)
  {
    int i=0;
    Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
    {
      if (item==item_to_select)
      {
        Q_EMIT setSelectedItem(i);
        break;
      }
      ++i;
    }
  }
  //! Sets the target list of indices as the selected indices.
  QList<int> setSelectedItemsList(QList<int> l )
  {
    Q_FOREACH(int i,l)
    {
       CGAL::Three::Scene_group_item* group =
               qobject_cast<CGAL::Three::Scene_group_item*>(item(i));
       if(group)
       {
         QList<int> list;
         Q_FOREACH(CGAL::Three::Scene_item* child, group->getChildren())
           list<<m_entries.indexOf(child);
         l << setSelectedItemsList(list);
       }

    }
    selected_items_list = l;
    return l;
  }

  // Accessors (setters)
  //!Sets the item at index i to visible or not visible.
  void setItemVisible(int, bool b);
  //!Sets the item_A as the item at index i .
  void setItemA(int i);
  //!Sets the item_B as the item at index i .
  void setItemB(int i);

Q_SIGNALS:
  //generated automatically by moc
  void itemPicked(const QModelIndex &);
  void newItem(int);
  void updated_bbox();
  void updated();
  void itemAboutToBeDestroyed(CGAL::Three::Scene_item*);
  void selectionRay(double, double, double, double, double, double);
  void selectionChanged(int i);
  void restoreCollapsedState();
private Q_SLOTS:
  //! Casts a selection ray and calls the item function select.
  void setSelectionRay(double, double, double, double, double, double);
  void callDraw(){  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin(); viewer->update();}

private:
  /*! Calls the drawing functions of each visible item according
   * to its current renderingMode. If with_names is true, uses
   * the OpenGL mode GL_WITH_NAMES, essentially used for the picking.*/
  void draw_aux(bool with_names, CGAL::Three::Viewer_interface*);
  //! Re-draw the hierarchy of the view.
  void organize_items(CGAL::Three::Scene_item* item, QStandardItem *root, int loop);
  //! List of Scene_items.
  typedef QList<CGAL::Three::Scene_item*> Entries;
  //!List containing all the scene_items.
  Entries m_entries;
  //! Index of the currently selected item.
  int selected_item;
  //!List containing all the scene_group_items.
  QList<CGAL::Three::Scene_group_item*> m_group_entries;
  //!List of indices of the currently selected items.
  QList<int> selected_items_list;
  //!Index of the item_A.
  int item_A;
  //!Index of the item_B.
  int item_B;
  bool picked;
  bool gl_init;
  static GlSplat::SplatRenderer* ms_splatting;
  static int ms_splattingCounter;
  QMap<QModelIndex, int> index_map;

public:
  static GlSplat::SplatRenderer* splatting();

}; // end class Scene

class QAbstractProxyModel;
/*!
 * \brief The SceneDelegate class
 * Handles the columns of the sceneView
 */
class SCENE_EXPORT SceneDelegate : public QItemDelegate
{
public:
  SceneDelegate(QObject * parent = 0)
    : QItemDelegate(parent),
      checkOnPixmap(":/cgal/icons/check-on.png"),
      checkOffPixmap(":/cgal/icons/check-off.png")
  {
  }
//! Handles the clicks on the sceneView
  bool editorEvent(QEvent *event, QAbstractItemModel *model,
                   const QStyleOptionViewItem &option,
                   const QModelIndex &index);
  //! Draws the content of the sceneView
  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const;
  void setProxy(QAbstractProxyModel* p_proxy){
      proxy = p_proxy;
  }
  void setScene(Scene* p_scene){
      scene = p_scene;
  }

private:
  QPixmap checkOnPixmap;
  QPixmap checkOffPixmap;
  QAbstractProxyModel *proxy;
  Scene *scene;
  mutable int size;
}; // end class SceneDelegate

#endif // SCENE_H


