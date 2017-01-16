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
  int addItem(CGAL::Three::Scene_item* item) Q_DECL_OVERRIDE;

  void changeGroup(CGAL::Three::Scene_item* item, CGAL::Three::Scene_group_item* target_group) Q_DECL_OVERRIDE;
  //!Sets item as the item at index and calls @ref Scene_item#changed().
  //!If emit_item_about_to_be_destroyed is set to true, emits
  //!an itemAboutToBeDestroyed signal.
  CGAL::Three::Scene_item* replaceItem(int index, CGAL::Three::Scene_item* item, bool emit_item_about_to_be_destroyed = false) Q_DECL_OVERRIDE;
  /*! Deletes the item with the target index.
   * @returns  the index of the polyhedra just before the
   * one that is erased, or just after. -1 if
   * the list is empty.
   */
  Q_INVOKABLE int erase(int) Q_DECL_OVERRIDE;

  /*! Deletes the items with the target indexes.
   * @returns the index of the polyhedra just before the
   * one that is erased, or just after. Returns -1 if
   * the list is empty.
   */
  int erase(QList<int>);

  /*! Duplicate a scene item.
   * @returns the ID of the new item (-1 on error).
   */
  int duplicate(int index) Q_DECL_OVERRIDE;

  // Accessors (getters)
  //! @returns the number of items.
  int numberOfEntries() const Q_DECL_OVERRIDE;
  //! @returns the list of items.
  const QList<CGAL::Three::Scene_item*>& entries() const { return m_entries; }
  //! @returns the item at the target index.
  Q_INVOKABLE CGAL::Three::Scene_item* item(int) const Q_DECL_OVERRIDE;
  //! @returns the id of the target item.
  int item_id(CGAL::Three::Scene_item*) const Q_DECL_OVERRIDE;
  
  //! \todo Replace Index based selection functionality with those
  //! functions.
  ///@{
  CGAL::Three::Scene_item* selectedItem() const;
  QList<CGAL::Three::Scene_item*> selectedItems() const;
  QList<CGAL::Three::Scene_item*> selectionA() const;
  QList<CGAL::Three::Scene_item*> selectionB() const;
  ///@}

  //!@returns the currently selected item's index.
  int mainSelectionIndex() const Q_DECL_OVERRIDE;
  //!@returns the list of currently selected items indices.
  QList<int> selectionIndices() const Q_DECL_OVERRIDE;
  //!@returns the index of the Item_A
  int selectionAindex() const Q_DECL_OVERRIDE;
  //!@returns the index of the Item_B
  int selectionBindex() const Q_DECL_OVERRIDE;

  /*! Is called by Viewer::initializeGL(). Allows all the initialization
   * of OpenGL code that needs a context.
   */
  void initializeGL() Q_DECL_OVERRIDE;
  /*! Sets the screen coordinates of the currently picked point.*/
  void setPickedPixel(const QPoint &p) Q_DECL_OVERRIDE {picked_pixel = p;}
  /*! Is called by Viewer::draw(Viewer_interface*). Calls draw_aux(false, viewer).
   * @see draw_aux(bool with_names, Viewer_interface).*/
  void draw(CGAL::Three::Viewer_interface*) Q_DECL_OVERRIDE;
  /*! Is called by Viewer::drawWithNames(Viewer_interface*). Calls draw_aux(true, viewer).
   * @see draw_aux(bool with_names, Viewer_interface).*/
  void drawWithNames(CGAL::Three::Viewer_interface*) Q_DECL_OVERRIDE;
  /*! Manages the key events.
   * @returns true if the keyEvent executed well.
   */
  bool keyPressEvent(QKeyEvent* e) Q_DECL_OVERRIDE;

  void printPrimitiveId(QPoint point,
                        CGAL::Three::Viewer_interface*) Q_DECL_OVERRIDE;
  void printPrimitiveIds(CGAL::Three::Viewer_interface*) Q_DECL_OVERRIDE;
  //!Re-computes the primitiveIds for `item`
  void updatePrimitiveIds(Viewer_interface *, Scene_item *item) Q_DECL_OVERRIDE;
  bool testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer) Q_DECL_OVERRIDE;
  //!@returns the scene bounding box
  Bbox bbox() const Q_DECL_OVERRIDE;
  float get_bbox_length() const Q_DECL_OVERRIDE;
  //!@returns the length of the bounding box's diagonal.
  double len_diagonal() const Q_DECL_OVERRIDE
  {
    Bbox box = bbox();
    double dx = box.xmax() - box.xmin();
    double dy = box.ymax() - box.ymin();
    double dz = box.zmax() - box.zmin();
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }


  // QStandardItemModel functions
  void moveRowUp();
  void moveRowDown();
  bool dropMimeData(const QMimeData *, Qt::DropAction, int, int,
                    const QModelIndex &parent) Q_DECL_OVERRIDE;
  QVariant data ( const QModelIndex & index,
                  int role = ::Qt::DisplayRole ) const Q_DECL_OVERRIDE;
  //!@returns the type of data correspondind to the role.
  QVariant headerData ( int section, ::Qt::Orientation orientation,
                        int role = ::Qt::DisplayRole ) const Q_DECL_OVERRIDE;
  //!@returns the flags for the item at the target index.
  ::Qt::ItemFlags flags ( const QModelIndex & index ) const Q_DECL_OVERRIDE;
  /*! Sets the column data for the target index. Returns false if index is not valid and
   * if role is not EditRole.*/
  bool setData(const QModelIndex &index, const QVariant &value,
               int role) Q_DECL_OVERRIDE;
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
  /*!
   *Calls itemChanged() on the sender if it's an item.

*/
  void itemChanged();
  /*! Notifies the scene that the item at index i was modified.
   * Called by the items. Calls @ref Scene_item#changed().
   * This function is called by the items.*/
  void itemChanged(int i) Q_DECL_OVERRIDE;
  /*! Notifies the scene that the item was modified.
   *  Calls @ref Scene_item#changed().
   * This function is called by the items.*/
  void itemChanged(CGAL::Three::Scene_item*) Q_DECL_OVERRIDE;
  //!Removes item from all the groups of the scene.
  void remove_item_from_groups(CGAL::Three::Scene_item* item);

  void addGroup(Scene_group_item* group) Q_DECL_OVERRIDE;
  //!Re-organizes the sceneView.
  void redraw_model();
  //! Sets the selected item to the target index. Emits a signal to notify
  //! that a new item is now selected.
  void setSelectedItemIndex(int i)
  {
    selected_item = i;
    Q_EMIT itemIndexSelected(i);
  }
  //! Clears the current selection then sets the selected item to the target index.
  //! Used to update the selection in the QTreeView.
  void setSelectedItem(int i ) Q_DECL_OVERRIDE
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
  //! Used to update the selected item in the QTreeView.
  void selectionChanged(int i);
  //! Used when you don't want to update the sleectedItem in the QTreeView.
  void itemIndexSelected(int i);
  void restoreCollapsedState();
  void drawFinished();
private Q_SLOTS:
  //! Casts a selection ray and calls the item function select.
  void setSelectionRay(double, double, double, double, double, double);
  void callDraw(){  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin(); viewer->update();}
  void s_itemAboutToBeDestroyed(CGAL::Three::Scene_item *);
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
  //!List of indices of the currently selected items.
  QList<int> selected_items_list;
  //!Index of the item_A.
  int item_A;
  //!Index of the item_B.
  int item_B;
  bool picked;
  QPoint picked_pixel;
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


