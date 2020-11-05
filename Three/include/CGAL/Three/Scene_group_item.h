// Copyright (c) 2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno


#ifndef SCENE_GROUP_ITEM_H
#define SCENE_GROUP_ITEM_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Scene_interface.h>
using namespace CGAL::Three;

#include <QtCore/qglobal.h>
#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
namespace CGAL {
namespace Three {
//!A Scene_group_item is a special Scene_item that does not draw anything,
//! but regroups other items as its children. It allows the
//! user to apply several actions to multiple items at the same time.
//! A custom Scene_item can derive from it to have children. They appear
//! hierarchically in the Geometric Objects list.
class DEMO_FRAMEWORK_EXPORT Scene_group_item : public Scene_item_rendering_helper
{
    Q_OBJECT
public :
    Scene_group_item(QString name = QString("New group"));
    ~Scene_group_item() {}
    //!Returns false to avoid disturbing the BBox of the scene.
    bool isFinite() const Q_DECL_OVERRIDE;
    //!Returns true to avoid disturbing the BBox of the scene.
    bool isEmpty() const Q_DECL_OVERRIDE;
    /*!
         * \brief Locks a child
         *
         * A locked child cannot be moved out of the group nor can it be deleted.
         * Use it to prevent a child to be destroyed without its parent.
         */
        void lockChild(CGAL::Three::Scene_item*);
        /*!
        * \brief Locks a child
        *
        * A locked child cannot be moved out of the group nor can it be deleted.
        * Use it to prevent a child to be destroyed without its parent.
        */
        void lockChild(Scene_interface::Item_id id);

        /*!
         * \brief Unlocks a child
         *
         * @see lockChild()
         */
        void unlockChild(CGAL::Three::Scene_item*);
        /*!
         * \brief Unlocks a child
         *
         * @see lockChild()
         */
        void unlockChild(Scene_interface::Item_id id);
        /*!
         * \brief Tells if a child is locked.
         * \return true if the child is locked.
         * @see lockChild()
         */
        bool isChildLocked(CGAL::Three::Scene_item*);
        /*!
             * \brief Tells if a child is locked.
             * \return true if the child is locked.
             * @see lockChild()
             */
        bool isChildLocked(Scene_interface::Item_id id);
    //!Returns if the group_item is currently expanded or collapsed in the Geometric Objects list.
    //! True means expanded, false means collapsed.
    //! @see isExpanded().
    bool isExpanded() const;
    //!Makes the group_item expanded or collapsed in the view.
    //! True means expanded, false means collapsed.
    //! @see isExpanded().
    void setExpanded(bool);
    //!Returns an empty Bbox to avoid disturbing the Bbox of the scene.
    Bbox bbox() const Q_DECL_OVERRIDE;
    //!Not supported.
    Scene_item* clone() const Q_DECL_OVERRIDE {return 0;}
    //! Indicates if the rendering mode is supported.
    //! \returns true for all rendering modes that are shared by
    //! all of the children.
    bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE;
    //!\returns a string containing the number of children.
    QString toolTip() const Q_DECL_OVERRIDE;

    /// Draw functions
    /// Scene_group_item's children are not drawn by the scene, they are drawn by the group.
    ///@{
    //!\brief draws all the children
    //!
    //! Calls `Scene_item::draw()`, then calls `Scene_item::drawEdges()`
    //! and `Scene_item::drawPoints` for each child if its current
    //! rendering mode is adequat.
    //! @see #RenderingMode
    virtual void draw(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    //!\brief draws all the children
    //!
    //! Calls `Scene_item::drawEdges()`, then calls `Scene_item::draw()`
    //! and `Scene_item::drawPoints` for each child if its current
    //! rendering mode is adequat.
    //! @see #RenderingMode
    virtual void drawEdges(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    //!\brief draws all the children
    //!
    //! Calls `Scene_item::drawPoints()`, then calls `Scene_item::draw()`
    //! and `Scene_item::drawEdges()` for each child if its current
    //! rendering mode is adequat.
    //! @see #RenderingMode
    virtual void drawPoints(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    //!
       //! \brief deals with the rendering, selecting and picking of
       //! the group's children.
       //!
       //! \param picked_item_IDs the depth-index map
       //! \param picked_pixel the screen point that has been picked.
       //! \param with_names should be `true` if a picking is being performed.
       //!
       virtual void renderChildren(Viewer_interface *,
                 QMap<float, int>& picked_item_IDs, const QPoint &picked_pixel,
                 bool with_names);

    ///@}

    //!Adds a CGAL::Three::Scene_item* to the list of children.
    //!@see getChildren() @see removeChild()
    void addChild(Scene_item* new_item);
    //!Adds a CGAL::Three::Scene_item* to the list of children from its id.
    //! //!@see getChildren() @see removeChild()
    void addChild(Scene_interface::Item_id new_id);
    //! \brief getChild gives access to the Scene_item of the wanted index.
    //! \returns the `CGAL::Three::Scene_item` which index is `id`.
    //!
    Scene_item* getChild(Scene_interface::Item_id id) { return scene->item(id);}
    Scene_item* getChild(Scene_interface::Item_id id) const{ return scene->item(id);}
    //!Sets all the children to the specified color.
    void setColor(QColor c) Q_DECL_OVERRIDE;
    //!Sets all the children in the specified rendering mode.
    void setRenderingMode(RenderingMode m) Q_DECL_OVERRIDE;
    //!Sets all the children to the specified visibility.
    void setVisible(bool b) Q_DECL_OVERRIDE;
    //!Sets all the children in points mode.
    void setPointsMode() {
      setRenderingMode(Points);
    }
    //!Sets all the children in wireframe rendering.
    void setWireframeMode() {
      setRenderingMode(Wireframe);
    }
    //!Sets all the children in wireframe rendering.
    void setWireframe() {
      setRenderingMode(Wireframe);
    }
    //!Sets all the children in flat rendering.
    void setFlat() {
      setRenderingMode(Flat);
    }
    //!Sets all the children in flat rendering.
    void setFlatMode() {
      setRenderingMode(Flat);
    }
    //!Sets all the children in flat rendering with edges.
    void setFlatPlusEdgesMode() {
      setRenderingMode(FlatPlusEdges);
    }
    //!Sets all the children in smooth rendering.
    void setGouraudMode() {
      setRenderingMode(Gouraud);
    }
    //!Sets all the children in point rendering with normals.
    void setPointsPlusNormalsMode(){
      setRenderingMode(PointsPlusNormals);
    }
    //!Sets the alpha value for the froup and all its children.
        virtual void setAlpha(int) Q_DECL_OVERRIDE;

    //! \brief Returns a list of all the direct children.
    //!
    //! Only returns children that have this item as a parent.
    //! Children of these children are not returned.
    QList<Scene_interface::Item_id> getChildren() const {return children;}

    //! \brief getChildrenForSelection returns the list of
    //! children to select along with the group.
    //!
    //! When a `Scene_group_item` is added to the selection of the scene,
    //! this function defines which of its children will be added too.
    //! Typically overriden to allow applying an operation from the
    //! Operation menu only to the parent item and not to its children.
    virtual QList<Scene_interface::Item_id> getChildrenForSelection() const {return children;}
    //!Removes a Scene_item from the list of children.
    //!@see getChildren() @see addChild()
    void removeChild( Scene_item* item)
    {
     if(isChildLocked(item))
      return;
     update_group_number(item,0);
     item->moveToGroup(0);
     children.removeOne(scene->item_id(item));
    }
    //!Removes a Scene_item from the list of children using its index.
    //!@see getChildren() @see addChild()
    void removeChild( Scene_interface::Item_id id)
    {
      removeChild(scene->item(id));
    }
    //!Moves a child up in the list.
    void moveUp(int);
    //!Moves a child down in the list.
    void moveDown(int);

    void compute_bbox() const Q_DECL_OVERRIDE{};
public Q_SLOTS:
    //!\brief Redraws children.
    //!
    //! As each drawing function of a group draws all parts of its children,
    //! once any of these functions is called, we skip all drawing calls
    //! until `resetDraw()` is called. This keeps children from being
    //! drawn several times. It is automatically called at the end of the scene's
    //! `draw()` function.
    void resetDraw() { already_drawn = false;}
    //!
    //! \brief adjustIds maintains the list of children up to date when an item has been erased.
    //! \param removed_id the index of the item that has been erased.
    //!
    void adjustIds(Scene_interface::Item_id removed_id)
    {
      for(int i = 0; i < children.size(); ++i)
      {
        if(children[i] > removed_id)
          --children[i];
        else if(children[i] == removed_id)//child has been removed from the scene, it doesn't exist anymore.
        {
          children.removeAll(removed_id);
        }
      }
    }
private:
    void update_group_number(Scene_item* new_item, int n);
    bool expanded;
    mutable bool already_drawn;
protected:
    Scene_interface *scene;
    //!Contains a reference to all the children of this group.
    QList<Scene_interface::Item_id> children;

}; //end of class Scene_group_item

}
}

#endif // SCENE_GROUP_ITEM_H
