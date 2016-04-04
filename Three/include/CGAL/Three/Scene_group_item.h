// Copyright (c) 2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Maxime Gimeno


#ifndef SCENE_GROUP_ITEM_H
#define SCENE_GROUP_ITEM_H

#include <CGAL/Three/Scene_item.h>
using namespace CGAL::Three;

#include <QtCore/qglobal.h>
#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
namespace CGAL {
namespace Three {
//!A Scene_group_item is a virtual Scene_item that does not draw anything,
//! but regroups other items as its children. It allows the
//! user to apply several actions to multiple items at the same time.
class DEMO_FRAMEWORK_EXPORT Scene_group_item : public Scene_item
{
    Q_OBJECT
public :
    Scene_group_item(QString name = QString("New group"));
    ~Scene_group_item() {}
    //!Returns false to avoid disturbing the BBox of the scene.
    bool isFinite() const;
    //!Returns true to avoid disturbing the BBox of the scene.
    bool isEmpty() const ;
    //!Returns if the group_item is currently expanded or collapsed in the view.
    //! True means expanded, false means collapsed.
    //! @see setExpanded.
    bool isExpanded() const;
    //!Makes the group_item expanded or collapsed in the view.
    //! True means expanded, false means collapsed.
    void setExpanded(bool);
    //! @see isExpanded.
    //!Returns an empty BBox to avoid disturbing the BBox of the scene.
    Bbox bbox() const;
    //!Not supported.
    Scene_item* clone() const {return 0;}
    //! Indicate if rendering mode is supported.
    bool supportsRenderingMode(RenderingMode m) const;
    //!Prints the number of children.
    QString toolTip() const;
    //!Adds a Scene_item* to the list of children.
    //!@see getChildren. @see removeChild.
    void addChild(Scene_item* new_item);
    //!Sets all the children to the specified color.
    void setColor(QColor c);
    //!Sets all the children in the specified rendering mode.
    void setRenderingMode(RenderingMode m);
    //!Sets all the children to the specified visibility.
    void setVisible(bool b);
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
    //!Sets all the children in splat rendering.
    void setSplattingMode(){
      setRenderingMode(Splatting);
    }
    //!Returns a list of all the children.
    QList<Scene_item*> getChildren() const {return children;}
    //!Removes a Scene_item from the list of children.
    //!@see getChildren @see addChild
    void removeChild( Scene_item* item)
    {
      Scene_group_item* group =
              qobject_cast<Scene_group_item*>(item);
      if(group)
        Q_FOREACH(Scene_item* child, group->getChildren())
            removeChild(child);
      item->has_group=0;
      children.removeOne(item);
    }
    //!Moves a child up in the list.
    void moveUp(int);
    //!Moves a child down in the list.
    void moveDown(int);

private:
    //!Contains a reference to all the children of this group.
    QList<Scene_item*> children;
    //!Updates the property has_group for each group and sub-groups containing new_item.
    void add_group_number(Scene_item*new_item);
    bool expanded;

}; //end of class Scene_group_item

}
}

#endif // SCENE_GROUP_ITEM_H
