#ifndef SCENE_GROUP_ITEM_H
#define SCENE_GROUP_ITEM_H

#include "Scene_item.h"
#include "Messages_interface.h"

class Q_DECL_EXPORT Scene_group_item : public Scene_item
{
    Q_OBJECT
public :
    Scene_group_item(QString name = QString());
    ~Scene_group_item() {}
    bool isFinite() const;

    bool isEmpty() const ;

    bool isExpanded() const;

    void setExpanded(bool);

    Bbox bbox() const;

    Scene_group_item* clone() const {return 0;}
    //! Indicate if rendering mode is supported.
    bool supportsRenderingMode(RenderingMode m) const;

    QString toolTip() const;

    void addChild(Scene_item* new_item);

    void setColor(QColor c);

    void setRenderingMode(RenderingMode m);

    void setVisible(bool b);

    void setPointsMode() {
      setRenderingMode(Points);
    }

    void setWireframeMode() {
      setRenderingMode(Wireframe);
    }
    void setWireframe() {
      setRenderingMode(Wireframe);
    }

    void setFlat() {
      setRenderingMode(Flat);
    }
    void setFlatMode() {
      setRenderingMode(Flat);
    }

    void setFlatPlusEdgesMode() {
      setRenderingMode(FlatPlusEdges);
    }

    void setGouraudMode() {
      setRenderingMode(Gouraud);
    }

    void setPointsPlusNormalsMode(){
      setRenderingMode(PointsPlusNormals);
    }

    void setSplattingMode(){
      setRenderingMode(Splatting);
    }

    QList<Scene_item*> getChildren() const {return children;}

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

Q_SIGNALS:
    void updated(int row, int column);

private:
    QList<Scene_item*> children;
    void add_group_number(Scene_item*new_item);
    bool expanded;

}; //end of class Scene_group_item


#endif // SCENE_GROUP_ITEM_H
