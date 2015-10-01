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

    Bbox bbox() const;

    Scene_group_item* clone() const {return 0;}
    //! Indicate if rendering mode is supported.
    bool supportsRenderingMode(RenderingMode m) const;

    QString toolTip() const;

    void addChild(Scene_item* new_item);

    QList<Scene_item*> getChildren() const {return children;}
private:
    QList<Scene_item*> children;

}; //end of class Scene_group_item


#endif // SCENE_GROUP_ITEM_H
