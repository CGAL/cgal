#ifndef POLYHEDRON_DEMO_GROUP_H
#define POLYHEDRON_DEMO_GROUP_H
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QList>
#include "Scene_group_item.h"
#include "CGAL/Three/Scene_interface.h"
#include "Scene.h"

class Scene_interface;
#include "CGAL/Three/Polyhedron_demo_plugin_helper.h"
class Polyhedron_demo_group_plugin  :
        public QObject,
        public CGAL::Three::Polyhedron_demo_plugin_helper
{
    //Configures CMake to use MOC correctly
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
    using Polyhedron_demo_plugin_helper::init;
    void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface*);

    bool applicable(QAction*) const
    {
        return true;
    }

    Scene* trueScene;
public Q_SLOTS:
    /*! Called by the action "add new group". Adds a new group to the Scene. If items were selected at the
     * moment this function is called, they are added to the group.
     * If all the selected items were in the same group, the new group
     * is added as a sub-group of this group.*/
    void add_group() {
        //Find the indices of the selected items
        QList<int> indices;
        QList<int> blacklist;
        Q_FOREACH(int id, scene->selectionIndices()){
            Scene_group_item* group =
                    qobject_cast<Scene_group_item*>(trueScene->item(id));
            if(group)
                Q_FOREACH(Scene_item *item, group->getChildren())
                    blacklist<<trueScene->item_id(item);

            if(!indices.contains(id) && !blacklist.contains(id))
                indices<<id;
}
        //checks if all the selected items are in the same group
        bool all_in_one = true;
        if(indices.isEmpty())
            all_in_one = false;
        // new group to create
        Scene_group_item * group = new Scene_group_item("new group");
        //group containing the selected item
        Scene_group_item * existing_group = 0;
        //for each selected item
        Q_FOREACH(int id, indices){
            //if the selected item is in a group
            if(trueScene->item(id)->has_group!=0){
                //for each group
                Q_FOREACH(Scene_group_item *item, trueScene->group_entries())
                {
                    //if the group contains the selected item
                    if(item->getChildren().contains(trueScene->item(id))){
                        //if it is the first one, we initialize existing_group
                        if(existing_group == 0)
                            existing_group = item;
                        //else we check if it is the same group as before.
                        //If not, all selected items are not in the same group
                        else if(existing_group != item)
                            all_in_one = false;
                        break;
                    }
                }//end for each group
            }
            //else it is impossible that all the selected items are in the same group
            else{
                all_in_one = false;
                break;
            }
        }//end foreach selected item

        //If all the selected items are in the same group, we put them in a sub_group of this group
        if(all_in_one)
        {
            Q_FOREACH(int id, indices)
                trueScene->changeGroup(trueScene->item(id),group);
            trueScene->changeGroup(group, existing_group);
            scene->addItem(group);
            trueScene->group_added();
        }
        //else wer create a new group
        else
        {
            Q_FOREACH(int id, indices)
                trueScene->changeGroup(trueScene->item(id),group);
            scene->addItem(group);
            trueScene->group_added();
        }

    }


private:
    //!Called when "Add new group" in the file menu is triggered.
    QAction* actionAddToGroup;
    void print_message(QString message) { messages->information(message); }
    Messages_interface* messages;

}; //end of class Polyhedron_demo_group_plugin
#endif
