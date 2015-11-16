#ifndef POLYHEDRON_DEMO_GROUP_H
#define POLYHEDRON_DEMO_GROUP_H
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QList>
#include <QAbstractListModel>
#include "Scene_group_item.h"
#include "CGAL_demo/Scene_interface.h"
#include "ui_Polyhedron_demo_group_plugin.h"
#include "Scene.h"

class Scene_interface;

class GroupModel : public QAbstractListModel
{

public:
    GroupModel(Scene* scene)
    {
        data_pool = scene;
    }

    QVariant data(const QModelIndex &index, int role) const;
    int rowCount(const QModelIndex &parent) const
    {
        Q_UNUSED (parent);
        return data_pool->group_entries().size();
    }
    QList<Scene_group_item*> groups() const
    {
        return data_pool->group_entries();
    }

    Qt::ItemFlags flags (const QModelIndex &index) const;
    bool setData(const QModelIndex &index, const QVariant &value, int role);

private:
    Scene* data_pool;
};

class ItemsModel : public QAbstractListModel
{
public:
    ItemsModel(Scene* scene)
    {
        data_pool = scene;
    }

    QVariant data(const QModelIndex &index, int role) const;

    int rowCount(const QModelIndex &parent = QModelIndex()) const
    {
        Q_UNUSED(parent);
        return data_pool->item_entries().size();
    }

    QList<Scene_item*> items() const
    {
        return data_pool->item_entries();
    }
private:
    Scene* data_pool;

};


#include "Polyhedron_demo_plugin_helper.h"
class Polyhedron_demo_group_plugin  :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    //Configures CMake to use MOC correctly
    Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
    using Polyhedron_demo_plugin_helper::init;
    void init(QMainWindow* mw, Scene_interface* sc, Messages_interface*);

    bool applicable(QAction*) const
    {
        return scene->numberOfEntries() > 0;
    }
    QList<QAction*> actions() const {
        return QList<QAction*>() << actionAddToGroup;
    }
    GroupModel *group_model;
    ItemsModel *items_model;
    void add_to_group() {
        Q_FOREACH(Scene_item* item, selected_items)
        {
            if(selected_group->getChildren().contains(item))
            {
                 //selected_group->removeChild(item);
                trueScene->remove_item_from_groups(item);
                item->has_group = 0;

            }
            else if(selected_group == item)
                 print_message("A group cannot contain itself.");
            else
            {
                trueScene->changeGroup(item, selected_group);

                if(!trueScene->item_entries().contains(selected_group))
                    scene->addItem(selected_group);
            }
        }
        selected_items.clear();
    }
    QList<Scene_item*> selected_items;
    Scene_group_item* selected_group;
    Scene* trueScene;
public Q_SLOTS:
    void GroupChoice();


private:
    QAction* actionAddToGroup;
    void print_message(QString message) { messages->information(message); }
    Messages_interface* messages;

}; //end of class Polyhedron_demo_group_plugin

class Polyhedron_demo_group_plugin_dialog : public QDialog, private Ui::GroupDialog
{
    Q_OBJECT
public:
    Polyhedron_demo_group_plugin_dialog(Polyhedron_demo_group_plugin  *p_plugin);

public Q_SLOTS:
    void selected_scene_items();
private :
    Polyhedron_demo_group_plugin  * m_plugin;
};
#endif
