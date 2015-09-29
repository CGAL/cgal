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
        return data_pool->group_entries().size();
    }
    QList<Scene_group_item*> groups() const
    {
        return data_pool->group_entries();
    }

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
    void init(QMainWindow* mainWindow,
              Scene_interface* scene_interface);

    bool applicable(QAction*) const
    {
        return true;
    }
    QList<QAction*> actions() const {
        return QList<QAction*>() << actionAddToGroup;
    }
    GroupModel *group_model;
    ItemsModel *items_model;
    void add_to_group() {
        Q_FOREACH(Scene_item* item, selected_items)
            selected_group->addChild(item);
    }
    QList<Scene_item*> selected_items;
    Scene_group_item* selected_group;
public Q_SLOTS:
    void GroupChoice();


private:
    QAction* actionAddToGroup;

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
