#ifndef POLYHEDRON_DEMO_GROUP_H
#define POLYHEDRON_DEMO_GROUP_H
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QList>
#include "Scene_group_item.h"
#include "CGAL_demo/Scene_interface.h"
#include "Scene.h"

class Scene_interface;
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
        return true;
    }
    QList<QAction*> actions() const {
        return QList<QAction*>() << actionAddToGroup;
    }
    Scene* trueScene;
public Q_SLOTS:
    void add_group() {
        scene->addItem(new Scene_group_item("new group"));
        trueScene->group_added();
    }


private:
    QAction* actionAddToGroup;
    void print_message(QString message) { messages->information(message); }
    Messages_interface* messages;

}; //end of class Polyhedron_demo_group_plugin
#endif
