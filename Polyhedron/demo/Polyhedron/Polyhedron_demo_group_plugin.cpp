#include "Polyhedron_demo_group_plugin.h"
#include "Polyhedron_demo_group_plugin.moc"
#include "Scene.h"

/****************
 * Group Plugin *
 ****************/

void Polyhedron_demo_group_plugin::init(QMainWindow* mainWindow,
          Scene_interface* scene_interface,Messages_interface* m ) {
    //get the references
    trueScene = dynamic_cast<Scene*>(scene_interface);
    this->scene = scene_interface;
    this->mw = mainWindow;
    messages = m;
    //creates and link the actions
    actionAddToGroup= new QAction("Add new group", mw);
    if(actionAddToGroup) {
        connect(actionAddToGroup, SIGNAL(triggered()),
                this, SLOT(add_group()));
    }
}
