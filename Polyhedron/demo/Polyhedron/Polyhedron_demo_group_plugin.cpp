#include "Polyhedron_demo_group_plugin.h"
#include "Polyhedron_demo_group_plugin.moc"
#include "Scene.h"
#include <QMenu>
/****************
 * Group Plugin *
 ****************/

void Polyhedron_demo_group_plugin::init(QMainWindow* mainWindow,
          CGAL::Three::Scene_interface* scene_interface,Messages_interface* m ) {
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

    QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
    if ( NULL != menuFile )
    {
      QList<QAction*> menuFileActions = menuFile->actions();

      // Look for action just after "Load..." action
      QAction* actionAfterLoad = NULL;
      for ( QList<QAction*>::iterator it_action = menuFileActions.begin(),
           end = menuFileActions.end() ; it_action != end ; ++ it_action ) //Q_FOREACH( QAction* action, menuFileActions)
      {
        if ( NULL != *it_action && (*it_action)->text().contains("Load") )
        {
          ++it_action;
          if ( it_action != end && NULL != *it_action )
          {
            actionAfterLoad = *it_action;
          }
        }
      }

      // Insert "Load implicit function" action
      if ( NULL != actionAfterLoad )
      {
        menuFile->insertAction(actionAfterLoad,actionAddToGroup);
      }
    }
}
