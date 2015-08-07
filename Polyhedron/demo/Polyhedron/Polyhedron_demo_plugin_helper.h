#ifndef POLYHEDRON_DEMO_OPERATION_HELPER_H
#define POLYHEDRON_DEMO_OPERATION_HELPER_H

#include "Scene_item_config.h" //defines SCENE_ITEM_EXPORT

#include <QString>
#include <QStringList>
#include <QMap>

class QAction;
struct QMetaObject;
class QMainWindow;
class QDockWidget;

#include "Polyhedron_demo_plugin_interface.h"
#include "Scene_interface.h"

class SCENE_ITEM_EXPORT Polyhedron_demo_plugin_helper
  : public Polyhedron_demo_plugin_interface
{
public:
  // get action object from its name
  static QAction* getActionFromMainWindow(QMainWindow*, QString action_name);
  
  // Init plugin
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  
  // Get list of actions supported by this plugin
  virtual QStringList actionsNames() const;
  virtual QList<QAction*> actions() const;

  // To get a selected item with the type of SceneType
  template<class SceneType>
  SceneType* get_selected_item() const {
    int item_id = scene->mainSelectionIndex();
    SceneType* scene_item = qobject_cast<SceneType*>(scene->item(item_id));
    if(!scene_item) {
      // no selected SceneType - if there is only one in list return it, otherwise NULL
      int counter = 0;
      int last_selected = 0;
      for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end && counter < 2; ++i) {
        if(SceneType* tmp = qobject_cast<SceneType*>(scene->item(i))) { 
          scene_item = tmp;
          counter++; 
          last_selected=i;
        }
      }
      if(counter != 1) { return NULL; }
      scene->setSelectedItem(last_selected);
    }
    return scene_item;
  }

  void add_dock_widget(QDockWidget* dock);

  // Auto-connect actions to slots. Called by init().
  void autoConnectActions();
  
protected:
  QMap<QString, QAction*> actions_map;
  Scene_interface* scene;
  QMainWindow* mw;
};

#endif // POLYHEDRON_DEMO_OPERATION_HELPER_H
