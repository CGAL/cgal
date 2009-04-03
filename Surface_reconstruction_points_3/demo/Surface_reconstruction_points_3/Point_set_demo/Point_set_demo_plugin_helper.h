#ifndef POINT_SET_DEMO_OPERATION_HELPER_H
#define POINT_SET_DEMO_OPERATION_HELPER_H

#include "Scene_item_config.h" //defines SCENE_ITEM_EXPORT

#include <QString>
#include <QStringList>
#include <QMap>

class QAction;
struct QMetaObject;
class QMainWindow;
class Scene_interface;

#include "Point_set_demo_plugin_interface.h"

class SCENE_ITEM_EXPORT Point_set_demo_plugin_helper
  : public Point_set_demo_plugin_interface
{
public:
  // get action object from its name
  static QAction* getActionFromMainWindow(QMainWindow*, QString action_name);
  
  // Init plugin
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  
  // Get list of actions supported by this plugin
  virtual QStringList actionsNames() const;
  virtual QList<QAction*> actions() const;

  // Auto-connect actions to slots. Called by init().
  void autoConnectActions();
  
protected:
  QMap<QString, QAction*> actions_map;
  Scene_interface* scene;
  QMainWindow* mw;
};

#endif // POINT_SET_DEMO_OPERATION_HELPER_H
