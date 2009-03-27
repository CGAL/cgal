#ifndef POLYHEDRON_DEMO_OPERATION_HELPER_H
#define POLYHEDRON_DEMO_OPERATION_HELPER_H

#include "Scene_item_config.h" //defines SCENE_ITEM_EXPORT

#include <QString>
#include <QStringList>
#include <QMap>

class QAction;
struct QMetaObject;
class QMainWindow;
class Scene_interface;

#include "Polyhedron_demo_plugin_interface.h"

class SCENE_ITEM_EXPORT Polyhedron_demo_plugin_helper
  : public Polyhedron_demo_plugin_interface
{
public:
  static QAction* getActionFromMainWindow(QMainWindow*, QString action_name);
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  virtual QStringList actionsNames() const;
  virtual QList<QAction*> actions() const;

  void autoConnectActions();
protected:
  QMap<QString, QAction*> actions_map;
  Scene_interface* scene;
  QMainWindow* mw;
};

#endif // POLYHEDRON_DEMO_OPERATION_HELPER_H
