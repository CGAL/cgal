#ifndef MESH_3_DEMO_PLUGIN_INTERFACE_H
#define MESH_3_DEMO_PLUGIN_INTERFACE_H

#include <QString>
#include <QList>
#include <QtPlugin>

class QAction;
class QMainWindow;
class Scene_interface;
class Messages_interface;

class Plugin_interface 
{
public:
  virtual ~Plugin_interface() {}
  virtual void init(QMainWindow*, Scene_interface*) {};
  virtual void init(QMainWindow* mw, Scene_interface* sc, Messages_interface*) {
    init(mw, sc);
  };
  virtual QList<QAction*> actions() const = 0;
};

Q_DECLARE_INTERFACE(Plugin_interface,
                    "com.geometryfactory.Mesh3Demo.PluginInterface/1.0")

#endif // POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
