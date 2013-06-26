#ifndef POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_PLUGIN_INTERFACE_H

#include <QString>
#include <QList>
#include <QtPlugin>

class QAction;
class QMainWindow;
class Scene_interface;
class Messages_interface;

class Polyhedron_demo_plugin_interface 
{
public:
  virtual ~Polyhedron_demo_plugin_interface() {}
  virtual void init(QMainWindow*, Scene_interface*) {}
  virtual void init(QMainWindow* mw, Scene_interface* sc, Messages_interface*) {
    init(mw, sc);
  }

  //! Checks the current state of the `Scene` or `MainWindow` and decides
  //! if the plugin can function, given that state.  It's actions are
  //! visible in contextmenus, if this returns true, not visible
  //! otherwise.  
  //!
  //! @returns \c true, if the plugin is applicable, \c false
  //! otherwise
  virtual bool applicable() const = 0;
  virtual QList<QAction*> actions() const = 0;
};

Q_DECLARE_INTERFACE(Polyhedron_demo_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

#endif // POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
