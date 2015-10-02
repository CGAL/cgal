   //! \file Polyhedron_demo_plugin_interface.h 
#ifndef POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_PLUGIN_INTERFACE_H

#include <QString>
#include <QList>
#include <QtPlugin>
#include <QDebug>

class QAction;
class QMainWindow;
class Messages_interface;
namespace CGAL {
namespace Three {
class Scene_interface;
  /*!
   * This class gives some virtual functions to help making a plugin
   */
class Polyhedron_demo_plugin_interface 
{
public:
  //! Destructor
  virtual ~Polyhedron_demo_plugin_interface() {}
  //!Initializes the plugin.
  virtual void init(QMainWindow*, CGAL::Three::Scene_interface*) {}
  //!Initializes the plugin.
  virtual void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface*) {
    init(mw, sc);
  }

  //! Checks the current state of the `Scene` or `MainWindow` and decides
  //! if the plugin can function, given that state.  Its actions are
  //! visible in contextmenus, if this returns true, not visible
  //! otherwise.  
  //!
  //! @returns \c true, if the plugin is applicable, \c false
  //! otherwise
  virtual bool applicable(QAction*) const = 0;
  //!The list of the actions of the plugin.
  virtual QList<QAction*> actions() const = 0;
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Polyhedron_demo_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

#endif // POLYHEDRON_DEMO_PLUGIN_INTERFACE_H
