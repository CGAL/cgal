#ifndef THREE_H
#define THREE_H

#include <CGAL/license/Three.h>

#include <QString>
#include <QObject>
#include <CGAL/Three/Scene_interface.h>
#include <QMainWindow>

#ifdef three_EXPORTS
#  define THREE_EXPORT Q_DECL_EXPORT
#else
#  define THREE_EXPORT Q_DECL_IMPORT
#endif

namespace CGAL{
namespace Three{
class Polyhedron_demo_plugin_interface;
class THREE_EXPORT Three{
public:

  Three();
  virtual ~Three(){}
  static QMainWindow* mainWindow();
  static Scene_interface* scene();
  static QObject* connectableScene();
  static Three* messages();
  /*! \brief Adds a dock widget to the interface
   *
   * Adds a dock widget in the left section of the MainWindow. If the slot is already taken, the dock widgets will be tabified.
   */
  void addDockWidget(QDockWidget* dock_widget);

  /*! \brief Gets an item of the templated type.
   * \returns the first `SceneType` item found in the scene's list of currently selected items;
   * \returns NULL if there is no `SceneType` in the list.
   */
  template<class SceneType>
  static SceneType* getSelectedItem();

  /*! \brief Automatically connects each action of `plugin` to the corresponding slot.
   *
   * \attention Each action named `ActionName` in the plugin's `actions()` list must have a corresponding slot named `on_ActionsName_triggered()`
   * in the plugin.
   */
  static void autoConnectActions(CGAL::Three::Polyhedron_demo_plugin_interface* plugin);
protected:
  static QMainWindow* s_mainwindow;
  static Scene_interface* s_scene;
  static QObject* s_connectable_scene;
  static Three* s_three;

};
}
}

#endif // THREE_H
