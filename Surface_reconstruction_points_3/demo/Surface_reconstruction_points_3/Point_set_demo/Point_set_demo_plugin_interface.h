#ifndef POINT_SET_DEMO_PLUGIN_INTERFACE_H
#define POINT_SET_DEMO_PLUGIN_INTERFACE_H

#include <QString>
#include <QList>
#include <QtPlugin>

class QAction;
class QMainWindow;
class Scene_interface;
class Messages_interface;

class Point_set_demo_plugin_interface 
{
public:
  virtual void init(QMainWindow*, Scene_interface*) {};
  virtual void init(QMainWindow* mw, Scene_interface* sc, Messages_interface*) {
    init(mw, sc);
  };
  virtual QList<QAction*> actions() const = 0;
};

Q_DECLARE_INTERFACE(Point_set_demo_plugin_interface,
                    "org.cgal.PointSetDemo.PluginInterface/1.0")

#endif // POINT_SET_DEMO_PLUGIN_INTERFACE_H
