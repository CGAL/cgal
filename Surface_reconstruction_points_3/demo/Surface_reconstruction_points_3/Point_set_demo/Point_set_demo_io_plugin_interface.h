#ifndef POINT_SET_DEMO_IO_PLUGIN_INTERFACE_H
#define POINT_SET_DEMO_IO_PLUGIN_INTERFACE_H

#include <QFileInfo>
#include <QStringList>

class Scene_item;

class Point_set_demo_io_plugin_interface 
{
public:
  virtual QStringList nameFilters() const = 0;

  virtual bool canLoad() const = 0;
  virtual Scene_item* load(QFileInfo fileinfo) = 0;

  virtual bool canSave(const Scene_item*) = 0;
  virtual bool save(const Scene_item*, QFileInfo fileinfo) = 0;
};

Q_DECLARE_INTERFACE(Point_set_demo_io_plugin_interface,
                    "org.cgal.PointSetDemo.IOPluginInterface/1.0")

#endif // POINT_SET_DEMO_IO_PLUGIN_INTERFACE_H
