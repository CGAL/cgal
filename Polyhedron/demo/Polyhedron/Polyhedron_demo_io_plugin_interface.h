#ifndef POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H

#include <QFileInfo>
#include <QStringList>

class Scene_item;

class Polyhedron_demo_io_plugin_interface 
{
public:
  virtual QString name() const = 0;
  virtual ~Polyhedron_demo_io_plugin_interface() {}
  
  virtual QString nameFilters() const = 0;

  virtual bool canLoad() const = 0;
  virtual Scene_item* load(QFileInfo fileinfo) = 0;

  virtual bool canSave(const Scene_item*) = 0;
  virtual bool save(const Scene_item*, QFileInfo fileinfo) = 0;
};

Q_DECLARE_INTERFACE(Polyhedron_demo_io_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

#endif // POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
