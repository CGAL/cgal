#ifndef POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
#define POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H

#include <QFileInfo>
#include <QStringList>

class Scene_item;
namespace CGAL{
namespace Three {
class Polyhedron_demo_io_plugin_interface 
{
public:
  //!Returns the name of the plugin
  virtual QString name() const = 0;
  virtual ~Polyhedron_demo_io_plugin_interface() {}
  /*! The filters for the names of the files that can be used
   * by the plugin.
   * Example : to filter OFF files : return "OFF files (*.off)"
*/
  virtual QString nameFilters() const = 0;
  //! Specifies if the io_plugin is able to load an item or not.
  virtual bool canLoad() const = 0;
  //!  Loads an item from a file.
  virtual Scene_item* load(QFileInfo fileinfo) = 0;
  //!Specifies if the io_plugin can save the item or not.
  virtual bool canSave(const Scene_item*) = 0;
  //!Saves the item in the file corresponding to the path
  //!contained in fileinfo. Returns false if error.
  virtual bool save(const Scene_item*, QFileInfo fileinfo) = 0;
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Polyhedron_demo_io_plugin_interface,
                    "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

#endif // POLYHEDRON_DEMO_IO_PLUGIN_INTERFACE_H
