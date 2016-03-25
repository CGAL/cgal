#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <fstream>

#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QDebug>
#include "Scene_surface_mesh_item.h"


//This plugin crates an action in Operations that displays "Hello World" in the 'console' dockwidet.
class SurfaceMeshIoPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
 QString name() const { return "surface_mesh_io_plugin"; }
 QString nameFilters() const { return "OFF files (*.off)"; }
 bool canLoad() const { return true; }
 CGAL::Three::Scene_item* load(QFileInfo fileinfo) {
     if(fileinfo.suffix().toLower() != "off") return 0;

     // Open file
     std::ifstream in(fileinfo.filePath().toUtf8());
     if(!in) {
       std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
       return NULL;
     }

     Scene_surface_mesh_item::SMesh *surface_mesh = new Scene_surface_mesh_item::SMesh();
     in >> *surface_mesh;
     Scene_surface_mesh_item* item = new Scene_surface_mesh_item(surface_mesh);
     return item;

 }
 bool canSave(const CGAL::Three::Scene_item* ) {
     return false;
 }

 bool save(const CGAL::Three::Scene_item* , QFileInfo ) {
     return false;
 }



private:
  QList<QAction*> _actions;
  //The reference to the main window
  QMainWindow* mw;
};
#include "Surface_mesh_io_plugin.moc"
