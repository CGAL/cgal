#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <fstream>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QDebug>
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Messages_interface.h"


//This plugin crates an action in Operations that displays "Hello World" in the 'console' dockwidet.
class SurfaceMeshIoPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")
public:
  void init(QMainWindow*, CGAL::Three::Scene_interface*, Messages_interface* m)
  {
    this->message = m;
  }
   bool applicable(QAction*) const
   {
     return false;
   }
   QList<QAction*> actions() const
   {
     return QList<QAction*>();
   }
   QString name() const { return "surface_mesh_io_plugin"; }
   QString nameFilters() const { return "OFF files to Surface_mesh (*.off)"; }
   bool canLoad() const { return true; }
   CGAL::Three::Scene_item* load(QFileInfo fileinfo) {
     if(fileinfo.suffix().toLower() != "off") return 0;

     // Open file
     std::ifstream in(fileinfo.filePath().toUtf8());
     if(!in) {
      message->error(QString("Cannot open file %1").arg((const char*)fileinfo.filePath().toUtf8()));
       return NULL;
     }

     Scene_surface_mesh_item::SMesh *surface_mesh = new Scene_surface_mesh_item::SMesh();
     in >> *surface_mesh;
     if(!in || surface_mesh->is_empty())
     {
       delete surface_mesh;
       in.close();
       // Try to read .off in a polygon soup
       Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item();
       soup_item->setName(fileinfo.completeBaseName());
       std::ifstream in2(fileinfo.filePath().toUtf8());
       if(!soup_item->load(in2)) {
         message->error(QString("Cannot open file %1").arg((const char*)fileinfo.filePath().toUtf8()));
         delete soup_item;
         return 0;
       }
       message->information("The facets don't seem to be oriented. Loading a Soup of polygons instead."
                            "To convert it to a Surface_mesh or a Polyhedron, use Polygon Mesh Processing -> Orient polygon soup");
       return soup_item;
     }
     Scene_surface_mesh_item* item = new Scene_surface_mesh_item(surface_mesh);
     item->setName(fileinfo.completeBaseName());
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
   Messages_interface* message;
};
#include "Surface_mesh_io_plugin.moc"
