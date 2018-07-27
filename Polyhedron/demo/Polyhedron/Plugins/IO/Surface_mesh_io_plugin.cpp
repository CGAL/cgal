#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <fstream>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QApplication>
#include <QDebug>
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Messages_interface.h"


class SurfaceMeshIoPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(
  CGAL::Three::Polyhedron_demo_plugin_interface
  CGAL::Three::Polyhedron_demo_io_plugin_interface
)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")
public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const 
  { 
    if(qobject_cast<const Scene_surface_mesh_item*>(item))
      return true; 
    return false;
  }
  
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
   QString loadNameFilters() const { return "OFF files to Surface_mesh (*.off);;Wavefront Surface_mesh OBJ (*.obj)"; }
   QString saveNameFilters() const { return "OFF files (*.off);;Wavefront OBJ (*.obj)"; }
   QString nameFilters() const { return QString(); }
   bool canLoad() const { return true; }
   CGAL::Three::Scene_item* load(QFileInfo fileinfo) {
     if(fileinfo.suffix().toLower() == "off")
     {
     // Open file
     std::ifstream in(fileinfo.filePath().toUtf8());
     if(!in) {
      message->error(QString("Cannot open file %1").arg((const char*)fileinfo.filePath().toUtf8()));
       return NULL;
     }

     SMesh *surface_mesh = new SMesh();
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
     std::size_t isolated_v = 0;
     BOOST_FOREACH(vertex_descriptor v, vertices(*surface_mesh))
     {
       if(surface_mesh->is_isolated(v))
       {
         ++isolated_v;
       }
     }
     if(isolated_v >0)
     {
       item->setNbIsolatedvertices(isolated_v);
       //needs two restore, it's not a typo
       QApplication::restoreOverrideCursor();
       QMessageBox::warning((QWidget*)NULL,
                      tr("Isolated vertices"),
                      tr("%1 isolated vertices found")
                      .arg(item->getNbIsolatedvertices()));
     }
     return item;
     }
     else if(fileinfo.suffix().toLower() == "obj")
     {
       // Open file
       std::ifstream in(fileinfo.filePath().toUtf8());
       if(!in) {
         std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
         return NULL;
       }
       Scene_surface_mesh_item* item = new Scene_surface_mesh_item();
       if(item->load_obj(in))
         return item;
     }

     return 0;

   }
   bool canSave(const CGAL::Three::Scene_item* item) {
       return qobject_cast<const Scene_surface_mesh_item*>(item) != 0;
   }

   bool save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo) {

     const Scene_surface_mesh_item* sm_item =
       qobject_cast<const Scene_surface_mesh_item*>(item);

     if(!sm_item)
       return false;

     std::ofstream out(fileinfo.filePath().toUtf8());
     out.precision (std::numeric_limits<double>::digits10 + 2);

     if(fileinfo.suffix().toLower() == "off"){
       return (sm_item && sm_item->save(out));
     }
     if(fileinfo.suffix().toLower() == "obj"){
       return (sm_item && sm_item->save_obj(out));
     }
     return false;
   }

private:
   QList<QAction*> _actions;
   //The reference to the main window
   Messages_interface* message;
};
#include "Surface_mesh_io_plugin.moc"
