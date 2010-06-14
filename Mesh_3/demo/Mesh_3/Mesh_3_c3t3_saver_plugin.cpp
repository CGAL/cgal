#include "Scene_c3t3_item.h"

#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Mesh_3_c3t3_saver_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface);

public:
  virtual QStringList nameFilters() const;
  
  virtual bool canLoad() const { return false; }
  virtual Scene_item* load(QFileInfo fileinfo) { return NULL; }

  virtual bool canSave(const Scene_item*);
  virtual bool save(const Scene_item*, QFileInfo fileinfo);
};


QStringList
Mesh_3_c3t3_saver_plugin::nameFilters() const
{ 
  return QStringList() << "Mesh (*.mesh)";
}


bool
Mesh_3_c3t3_saver_plugin::canSave(const Scene_item* item)
{
  return true;
}

bool
Mesh_3_c3t3_saver_plugin::save(const Scene_item* item, QFileInfo fileInfo)
{
  const Scene_c3t3_item* c3t3_item = qobject_cast<const Scene_c3t3_item*>(item);
  if ( NULL == c3t3_item )
  {
    return false;
  }
  
  QString path = fileInfo.absoluteFilePath();
  std::ofstream medit_file (qPrintable(path));
  c3t3_item->c3t3().output_to_medit(medit_file,true,true);
  
  return true;
}


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Mesh_3_c3t3_saver_plugin, Mesh_3_c3t3_saver_plugin);
#include "Mesh_3_c3t3_saver_plugin.moc"
