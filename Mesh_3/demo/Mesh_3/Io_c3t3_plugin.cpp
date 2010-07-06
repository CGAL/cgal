#include "Scene_c3t3_item.h"

#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Io_c3t3_plugin :
  public QObject,
  public Io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Io_plugin_interface)

public:
  virtual QStringList nameFilters() const;
  
  virtual bool canLoad() const { return false; }
  virtual Scene_item* load(QFileInfo) { return NULL; }

  virtual bool canSave(const Scene_item*);
  virtual bool save(const Scene_item*, QFileInfo fileinfo);
};


QStringList
Io_c3t3_plugin::nameFilters() const
{ 
  return QStringList() << "Mesh (*.mesh)";
}


bool
Io_c3t3_plugin::canSave(const Scene_item*)
{
  return true;
}

bool
Io_c3t3_plugin::save(const Scene_item* item, QFileInfo fileInfo)
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
Q_EXPORT_PLUGIN2(Io_c3t3_plugin, Io_c3t3_plugin)
#include "Io_c3t3_plugin.moc"
