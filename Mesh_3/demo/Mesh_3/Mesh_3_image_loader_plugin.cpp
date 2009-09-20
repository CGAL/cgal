#include "Scene_segmented_image_item.h"
#include "Image_type.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

class Mesh_3_image_loader_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface);

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo) { return false; }
};

QStringList Mesh_3_image_loader_plugin::nameFilters() const {
  return QStringList() << "Inrimage files (*.inr *.inr.gz)"
                       << "Analize files (*.hdr *.img *img.gz)"
                       << "All files (*.*)";
};

bool Mesh_3_image_loader_plugin::canLoad() const {
  return true;
}

Scene_item* 
Mesh_3_image_loader_plugin::load(QFileInfo fileinfo) {
  Image* image = new Image;
  if(!image->read(fileinfo.filePath().toUtf8()))
  {
    delete image;
    return NULL;
  }

  Scene_segmented_image_item* image_item = 
    new Scene_segmented_image_item(image);
  image_item->setName(fileinfo.baseName());

  return image_item;
}

bool Mesh_3_image_loader_plugin::canSave(const Scene_item* item)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Mesh_3_image_loader_plugin, Mesh_3_image_loader_plugin);
#include "Mesh_3_image_loader_plugin.moc"
