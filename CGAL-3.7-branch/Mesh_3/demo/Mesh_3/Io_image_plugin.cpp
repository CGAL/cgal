#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
#  include <GL/glew.h>
#endif
#include "Scene_segmented_image_item.h"
#include "Image_type.h"

#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Io_image_plugin :
  public QObject,
  public Io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Io_plugin_interface)

public:
  Io_image_plugin() {
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
    glewInit();
#endif
  }

  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo) { return false; }
};

QStringList Io_image_plugin::nameFilters() const {
  return QStringList() << "Inrimage files (*.inr *.inr.gz)"
                       << "Analyze files (*.hdr *.img *img.gz)"
                       << "All files (*.*)";
}

bool Io_image_plugin::canLoad() const {
  return true;
}

Scene_item* 
Io_image_plugin::load(QFileInfo fileinfo) {
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

bool Io_image_plugin::canSave(const Scene_item*)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Io_image_plugin, Io_image_plugin)
#include "Io_image_plugin.moc"
