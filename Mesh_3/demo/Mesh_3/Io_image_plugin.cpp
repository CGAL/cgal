#include "config.h"
#include "Scene_segmented_image_item.h"
#include "Image_type.h"
#include "ui_Image_res_dialog.h"

#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Io_image_plugin :
  public QObject,
  public Io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  Io_image_plugin() {
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  //  glewInit();
#endif
  }

  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo, QString) { return false; }
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
  
  // Get display precision
  QDialog dialog;
  Ui::ImagePrecisionDialog ui;
  ui.setupUi(&dialog);
  
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  
  // Add precision values to the dialog
  for ( int i=1 ; i<9 ; ++i )
  {
    QString s = tr("1:%1").arg(i*i*i);
    ui.precisionList->addItem(s);
  }
  
  // Open window
  int return_code = dialog.exec();
  if(return_code != QDialog::Accepted)
  {
    delete image;
    return NULL;
  }
  
  // Get selected precision
  int voxel_scale = ui.precisionList->currentIndex() + 1;

  Scene_segmented_image_item* image_item = 
    new Scene_segmented_image_item(image,voxel_scale);
  image_item->setName(fileinfo.baseName());

  return image_item;
}

bool Io_image_plugin::canSave(const Scene_item*)
{
  return false;
}

#include "Io_image_plugin.moc"
