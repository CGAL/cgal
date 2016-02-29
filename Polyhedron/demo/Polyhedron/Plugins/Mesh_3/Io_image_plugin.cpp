#include "config.h"
#include "Scene_segmented_image_item.h"
#include "Image_type.h"
#include "ui_Image_res_dialog.h"

#include <QMessageBox>
#include <QSettings>
#include <QUrl>
#include "Raw_image_dialog.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <fstream>
class Io_image_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  Io_image_plugin() {}

  QString nameFilters() const;
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo) { return false; }
  QString name() const { return "segmented images"; }
};


QString Io_image_plugin::nameFilters() const {
  return QString("Inrimage files (*.inr *.inr.gz) ;; Analyze files (*.hdr *.img *img.gz)");
}


bool Io_image_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Io_image_plugin::load(QFileInfo fileinfo) {
  Image* image = new Image;
  QApplication::restoreOverrideCursor();
  if(!image->read(fileinfo.filePath().toUtf8()))
    {
      QMessageBox qmb(QMessageBox::NoIcon,
                      "Raw Dialog",
                      tr("Error with file <tt>%1</tt>:\n"
                         "unknown file format!\n"
                         "\n"
                         "Open it as a raw image?").arg(fileinfo.fileName()),
                      QMessageBox::Yes|QMessageBox::No);

      bool success = true;
      if(qmb.exec() == QMessageBox::Yes) {
        Raw_image_dialog raw_dialog;
        raw_dialog.label_file_size->setText(QString("%1 B").arg(fileinfo.size()));
        if( raw_dialog.exec() ){
          QApplication::setOverrideCursor(Qt::WaitCursor);

          if(image->read_raw(fileinfo.filePath().toUtf8(),
			     raw_dialog.dim_x->value(),
			     raw_dialog.dim_y->value(),
			     raw_dialog.dim_z->value(),
			     raw_dialog.spacing_x->value(),
			     raw_dialog.spacing_y->value(),
			     raw_dialog.spacing_z->value(),
                             raw_dialog.offset->value())){
            QSettings settings;
            settings.beginGroup(QUrl::toPercentEncoding(fileinfo.absoluteFilePath()));
            settings.setValue("is_raw", true);
            settings.setValue("dim_x", raw_dialog.dim_x->value());
            settings.setValue("dim_y", raw_dialog.dim_y->value());
            settings.setValue("dim_z", raw_dialog.dim_z->value());
            settings.setValue("spacing_x", raw_dialog.spacing_x->value());
            settings.setValue("spacing_y", raw_dialog.spacing_y->value());
            settings.setValue("spacing_z", raw_dialog.spacing_z->value());
            settings.setValue("offset", raw_dialog.offset->value());
            settings.endGroup();
          }else {
            success = false;
          }
        }else {
          success = false;
        }
      }else {
        success = false;   
      }
      if(!success){
        delete image;
        return NULL;
      }
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
  QApplication::restoreOverrideCursor();
  int return_code = dialog.exec();
  if(return_code != QDialog::Accepted)
  {
    delete image;
    return NULL;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
  // Get selected precision
  int voxel_scale = ui.precisionList->currentIndex() + 1;

  Scene_segmented_image_item* image_item = 
    new Scene_segmented_image_item(image,voxel_scale);
  image_item->setName(fileinfo.baseName());

  return image_item;
}

bool Io_image_plugin::canSave(const CGAL::Three::Scene_item*)
{
  return false;
}

#include "Io_image_plugin.moc"
