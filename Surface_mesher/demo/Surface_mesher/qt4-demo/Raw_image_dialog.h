#include "ui_raw_image.h"

struct Raw_image_dialog : public QDialog, public Ui::Raw_image_dialog {

  Raw_image_dialog(QWidget* parent = 0) : QDialog(parent) 
  {
    setupUi(this);
  }
};
