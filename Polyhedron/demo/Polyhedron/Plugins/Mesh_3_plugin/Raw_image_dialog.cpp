#include "Raw_image_dialog.h"

Raw_image_dialog::Raw_image_dialog(QWidget* parent)
  : QDialog(parent) 
{
  setupUi(this);
}

void Raw_image_dialog::update_image_size() {
  label_image_size->setNum((int)image_word_size() *
			   dim_x->value() *
			   dim_y->value() *
			   dim_z->value());
}

unsigned int Raw_image_dialog::image_word_size() const {
  if(short_bt->isChecked())
    return 2;
  if(int_bt->isChecked())
    return 4;
  else if(float_bt->isChecked())
    return 4;
  else if(double_bt->isChecked())
    return 8;
  else
    return 1;
}

