#include "Raw_image_dialog.h"
#include <QPushButton>

Raw_image_dialog::Raw_image_dialog(QWidget* parent)
  : QDialog(parent)
{
  setupUi(this);
}

void Raw_image_dialog::update_image_size() {
  int size = (int)image_word_size() *
      dim_x->value() *
      dim_y->value() *
      dim_z->value() +
      offset->value();
  label_image_size->setText(QString("%1 B").arg(size));

  if(label_image_size->text() == label_file_size->text())
    buttonBox->button(QDialogButtonBox::Open)->setEnabled(true);
  else
    buttonBox->button(QDialogButtonBox::Open)->setEnabled(false);
}

std::size_t Raw_image_dialog::image_word_size() const {
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

WORD_KIND Raw_image_dialog::image_word_kind() const
{
  if(float_bt->isChecked() ||
     double_bt->isChecked())
    return WK_FLOAT;
  else
    return WK_FIXED;
}
SIGN Raw_image_dialog::image_sign() const
{
  if(signed_bt->isChecked())
    return SGN_SIGNED;
  else
    return SGN_UNSIGNED;
}
