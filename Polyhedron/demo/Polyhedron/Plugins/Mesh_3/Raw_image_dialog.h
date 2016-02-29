#ifndef RAW_IMAGE_DIALOG_H
#define RAW_IMAGE_DIALOG_H

#include "ui_raw_image.h"

class Raw_image_dialog : public QDialog, public Ui::Raw_image_dialog 
{
  Q_OBJECT

public:
  Raw_image_dialog(QWidget* parent = 0);

  unsigned int image_word_size() const;

private Q_SLOTS:
  void update_image_size();
};

#endif // RAW_IMAGE_DIALOG
