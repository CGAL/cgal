#ifndef CGAL_IMAGE_INTERFACE_H
#define CGAL_IMAGE_INTERFACE_H
#include <QDialog>
#include <QWidget>
#include <QApplication>

#include "ui_ImageInterface.h"
class ImageInterface: public QDialog, public Ui::ImageInterface
{
  Q_OBJECT
  qreal ratio;
  QWidget *currentlyFocused;
public:
  ImageInterface(QWidget *parent, qreal ratio)
    : QDialog(parent), ratio(ratio)
  {
    currentlyFocused = NULL;
    setupUi(this);
    connect(imgHeight, SIGNAL(valueChanged(int)),
            this, SLOT(imgHeightValueChanged(int)));

    connect(imgWidth, SIGNAL(valueChanged(int)),
            this, SLOT(imgWidthValueChanged(int)));

    connect(qApp, SIGNAL(focusChanged(QWidget*, QWidget*)),
            this, SLOT(onFocusChanged(QWidget*, QWidget*)));
  }
private Q_SLOTS:
  void imgHeightValueChanged(int i)
  {
    if(currentlyFocused == imgHeight
       && ratioCheckBox->isChecked())
    {imgWidth->setValue(i*ratio);}
  }

  void imgWidthValueChanged(int i)
  {
    if(currentlyFocused == imgWidth
       && ratioCheckBox->isChecked())
    {imgHeight->setValue(i/ratio);}
  }

  void onFocusChanged(QWidget*, QWidget* now)
  {
    currentlyFocused = now;
  }
};
#endif // CGAL_IMAGE_INTERFACE_H
