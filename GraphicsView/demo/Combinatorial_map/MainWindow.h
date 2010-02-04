#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>
#include <QSlider>
#include <QFileDialog>
class QWidget;

class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

  public:
  MainWindow(QWidget* parent = 0);

  void connectActions();


  Scene scene;
  Timer timer;

public slots:
  void import_off();
  void add_off();
  void load_off(const QString& fileName);

  void import_3DTDS();
  void load_3DTDS(const QString& fileName);
  
  void subdivide();
  void create_cube();

 signals:
  void sceneChanged();

 private:
  unsigned int nbcube;
};




#endif
