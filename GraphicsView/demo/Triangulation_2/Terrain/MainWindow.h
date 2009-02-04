#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>

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
  void open(const QString& fileName);
  void open_file();

  signals:
  void sceneChanged();
};




#endif
