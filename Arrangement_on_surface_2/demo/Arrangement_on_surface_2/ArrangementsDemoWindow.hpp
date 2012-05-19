#ifndef ARRANGEMENTS_DEMO_WINDOW_HPP
#define ARRANGEMENTS_DEMO_WINDOW_HPP
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementsDemoWindow.h"
//#include <QFileDialog>
//#include <QInputDialog>
//#include <QMessageBox>
//#include <QtGui>


class ArrangementsDemoWindow : public CGAL::Qt::DemosMainWindow,
    private Ui::ArrangementsDemoWindow
{
Q_OBJECT
public:
  ArrangementsDemoWindow(QWidget* parent = 0);

  ~ArrangementsDemoWindow();

public slots:
    void on_actionQuit_triggered( );

};
#endif
