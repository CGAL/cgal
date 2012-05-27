#ifndef ARRANGEMENTS_DEMO_WINDOW_HPP
#define ARRANGEMENTS_DEMO_WINDOW_HPP
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementsDemoWindow.h"
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPointInput.h>
#include "ArrangementTypes.h"

//#include <QFileDialog>
//#include <QInputDialog>
//#include <QMessageBox>
//#include <QtGui>


class ArrangementsDemoWindow : public CGAL::Qt::DemosMainWindow,
    private Ui::ArrangementsDemoWindow
{
Q_OBJECT
public:
    typedef Seg_traits::Point_2 Point;
    
    ArrangementsDemoWindow(QWidget* parent = 0);

    ~ArrangementsDemoWindow();

public slots:
    void processInput( CGAL::Object o );
    void on_actionQuit_triggered( );

private:
    CGAL::Qt::ArrangementGraphicsItem< Seg_arr >* agi;
    CGAL::Qt::GraphicsViewPointInput< Seg_traits >* pointInputCallback;
    Seg_arr arrangement;
    QGraphicsScene scene;

};
#endif // ARRANGEMENTS_DEMO_WINDOW_HPP
