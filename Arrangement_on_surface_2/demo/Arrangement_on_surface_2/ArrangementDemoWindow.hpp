#ifndef ARRANGEMENT_DEMO_WINDOW_HPP
#define ARRANGEMENT_DEMO_WINDOW_HPP
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementDemoWindow.h"
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPointInput.h>
#include <CGAL/IO/pixmaps/hand.xpm>
#include "ArrangementTypes.h"

//#include <QFileDialog>
//#include <QInputDialog>
//#include <QMessageBox>
//#include <QtGui>


class ArrangementDemoWindow : public CGAL::Qt::DemosMainWindow,
    private Ui::ArrangementDemoWindow
{
Q_OBJECT
public:
    typedef Seg_traits::Point_2 Point;
    
    ArrangementDemoWindow(QWidget* parent = 0);

    ~ArrangementDemoWindow();
    
public slots:
    void processInput( CGAL::Object o );
    void on_actionQuit_triggered( );

protected:
    void setup( );

    CGAL::Qt::ArrangementGraphicsItem< Seg_arr >* agi;
    CGAL::Qt::GraphicsViewPointInput< Seg_traits >* pointInputCallback;
    Seg_arr arrangement;
    QGraphicsScene scene;
};
#endif // ARRANGEMENT_DEMO_WINDOW_HPP
