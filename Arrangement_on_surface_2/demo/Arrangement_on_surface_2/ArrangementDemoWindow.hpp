#ifndef ARRANGEMENT_DEMO_WINDOW_HPP
#define ARRANGEMENT_DEMO_WINDOW_HPP
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementDemoWindow.h"
#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPointInput.h>
#include <CGAL/Qt/GraphicsViewSegmentInput.h>
#include <CGAL/IO/pixmaps/hand.xpm>
#include "ArrangementTypes.h"
#include "DeleteCurveCallback.hpp"

#include <Qt>

//#include <QFileDialog>
//#include <QInputDialog>
//#include <QMessageBox>
//#include <QtGui>

class ArrangementDemoWindow : public CGAL::Qt::DemosMainWindow
{
Q_OBJECT
public:
    typedef Seg_traits::Point_2 Point;
    typedef Seg_traits::Segment_2 Segment;
    
    ArrangementDemoWindow(QWidget* parent = 0);

    ~ArrangementDemoWindow();
    
public slots:
    void processInput( CGAL::Object o );
    void updateMode( QAction* a );
    void on_actionQuit_triggered( );

signals:
    void modelChanged( );

protected:
    void setupUi( );

    CGAL::Qt::ArrangementGraphicsItem< Seg_arr >* agi;
    CGAL::Qt::GraphicsViewSegmentInput< Seg_traits >* segmentInputCallback;
    DeleteCurveCallback< Seg_arr >* deleteCurveCallback;
    Seg_arr arrangement;
    QGraphicsScene scene;
    Ui::ArrangementDemoWindow* ui;
    QActionGroup* modeGroup;
    QAction* activeMode;
};
#endif // ARRANGEMENT_DEMO_WINDOW_HPP
