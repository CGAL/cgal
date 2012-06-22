#ifndef ARRANGEMENT_DEMO_WINDOW_H
#define ARRANGEMENT_DEMO_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementDemoWindow.h"
#include <CGAL/Qt/ArrangementGraphicsItem.h>
#include <CGAL/IO/pixmaps/hand.xpm>
#include "ArrangementTypes.h"
#include "ArrangementSegmentInputCallback.h"
#include "DeleteCurveCallback.h"
#include "PointLocationCallback.h"
#include "VerticalRayShootCallback.h"
#include "MergeEdgeCallback.h"
#include "SplitEdgeCallback.h"
#include "EnvelopeCallback.h"
#include "ArrangementDemoTab.h"

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
    typedef enum TraitsType {
        SEGMENT_TRAITS,
        POLYLINE_TRAITS,
        CONIC_TRAITS
    } TraitsType;
    
    ArrangementDemoWindow(QWidget* parent = 0);
    ~ArrangementDemoWindow();

    ArrangementDemoTabBase* makeTab( TraitsType tt );
    
public slots:
    void updateMode( QAction* a );
    void updateEnvelope( QAction* a );
    void updateSnapping( QAction* a );
    void on_actionNewTab_triggered( );
    void on_actionQuit_triggered( );

signals:
    void modelChanged( );

protected:
    void setupUi( );

    std::vector< ArrangementDemoTabBase* > tabs;
    std::vector< CGAL::Object > arrangements;

    Ui::ArrangementDemoWindow* ui;
    QActionGroup* modeGroup;
    QActionGroup* envelopeGroup;
    QActionGroup* snapGroup;
    QAction* activeMode;
};
#endif // ARRANGEMENT_DEMO_WINDOW_H
