#ifndef ARRANGEMENT_DEMO_WINDOW_H
#define ARRANGEMENT_DEMO_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_ArrangementDemoWindow.h"
#include "ArrangementGraphicsItem.h"
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
    void on_tabWidget_currentChanged( );

signals:
    void modelChanged( );

protected:
    void setupUi( );
    void resetCallbackState( int tabIndex );
    void removeCallback( int tabIndex );

    std::vector< ArrangementDemoTabBase* > tabs;
    std::vector< CGAL::Object > arrangements;
    std::vector< QAction* > activeModes; // for the current tab
    int lastTabIndex;

    Ui::ArrangementDemoWindow* ui;
    QActionGroup* modeGroup;
    QActionGroup* envelopeGroup;
    QActionGroup* snapGroup;
};
#endif // ARRANGEMENT_DEMO_WINDOW_H
