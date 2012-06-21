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
    
    ArrangementDemoWindow(QWidget* parent = 0);
    ~ArrangementDemoWindow();
    
public slots:
    void updateMode( QAction* a );
    void updateEnvelope( QAction* a );
    void updateSnapping( QAction* a );
    void on_actionQuit_triggered( );

signals:
    void modelChanged( );

protected:
    void setupUi( );

    ArrangementDemoTabBase* tab;

    CGAL::Qt::GraphicsItem* agi;
    CGAL::Qt::GraphicsViewSegmentInputBase* segmentInputCallback;
    CGAL::Qt::Callback* deleteCurveCallback;
    CGAL::Qt::Callback* pointLocationCallback;
    VerticalRayShootCallbackBase* verticalRayShootCallback;
    CGAL::Qt::Callback* mergeEdgeCallback;
    SplitEdgeCallbackBase* splitEdgeCallback;
    EnvelopeCallbackBase* envelopeCallback;
    Seg_arr arrangement;
    QGraphicsScene* scene;
    Ui::ArrangementDemoWindow* ui;
    QActionGroup* modeGroup;
    QActionGroup* envelopeGroup;
    QActionGroup* snapGroup;
    QAction* activeMode;
};
#endif // ARRANGEMENT_DEMO_WINDOW_H
