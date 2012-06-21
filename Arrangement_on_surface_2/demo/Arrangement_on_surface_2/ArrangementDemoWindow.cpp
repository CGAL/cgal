#include "ArrangementDemoWindow.h"
#include <QActionGroup>

ArrangementDemoWindow::
ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow( parent ),
    arrangement( Seg_arr( ) ),
    agi( new CGAL::Qt::ArrangementGraphicsItem< Seg_arr >( &( this->arrangement ) ) ),
    scene( new QGraphicsScene( -100, -100, 100, 100 ) ),
    ui( new Ui::ArrangementDemoWindow ),
    segmentInputCallback( new ArrangementSegmentInputCallback< Seg_arr >( &( this->arrangement ), this ) ),
    deleteCurveCallback( new DeleteCurveCallback< Seg_arr >( &( this->arrangement ), this ) ),
    pointLocationCallback( new PointLocationCallback< Seg_arr >( &( this->arrangement ), this ) ),
    verticalRayShootCallback( new VerticalRayShootCallback< Seg_arr >( &( this->arrangement ), this ) ),
    mergeEdgeCallback( new MergeEdgeCallback< Seg_arr >( &( this->arrangement ), this ) ),
    splitEdgeCallback( new SplitEdgeCallback< Seg_arr >( &( this->arrangement ), this ) ),
    envelopeCallback( new EnvelopeCallback< Seg_arr >( &( this->arrangement ), this ) ),
    tab( new ArrangementDemoTab< Seg_arr >( &( this->arrangement ), 0 ) )
{
    // set up the demo window
    this->setupUi( );
    this->setupStatusBar( );
    this->addNavigation( this->ui->graphicsView );
    this->setupOptionsMenu( );
    this->addAboutDemo( ":/help/about.html" );
    this->addAboutCGAL( );

    // set up the scene
    this->ui->graphicsView->setScene( this->scene );
    this->ui->graphicsView->setMouseTracking( true );
    this->scene->addItem( this->agi );

    // set up demo components
    this->segmentInputCallback->setScene( this->scene );
    this->deleteCurveCallback->setScene( this->scene );
    this->pointLocationCallback->setScene( this->scene );
    this->verticalRayShootCallback->setScene( this->scene );
    this->mergeEdgeCallback->setScene( this->scene );
    this->splitEdgeCallback->setScene( this->scene );
    this->envelopeCallback->setScene( this->scene );

    this->ui->tabWidget->addTab( tab, QString( ) );
    
    // set up callbacks
    this->scene->installEventFilter( this->segmentInputCallback );
    QObject::connect( this->modeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateMode( QAction* ) ) );
    QObject::connect( this->envelopeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateEnvelope( QAction* ) ) );
    QObject::connect( this->snapGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateSnapping( QAction* ) ) );
    QObject::connect( this->segmentInputCallback, SIGNAL( modelChanged( ) ), this, SIGNAL( modelChanged( ) ) );
    QObject::connect( this->deleteCurveCallback, SIGNAL( modelChanged( ) ), this, SIGNAL( modelChanged( ) ) );
    QObject::connect( this, SIGNAL( modelChanged( ) ), this->agi, SLOT( modelChanged( ) ) );
    QObject::connect( this, SIGNAL( modelChanged( ) ), this->envelopeCallback, SLOT( slotModelChanged( ) ) );
}

ArrangementDemoWindow::
~ArrangementDemoWindow( )
{
    //delete this->modeGroup;
}

void
ArrangementDemoWindow::
setupUi( )
{
    this->ui->setupUi( this );

    this->modeGroup = new QActionGroup( this );
    this->modeGroup->addAction( this->ui->actionDrag );
    this->modeGroup->addAction( this->ui->actionInsert );
    this->modeGroup->addAction( this->ui->actionDelete );
    this->modeGroup->addAction( this->ui->actionPointLocation );
    this->modeGroup->addAction( this->ui->actionRayShootingUp );
    this->modeGroup->addAction( this->ui->actionRayShootingDown );
    this->modeGroup->addAction( this->ui->actionMerge );
    this->modeGroup->addAction( this->ui->actionSplit );
    this->activeMode = this->ui->actionInsert;

    this->envelopeGroup = new QActionGroup( this );
    this->envelopeGroup->addAction( this->ui->actionLowerEnvelope );
    this->envelopeGroup->addAction( this->ui->actionUpperEnvelope );
    this->envelopeGroup->setExclusive( false );

    this->snapGroup = new QActionGroup( this );
    this->snapGroup->addAction( this->ui->actionSnapMode );
    this->snapGroup->addAction( this->ui->actionGridSnapMode );
    this->snapGroup->setExclusive( false );
    this->ui->actionGridSnapMode->setEnabled( false );
}

void
ArrangementDemoWindow::
updateMode( QAction* newMode )
{
    QWidget* widget = this->ui->tabWidget->currentWidget( );
    //ArrangementDemoTabBase* demoTab = static_cast< ArrangementDemoTabBase* >( widget );
    QGraphicsScene* activeScene = this->scene;

    // unhook the old active mode
    if ( this->activeMode == this->ui->actionInsert )
    {
        activeScene->removeEventFilter( this->segmentInputCallback );
    }
    else if ( this->activeMode == this->ui->actionDrag )
    {
        this->ui->graphicsView->setDragMode( QGraphicsView::NoDrag );
    }
    else if ( this->activeMode == this->ui->actionDelete )
    {
        this->deleteCurveCallback->reset( );
        activeScene->removeEventFilter( this->deleteCurveCallback );
    }
    else if ( this->activeMode == this->ui->actionPointLocation )
    {
        this->pointLocationCallback->reset( );
        activeScene->removeEventFilter( this->pointLocationCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingUp )
    {
        this->verticalRayShootCallback->reset( );
        activeScene->removeEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingDown )
    {
        this->verticalRayShootCallback->reset( );
        activeScene->removeEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionMerge )
    {
        this->mergeEdgeCallback->reset( );
        activeScene->removeEventFilter( this->mergeEdgeCallback );
    }
    else if ( this->activeMode == this->ui->actionSplit )
    {
        this->splitEdgeCallback->reset( );
        activeScene->removeEventFilter( this->splitEdgeCallback );
    }

    // update the active mode
    this->activeMode = newMode;

    // hook up the new active mode
    if ( this->activeMode == this->ui->actionInsert )
    {
        activeScene->installEventFilter( this->segmentInputCallback );
    }
    else if ( this->activeMode == this->ui->actionDrag )
    {
        this->ui->graphicsView->setDragMode( QGraphicsView::ScrollHandDrag );
    }
    else if ( this->activeMode == this->ui->actionDelete )
    {
        activeScene->installEventFilter( this->deleteCurveCallback );
    }
    else if ( this->activeMode == this->ui->actionPointLocation )
    {
        activeScene->installEventFilter( this->pointLocationCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingUp )
    {
        // -y is up for Qt, so we shoot down
        this->verticalRayShootCallback->setShootingUp( false );
        activeScene->installEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingDown )
    {
        // the bottom of the viewport for Qt is +y, so we shoot up
        this->verticalRayShootCallback->setShootingUp( true );
        activeScene->installEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionMerge )
    {
        activeScene->installEventFilter( this->mergeEdgeCallback );
    }
    else if ( this->activeMode == this->ui->actionSplit )
    {
        activeScene->installEventFilter( this->splitEdgeCallback );
    }
}

void
ArrangementDemoWindow::
updateEnvelope( QAction* newMode )
{
    bool show = newMode->isChecked( );
    if ( newMode == this->ui->actionLowerEnvelope )
    {
        this->envelopeCallback->showUpperEnvelope( show );
    }
    else if ( newMode == this->ui->actionUpperEnvelope )
    {
        this->envelopeCallback->showLowerEnvelope( show );
    }
}

void
ArrangementDemoWindow::
updateSnapping( QAction* newMode )
{
    bool enabled = newMode->isChecked( );
    if ( newMode == this->ui->actionSnapMode )
    {
        this->segmentInputCallback->setSnappingEnabled( enabled );
        this->splitEdgeCallback->setSnappingEnabled( enabled );
        if ( ! enabled )
        {
            this->ui->actionGridSnapMode->setChecked( false );
            this->ui->actionGridSnapMode->setEnabled( false );
            this->segmentInputCallback->setSnapToGridEnabled( false );
            this->splitEdgeCallback->setSnapToGridEnabled( false );
        }
        else
        {
            this->ui->actionGridSnapMode->setEnabled( true );
        }
    }
    else if ( newMode == this->ui->actionGridSnapMode )
    {
        this->segmentInputCallback->setSnapToGridEnabled( enabled );
        this->splitEdgeCallback->setSnapToGridEnabled( enabled );
        this->ui->graphicsView->setShowGrid( enabled );
    }
    this->scene->update( );
}

void 
ArrangementDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
