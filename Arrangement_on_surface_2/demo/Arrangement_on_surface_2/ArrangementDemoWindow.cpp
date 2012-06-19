#include "ArrangementDemoWindow.h"
#include <QActionGroup>

ArrangementDemoWindow::
ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow( parent ),
    arrangement( Seg_arr( ) ),
    agi( new CGAL::Qt::ArrangementGraphicsItem< Seg_arr >( &( this->arrangement ) ) ),
    ui( new Ui::ArrangementDemoWindow ),
    segmentInputCallback( new ArrangementSegmentInputCallback< Seg_arr >( &( this->arrangement ), this ) ),
    deleteCurveCallback( new DeleteCurveCallback< Seg_arr >( &( this->arrangement ), this ) ),
    pointLocationCallback( new PointLocationCallback< Seg_arr >( &( this->arrangement ), this ) ),
    verticalRayShootCallback( new VerticalRayShootCallback< Seg_arr >( &( this->arrangement ), this ) ),
    mergeEdgeCallback( new MergeEdgeCallback< Seg_arr >( &( this->arrangement ), this ) ),
    splitEdgeCallback( new SplitEdgeCallback< Seg_arr >( &( this->arrangement ), this ) ),
    envelopeCallback( new EnvelopeCallback< Seg_arr >( &( this->arrangement ), this ) )
{
    // set up the demo window
    this->setupUi( );
    this->setupStatusBar( );
    this->addNavigation( this->ui->graphicsView );
    this->setupOptionsMenu( );
    this->addAboutDemo( ":/help/about.html" );
    this->addAboutCGAL( );

    // set up demo components
    this->segmentInputCallback->setScene( &( this->scene ) );
    this->deleteCurveCallback->setScene( &( this->scene ) );
    this->pointLocationCallback->setScene( &( this->scene ) );
    this->verticalRayShootCallback->setScene( &( this->scene ) );
    this->mergeEdgeCallback->setScene( &( this->scene ) );
    this->splitEdgeCallback->setScene( &( this->scene ) );
    this->envelopeCallback->setScene( &( this->scene ) );

    // set up the scene
    this->scene.setSceneRect( -100, -100, 100, 100 );
    this->ui->graphicsView->setScene( &( this->scene ) );
    this->ui->graphicsView->setMouseTracking( true );
    this->scene.addItem( this->agi );
    
    // set up callbacks
    this->scene.installEventFilter(  this->segmentInputCallback );
    QObject::connect( this->modeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateMode( QAction* ) ) );
    QObject::connect( this->envelopeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateEnvelope( QAction* ) ) );
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
}

void
ArrangementDemoWindow::
updateMode( QAction* newMode )
{
    // unhook the old active mode
    if ( this->activeMode == this->ui->actionInsert )
    {
        this->scene.removeEventFilter( this->segmentInputCallback );
    }
    else if ( this->activeMode == this->ui->actionDrag )
    {
        this->ui->graphicsView->setDragMode( QGraphicsView::NoDrag );
    }
    else if ( this->activeMode == this->ui->actionDelete )
    {
        this->deleteCurveCallback->reset( );
        this->scene.removeEventFilter( this->deleteCurveCallback );
    }
    else if ( this->activeMode == this->ui->actionPointLocation )
    {
        this->pointLocationCallback->reset( );
        this->scene.removeEventFilter( this->pointLocationCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingUp )
    {
        this->verticalRayShootCallback->reset( );
        this->scene.removeEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingDown )
    {
        this->verticalRayShootCallback->reset( );
        this->scene.removeEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionMerge )
    {
        this->mergeEdgeCallback->reset( );
        this->scene.removeEventFilter( this->mergeEdgeCallback );
    }
    else if ( this->activeMode == this->ui->actionSplit )
    {
        this->splitEdgeCallback->reset( );
        this->scene.removeEventFilter( this->splitEdgeCallback );
    }

    // update the active mode
    this->activeMode = newMode;

    // hook up the new active mode
    if ( this->activeMode == this->ui->actionInsert )
    {
        this->scene.installEventFilter( this->segmentInputCallback );
    }
    else if ( this->activeMode == this->ui->actionDrag )
    {
        this->ui->graphicsView->setDragMode( QGraphicsView::ScrollHandDrag );
    }
    else if ( this->activeMode == this->ui->actionDelete )
    {
        this->scene.installEventFilter( this->deleteCurveCallback );
    }
    else if ( this->activeMode == this->ui->actionPointLocation )
    {
        this->scene.installEventFilter( this->pointLocationCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingUp )
    {
        // -y is up for Qt, so we shoot down
        this->verticalRayShootCallback->setShootingUp( false );
        this->scene.installEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionRayShootingDown )
    {
        // the bottom of the viewport for Qt is +y, so we shoot up
        this->verticalRayShootCallback->setShootingUp( true );
        this->scene.installEventFilter( this->verticalRayShootCallback );
    }
    else if ( this->activeMode == this->ui->actionMerge )
    {
        this->scene.installEventFilter( this->mergeEdgeCallback );
    }
    else if ( this->activeMode == this->ui->actionSplit )
    {
        this->scene.installEventFilter( this->splitEdgeCallback );
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
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
