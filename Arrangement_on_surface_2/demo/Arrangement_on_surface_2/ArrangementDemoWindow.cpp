#include "ArrangementDemoWindow.hpp"
#include <QActionGroup>

ArrangementDemoWindow::
ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow( parent ),
    arrangement( Seg_arr( ) ),
    agi( new CGAL::Qt::ArrangementGraphicsItem< Seg_arr >( &( this->arrangement ) ) ),
    ui( new Ui::ArrangementDemoWindow ),
    segmentInputCallback( new CGAL::Qt::GraphicsViewSegmentInput< Seg_traits >( this ) ),
    deleteCurveCallback( new DeleteCurveCallback< Seg_arr >( &( this->arrangement ), this ) )
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

    // set up the scene
    this->scene.setSceneRect( -100, -100, 100, 100 );
    this->ui->graphicsView->setScene( &( this->scene ) );
    this->ui->graphicsView->setMouseTracking( true );
    this->scene.addItem( this->agi );
    
    // set up callbacks
    this->scene.installEventFilter( this->segmentInputCallback );
    QObject::connect( this->modeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateMode( QAction* ) ) );
    QObject::connect( this->segmentInputCallback, SIGNAL( generate( CGAL::Object ) ),
        this, SLOT( processInput( CGAL::Object ) ) );
    QObject::connect( this, SIGNAL( modelChanged( ) ), this->agi, SLOT( modelChanged( ) ) );
}

ArrangementDemoWindow::
~ArrangementDemoWindow( )
{
    delete this->modeGroup;
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
}

void
ArrangementDemoWindow::
processInput( CGAL::Object o )
{
    Segment segment;
    if ( CGAL::assign( segment, o ) )
    {
        // insert a segment
        Point p1 = segment.source( );
        Point p2 = segment.target( );
        Arr_xseg_2 curve( p1, p2 );
        CGAL::insert( this->arrangement, curve );
    }

    emit modelChanged( );
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
}

void 
ArrangementDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
