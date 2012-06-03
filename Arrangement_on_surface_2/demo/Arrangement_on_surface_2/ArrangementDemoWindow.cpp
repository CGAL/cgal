#include "ArrangementDemoWindow.hpp"
#include <QActionGroup>

ArrangementDemoWindow::
ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow(parent)
//    demoMode( new QActionGroup( this ) )
{
    this->setupUi( this );
    QActionGroup* actionGroup = new QActionGroup( this );
    actionGroup->addAction( this->actionDrag );
    actionGroup->addAction( this->actionInsert );
    actionGroup->addAction( this->actionDelete );
    actionGroup->addAction( this->actionPointLocation );
    actionGroup->addAction( this->actionRayShootingUp );
    actionGroup->addAction( this->actionRayShootingDown );
    actionGroup->addAction( this->actionMerge );
    actionGroup->addAction( this->actionSplit );

    this->agi = new CGAL::Qt::ArrangementGraphicsItem< Seg_arr >( &( this->arrangement ) );
    this->pointInputCallback = new CGAL::Qt::GraphicsViewPointInput< Seg_traits >( this );

    // set up the scene
    this->scene.setSceneRect( -100, -100, 100, 100 );
    this->graphicsView->setScene( &( this->scene ) );
    this->graphicsView->setMouseTracking( true );

    this->setupStatusBar( );
    this->addNavigation( this->graphicsView );
    this->setupOptionsMenu( );
    this->addAboutDemo( ":/help/about.html" );
    this->addAboutCGAL( );
    
    // set up callbacks
    this->scene.installEventFilter( this->pointInputCallback );
    this->connect( this->pointInputCallback, SIGNAL( generate( CGAL::Object ) ),
        this, SLOT( processInput( CGAL::Object ) ) );
}

ArrangementDemoWindow::
~ArrangementDemoWindow( )
{

}


void
ArrangementDemoWindow::
setup( )
{

    /*
    this->demoMode = new QActionGroup( o );

    this->demoMode->addAction( this->actionDrag );
    this->demoMode->addAction( this->actionInsert );
    this->demoMode->addAction( this->actionDelete );
    this->demoMode->addAction( this->actionPointLocation );
    this->demoMode->addAction( this->actionRayShootingUp );
    this->demoMode->addAction( this->actionRayShootingDown );
    this->demoMode->addAction( this->actionMerge );
    this->demoMode->addAction( this->actionSplit );
    */
}

void
ArrangementDemoWindow::
processInput( CGAL::Object o )
{
    Point point;
    if ( CGAL::assign( point, o ) )
    {
        std::cout << "point generated (" 
            << point.x( ) << "," << point.y( ) << ")"
            << "... but we actually need segments to add to the arrangement!" 
            << std::endl;
    }
}

void 
ArrangementDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
