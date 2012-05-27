#include "ArrangementsDemoWindow.hpp"

ArrangementsDemoWindow::
ArrangementsDemoWindow(QWidget* parent)
: CGAL::Qt::DemosMainWindow(parent)
{
    this->agi = new CGAL::Qt::ArrangementGraphicsItem< Seg_arr >( &( this->arrangement ) );
    this->pointInputCallback = new CGAL::Qt::GraphicsViewPointInput< Seg_traits >( this );

    this->setupUi( this );

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

ArrangementsDemoWindow::
~ArrangementsDemoWindow()
{

}

void
ArrangementsDemoWindow::
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
ArrangementsDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
