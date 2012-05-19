#include "ArrangementsDemoWindow.hpp"

ArrangementsDemoWindow::
ArrangementsDemoWindow(QWidget* parent)
: CGAL::Qt::DemosMainWindow(parent)
{
    // construct a QGraphicsView
    this->setupUi( this );

    this->setupStatusBar( );
    this->addNavigation( this->graphicsView );
    this->setupOptionsMenu( );
    this->addAboutDemo( ":/help/about.html" );
    this->addAboutCGAL( );
}

ArrangementsDemoWindow::
~ArrangementsDemoWindow()
{

}

void 
ArrangementsDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}
