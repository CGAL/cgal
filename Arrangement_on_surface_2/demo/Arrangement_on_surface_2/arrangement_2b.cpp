#include "ArrangementsDemoWindow.hpp"
//#include <CGAL/Qt/DemosMainWindow.h>
#include <QApplication>

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );

    ArrangementsDemoWindow demoWindow;
    demoWindow.show( );

    return app.exec( );
}
