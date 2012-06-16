#include "ArrangementDemoWindow.h"
#include <QApplication>

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );

    ArrangementDemoWindow demoWindow;
    demoWindow.show( );

    return app.exec( );
}
