#include <QtGui>
#include "window.h"

int main(int argv, char **args)
{	
    srand(1);
	QApplication app(argv, args);
	app.setApplicationName("Reconstruction_simplification_2 Demo");
	MainWindow window;
	window.show();
	return app.exec();
}
