#include "mainwindow.h"

int main(int argc, char** argv)
{
  QApplication application(argc,argv);

  MainWindow w;

  if(argc>1)
    w.surface_open(argv[1]);

  w.show();

  return application.exec();
}
