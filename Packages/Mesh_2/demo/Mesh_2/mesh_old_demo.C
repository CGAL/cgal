class Msh;
class Point;
class Line;

#include "viewer.h"
#include <iostream>
#include <qapplication.h>

int main(int argc, char **argv) {
  QApplication a(argc, argv);
  TrFrame *fr=new TrFrame();
  fr->setCaption("Triangulation Viewer");
  fr->show();
  a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
  return a.exec();
}
