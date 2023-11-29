// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/Qt/utility.h>

typedef QRectF (*mapToSceneFunction)(const QGraphicsView* , const QRect);

int main()
{
  mapToSceneFunction f = CGAL::Qt::mapToScene;
  return (&f > 0) ? 0 : 1;
}
