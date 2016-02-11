 #include "Viewer.h"
 #include <CGAL/Qt/CreateOpenGLContext.h>

Viewer::Viewer(QWidget *parent)
: QGLViewer(CGAL::Qt::createOpenGLContext(),parent)
{}
Viewer::~Viewer()
{}
