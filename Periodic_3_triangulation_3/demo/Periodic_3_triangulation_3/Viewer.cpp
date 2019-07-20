 #include "Viewer.h"
 #include <CGAL/Qt/CreateOpenGLContext.h>

Viewer::Viewer(QWidget *parent)
: CGAL::QGLViewer(parent)
{}
Viewer::~Viewer()
{}
