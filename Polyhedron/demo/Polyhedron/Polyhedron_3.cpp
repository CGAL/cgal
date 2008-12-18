#include "MainWindow.h"
#include <QApplication>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polyhedron_3 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2); 
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  if(!args.empty() && args[0] == "--use-meta")
  {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
    args.removeAt(0);
  }

  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }
  return app.exec();
}

#ifndef USE_FORWARD_DECL
#  include "Scene.cpp"
#  include "Scene_moc.cpp"
#  include "Scene_polyhedron_operations.cpp"
#  include "Scene_tex_polyhedron_operations.cpp"
#  include "Scene_nef_polyhedron_operations.cpp"
#  include "Scene_nef_and_polyhedron_operations.cpp"
#  include "Scene_rendering.cpp"
#  include "Scene_nef_rendering.cpp"
#  include "Scene_tex_rendering.cpp"
#  include "Viewer.cpp"
#  include "Viewer_moc.cpp"
#  include "MainWindow.cpp"
#  include "MainWindow_moc.cpp"
#  include "MainWindow_boolean_operations.cpp"
#  include "MainWindow_convex_hull.cpp"
#  include "MainWindow_curvature_estimation.cpp"
#  include "MainWindow_inside_out.cpp"
#  include "MainWindow_kernel.cpp"
#  include "MainWindow_pca.cpp"
#  include "MainWindow_remeshing.cpp"
#  include "MainWindow_self_intersection.cpp"
#  include "MainWindow_simplify.cpp"
#  include "MainWindow_subdivision_methods.cpp"
#  include "MainWindow_parameterization.cpp"
#  include "texture.cpp"
#endif
