#include <CGAL/Three/Three.h>


using namespace CGAL::Three;

QMainWindow* Three::s_mainwindow = NULL;
Scene_interface* Three::s_scene = NULL;
Three* Three::s_three = NULL;

QMainWindow* Three::mainWindow()
{
  return s_mainwindow;
}

Scene_interface* Three::scene()
{
  return s_scene;
}

Three* Three::messages()
{
  return s_three;
}

Three::Three()
{
  Three::s_three = this;
}

Three::~Three()
  {
    delete s_three;
  }
