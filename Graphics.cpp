#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>

#include "Graphics.h"
//////////////
//	globals //
//////////////
Graphics* global_graphics = NULL;

/////////////////
//CTR's & DTR's//
/////////////////

Graphics::Graphics(int argc, char* argv[],double min_x, double max_x , double min_y , double max_y )
    :app(argc, argv), _min_x(min_x), _max_x(max_x), _min_y(min_y), _max_y(max_y)
{
  //const int SCREEN_SCENE_SIZE = 1000;
  int SCREEN_SCENE_SIZE = 1000;
  _scene.setSceneRect(-SCREEN_SCENE_SIZE,-SCREEN_SCENE_SIZE, 2*SCREEN_SCENE_SIZE, 2*SCREEN_SCENE_SIZE);

  //we assume that the absolute value in the workspace is in the square (-0.5-0.5) -> (0.5,0.5)
  //const int SCREEN_SIZE = 100;
  int FILE_SCENE_SIZE = 1;
  _file_scene.setSceneRect(-FILE_SCENE_SIZE ,-FILE_SCENE_SIZE ,2*FILE_SCENE_SIZE , 2*FILE_SCENE_SIZE );

  _ratio_x = SCREEN_SCENE_SIZE / (2.5*(_max_x - _min_x));
  _ratio_y = SCREEN_SCENE_SIZE / (2.5*(_max_y - _min_y));

  global_graphics = this;
}

Graphics::~Graphics(void)
{}

/////////////////
//DISPLAY SCENE//
/////////////////
void Graphics::display()
{
  QGraphicsView* view = new QGraphicsView(&_scene);
  CGAL::Qt::GraphicsViewNavigation navigation;
  view->installEventFilter(&navigation);
  view->viewport()->installEventFilter(&navigation);
  view->setRenderHint(QPainter::Antialiasing);
  view->show();

  app.exec();

  delete (view);
  return;
}
/////////////////
//EXPORT  SCENE//
/////////////////
void Graphics::save_scene(const std::string& path_name)
{
  static int count = 0;
  static int max_num_of_frames = 1000;
  //create filename
  std::stringstream str;
  //str << "../images/";
  str << path_name.c_str();
  str << "img";
  str << max_num_of_frames + count++;
  //str << ".png";
  str << ".jpg";
  
  QGraphicsView* view = new QGraphicsView(&this->_file_scene);
  CGAL::Qt::GraphicsViewNavigation navigation;
  view->installEventFilter(&navigation);
  view->viewport()->installEventFilter(&navigation);
  view->setRenderHint(QPainter::Antialiasing);
  
  // a white background
  _file_scene.setBackgroundBrush(Qt::white);
  
  QImage img(1024,768,QImage::Format_ARGB32_Premultiplied);
  QPainter p(&img);
  _file_scene.render(&p);
  p.end();
  
  img.save(str.str().c_str());
  delete (view);
 
  return;
}

/////////////////
//CLEAR SCENE//
/////////////////
void Graphics::clear()
{
  _scene.clear();
  _file_scene.clear();
  return;
}