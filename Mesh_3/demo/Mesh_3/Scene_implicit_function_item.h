#ifndef SCENE_IMPLICIT_FUNCTION_ITEM_H
#define SCENE_IMPLICIT_FUNCTION_ITEM_H

#include <CGAL_demo/Scene_item_with_display_list.h>
#include <CGAL_demo/Scene_interface.h>
#include "Scene_implicit_function_item_config.h"
#include "implicit_functions/Implicit_function_interface.h"
#include "Color_ramp.h"

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#define SCENE_IMPLICIT_GRID_SIZE 120


class SCENE_IMPLICIT_FUNCTION_ITEM_EXPORT Scene_implicit_function_item 
  : public Scene_item_with_display_list
{
  Q_OBJECT
  
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  
public:
  Scene_implicit_function_item(Implicit_function_interface*);
  virtual ~Scene_implicit_function_item();
  
  Implicit_function_interface* function() const { return function_; }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  Bbox bbox() const;

  Scene_implicit_function_item* clone() const { return NULL; }

  // rendering mode
  virtual bool supportsRenderingMode(RenderingMode m) const;
  virtual bool manipulatable() const { return true; }
  virtual ManipulatedFrame* manipulatedFrame() { return frame_; }
  
  // draw (overload only direct_draw() to use display list of base class)
  virtual void direct_draw() const;
  
  virtual QString toolTip() const;

public slots:
  void compute_function_grid();

private:
  typedef qglviewer::Vec                  Point;
  typedef std::pair <Point,double>        Point_value;
  
  void draw_bbox() const;
  void draw_function_grid(const Color_ramp&, const Color_ramp&) const;
  void draw_grid_vertex(const Point_value&,
                        const Color_ramp&, const Color_ramp&) const;
  
  void compute_min_max();
  
private:
  Implicit_function_interface* function_;
  ManipulatedFrame* frame_;
  
  int grid_size_;
  double max_value_;
  double min_value_;
  Point_value implicit_grid_[SCENE_IMPLICIT_GRID_SIZE][SCENE_IMPLICIT_GRID_SIZE];
  
  Color_ramp blue_color_ramp_;
  Color_ramp red_color_ramp_;
};

#endif // SCENE_IMPLICIT_FUNCTION_ITEM
