#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include "grid.xpm"

class show_segments_layer : public CGAL::Qt_widget_layer {
public:
  bool show_hp;
  bool show_output;
  bool show_input;
  bool show_grid;

  /*! draw_grid - draw the grid
   */
  void draw_grid()
  {
    *widget << CGAL::WHITE << CGAL::LineWidth(1);

    // get the edge coordinate
    int min_x = static_cast<int>(widget->x_min());
    int min_y = static_cast<int>(widget->y_min());
    int max_x = static_cast<int>(widget->x_max());
    int max_y = static_cast<int>(widget->y_max());

    int i;
    for (i = min_x; i <= max_x; i++) {
      Segment_2 seg(Point_2(i, min_y), Point_2(i, max_y));
      *widget << seg;
    }
    for (i = min_y; i <= max_y; i++) {
      Segment_2 seg(Point_2(min_x, i), Point_2(max_x, i));
      *widget << seg;
    }
  }

  /*!
   */
  void draw()
  {
    if (show_grid) draw_grid();
    
    widget->lock();
    widget->setRasterOp(CopyROP);
    if(show_input) {
      *widget << CGAL::WHITE << CGAL::LineWidth(1);
      for(Segment_2_list_const_iterator i1 = seg_list.begin();
        i1 != seg_list.end();
        ++i1)
      *widget << *i1;
    }

    *widget << CGAL::LineWidth(2);
    for(Polyline_2_list_const_iterator i = output_list.begin();
        i != output_list.end();
        ++i) {
      Point_2_list_const_iterator prev = i->begin();
      Point_2_list_const_iterator i2 = prev;
      bool seg_painted = false;

      if(show_hp)
        *widget << CGAL::GREEN
                << Iso_rectangle_2(Point_2(i2->x() -
                                           prec / Number_type(2.0),
                                           i2->y() -
                                           prec / Number_type(2.0)),
                                   Point_2(i2->x() +
                                           prec / Number_type(2.0),
                                           i2->y() +
                                           prec / Number_type(2.0)));
      for(++i2; i2 != i->end(); ++i2) {
        seg_painted = true;
        if(show_output)
          *widget << CGAL::RED << Segment_2(*prev,*i2);
        if(show_hp)
          *widget << CGAL::GREEN << 
            Iso_rectangle_2(Point_2(i2->x() - prec / Number_type(2.0),
                                    i2->y() - prec / Number_type(2.0)),
                        Point_2(i2->x() + prec / Number_type(2.0),
                                i2->y() + prec / Number_type(2.0)));
        prev = i2;
      }

      if(!seg_painted && show_output) // segment entirely inside hot pixel
        *widget << CGAL::RED << *(i->begin());
    }
    widget->unlock();  
  }
};

/* XPM */
static char *show_hot_points_small_xpm[] = {
/* columns rows colors chars-per-pixel */
"16 16 3 1",
"  c black",
". c green",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXX",
"XX...XXXXXXXXXXX",
"XX.X.XXXXXXXXXXX",
"XX...XXXXXXXXXXX",
"XXXXXXXXXX...XXX",
"XXXXXXXXXX.X.XXX",
"XXXXXXXXXX...XXX",
"XXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXX",
"XXX...XXXXXXXXXX",
"XXX.X.XXX...XXXX",
"XXX...XXX.X.XXXX",
"XXXXXXXXX...XXXX",
"XXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXX"
};

/* XPM */
static char *show_hot_points_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c black",
". c green",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXX......XXXXXXXXXXXXXXXXXXXXXXX",
"XXX.XXXX.XXXXXXXXXXXXXXXXXXXXXXX",
"XXX.XXXX.XXXXXXXXXXXXXXXXXXXXXXX",
"XXX.XXXX.XXXXXXXXXX......XXXXXXX",
"XXX.XXXX.XXXXXXXXXX.XXXX.XXXXXXX",
"XXX......XXXXXXXXXX.XXXX.XXXXXXX",
"XXXXXXXXXXXXXXXXXXX.XXXX.XXXXXXX",
"XXXXXXXXXXXXXXXXXXX.XXXX.XXXXXXX",
"XXXXXXXXXXXXXXXXXXX......XXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXX......XXXXXXXXXXXXXXXXXXXXX",
"XXXXX.XXXX.XXXXXXXXXXXXXXXXXXXXX",
"XXXXX.XXXX.XXXXXXXXXXXXXXXXXXXXX",
"XXXXX.XXXX.XXXXXXXXX......XXXXXX",
"XXXXX.XXXX.XXXXXXXXX.XXXX.XXXXXX",
"XXXXX......XXXXXXXXX.XXXX.XXXXXX",
"XXXXXXXXXXXXXXXXXXXX.XXXX.XXXXXX",
"XXXXXXXXXXXXXXXXXXXX.XXXX.XXXXXX",
"XXXXXXXXXXXXXXXXXXXX......XXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static char *show_inputs_small_xpm[] = {
/* columns rows colors chars-per-pixel */
"16 16 3 1",
"  c black",
". c #C0C0C0",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXX",
"XXXXXXX .XXXXXXX",
"XXXX. XX .XXXX.X",
"X.XX. XX .XXX. X",
"X .. XXXX .X. XX",
"XX . XXXXX . XXX",
"XXX. XXXXX .XXXX",
"XX.  .XXX.  .XXX",
"XX. X .X. XX XXX",
"XX. XX . XXX .XX",
"XX. XX. .XXXX XX",
"X. XX. X .XXXXXX",
"X. X. XXX .XXXXX",
"X.XX XXXXX .XXXX",
"XXXXXXXXXXX XXXX",
"XXXXXXXXXXXXXXXX"
};

/* XPM */
static char *show_inputs_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c black",
". c #C0C0C0",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXX.XX",
"XXXXXXXX. XXXXXXXXXXXXXXXXXX. XX",
"XXXXXXXX. XXXXXXXXXXXXXXXXXX XXX",
"XXXXXXXX. XXXXXXXXXXXXXXXXX. XXX",
"XXXXXXX. XXXXX .XXXXXXXXXX. XXXX",
"XXXXXXX. XXXXXX .XXXXXXXX. XXXXX",
"XXXXXXX. XXXXXXX .XXXXXXX. XXXXX",
"XXXXXXX. XXXXXXXX .XXXXX. XXXXXX",
"XXXXXXX. XXXXXXXXX .XXX. XXXXXXX",
"XXXXXX. XXXXXXXXXXX .XX. XXXXXXX",
"XXXXXX. XXXXXXXXXXXX .. XXXXXXXX",
"XXXXXX. XXXXXXXXXXXXX .XXXXXXXXX",
"XXXXXX. XXXXXXXXXXXXX. .XXXXXXXX",
"XXXXXX. XXXXXXXXXXXX. X .XXXXXXX",
"XXXXX.X XXXXXXXXXXX. XXX .XXXXXX",
"XXXXX. XXXXXXXXXXX. XXXXX .XXXXX",
"XXXXX. XXXXXXXXXXX. XXXXXX .XXXX",
"XX..X. XXXXXXXXXX. XXXXXXXX .XXX",
"XX  ...XXXXXXXXX. XXXXXXXXXX .XX",
"XXXX.   .XXXXXXX XXXXXXXXXXXX XX",
"XXXX. XX  ..XXX. XXXXXXXXXXXXXXX",
"XXXX. XXXX  ... XXXXXXXXXXXXXXXX",
"XXXX. XXXXXX  ....XXXXXXXXXXXXXX",
"XXXX. XXXXXXX.  XX...XXXXXXXXXXX",
"XXX. XXXXXXX. XX   XX..XXXXXXXXX",
"XXX. XXXXXX. XXXXXX    ...XXXXXX",
"XXX. XXXXXX. XXXXXXXXXX   ..XXXX",
"XXXX XXXXX. XXXXXXXXXXXXXX  XXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static char *show_outputs_small_xpm[] = {
/* columns rows colors chars-per-pixel */
"16 16 3 1",
"  c black",
". c red",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXX",
"XX. XXXXXXXXXXXX",
"XX. XXXXXXXXX .X",
"XX. XXXXXXXX .XX",
"XX. XXXXXXX .XXX",
"XX. XXXXXX .XXXX",
"XX. XXXXX .XXXXX",
"XX. XXXX .XXXXXX",
"XX.     .    XXX",
"XX...........XXX",
"XX. X. XXXXXXXXX",
"XX. XX. XXXXXXXX",
"XX. XXX. XXXXXXX",
"XXXXXXXX. XXXXXX",
"XXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXX"
};

/* XPM */
static char *show_outputs_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c black",
". c red",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XX. XXXXXXXXXXXXXXXXXXXX. XXXXXX",
"XXX. XXXXXXXXXXXXXXXXXXX. XXXXXX",
"XXXX. XXXXXXXXXXXXXXXXXX. XXXXXX",
"XXXX.. XXXXXXXXXXXXXXXXX. XXXXXX",
"XXXX. . XXXXXXXXXXXXXXXX. XXXXXX",
"XXXX. X. XXXXXXXXXXXXXXX. XXXXXX",
"XXXX. XX. XXXXXXXXXXXXXX. XXXXXX",
"XXXX. XXX. XXXXXXXXXXXXX. XXXXXX",
"XXXX. XXXX.             . XXXXXX",
"XXXX. XXXXX...................XX",
"XXXX. XXXXXX. XXXXXXXXXXXXXXXXXX",
"XXXX. XXXXXXX. XXXXXXXXXXXXXXXXX",
"XXXX. XXXXXXXX. XXXXXXXXXXXXXXXX",
"XXXX. XXXXXXXXX. XXXXXXXXXXXXXXX",
"XXXX. XXXXXXXXXX. XXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXX. XXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX.. XXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. . XXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. X. XXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XX. XXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXX.XXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXXXXXXXXXXX",
"XX................. XXXXXXXXXXXX",
"XX                . XXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXX. XXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};
