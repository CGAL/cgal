// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken, Algorithmic Solutions

#ifndef CGAL_WINDOW_X_WINDOW_H
#define CGAL_WINDOW_X_WINDOW_H


#include <CGAL/LEDA/X11/cursorfont.h>

namespace CGAL {

// windows, events, colors ...

enum { 
  key_press_event, 
  key_release_event, 
  button_press_event, 
  button_release_event,
  panel_press_event,
  configure_event, 
  exposure_event, 
  motion_event, 
  destroy_event, 
  timer_event,
  no_event 
};


extern const char* event_name[];


enum {
 KEY_BACKSPACE,
 KEY_RETURN,
 KEY_ESCAPE,
 KEY_LEFT,
 KEY_RIGHT,
 KEY_UP,
 KEY_DOWN,
 KEY_HOME,
 KEY_END,
 KEY_PAGEUP,
 KEY_PAGEDOWN,
 KEY_PRINT,
 KEY_F1,
 KEY_F2,
 KEY_F3,
 KEY_F4,
 KEY_F5,
 KEY_F6,
 KEY_F7,
 KEY_F8,
 KEY_F9,
 KEY_F10,
 KEY_F11,
 KEY_F12,
 KEY_PRINT1
};


enum {
  DEF_COLOR = -16,
  invisible = -1,
  white  =  0,
  black  =  1,
  red    =  2,
  green  =  3,
  blue   =  4,
  yellow =  5,
  violet =  6,
  orange =  7,
  cyan   =  8,
  brown  =  9,
  pink   = 10,
  green2 = 11,
  blue2  = 12,
  grey1  = 13,
  grey2  = 14,
  grey3  = 15,
  ivory  = 16 
};

enum point_style  { pixel_point, cross_point, plus_point, circle_point, 
                    disc_point, rect_point, box_point };

enum line_style   {solid, dashed, dotted, dashed_dotted};
enum text_mode    {transparent, opaque};
enum drawing_mode {src_mode, xor_mode, or_mode, and_mode, diff_mode};
enum grid_style   {invisible_grid, point_grid, line_grid};

}


#endif
