// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : include/CGAL/LEDA/impl/x_window.h
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

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
