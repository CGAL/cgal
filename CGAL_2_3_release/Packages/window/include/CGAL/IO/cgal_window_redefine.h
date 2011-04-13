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
// release       : 
// release_date  : 
// 
// file          : include/CGAL/IO/cgal_window_redefine.h
// package       : window (2.8.6)
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// revision      : 2.8.6
// revision_date : 27 June 2001
// author(s)     : Matthias Baesken
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#define leda_color      CGAL::color
#define leda_string     std::string
#define leda_panel      CGAL::panel
#define leda_window     CGAL::window

#define leda_invisible  CGAL::invisible
#define leda_white      CGAL::white
#define leda_black      CGAL::black
#define leda_red        CGAL::red
#define leda_green      CGAL::green
#define leda_blue       CGAL::blue
#define leda_yellow     CGAL::yellow
#define leda_violet     CGAL::violet
#define leda_orange     CGAL::orange
#define leda_cyan       CGAL::cyan
#define leda_brown      CGAL::brown
#define leda_pink       CGAL::pink
#define leda_green2     CGAL::green2
#define leda_blue2      CGAL::blue2
#define leda_grey1      CGAL::grey1
#define leda_grey2      CGAL::grey2
#define leda_grey3      CGAL::grey3
#define leda_ivory      CGAL::ivory


#define leda_line_style CGAL::line_style
#define leda_solid      CGAL::solid
#define leda_dashed     CGAL::dashed
#define leda_dotted     CGAL::dotted


#define leda_point_style   CGAL::point_style
#define leda_pixel_point   CGAL::pixel_point
#define leda_cross_point   CGAL::cross_point
#define leda_plus_point    CGAL::plus_point
#define leda_circle_point  CGAL::circle_point
#define leda_disc_point    CGAL::disc_point
#define leda_rect_point    CGAL::rect_point
#define leda_box_point     CGAL::box_point


#define leda_text_mode     CGAL::text_mode
#define leda_transparent   CGAL::transparent
#define leda_opaque        CGAL::opaque


#define leda_drawing_mode  CGAL::drawing_mode
#define leda_src_mode      CGAL::src_mode 
#define leda_xor_mode      CGAL::xor_mode
#define leda_or_mode       CGAL::or_mode
#define leda_and_mode      CGAL::and_mode

#define grid_style         CGAL::grid_style
#define invisible_grid     CGAL::invisible_grid
#define point_grid         CGAL::point_grid 
#define line_grid          CGAL::line_grid

// events
#define key_press_event      CGAL::key_press_event 
#define key_release_event    CGAL::key_release_event
#define button_press_event   CGAL::button_press_event
#define button_release_event CGAL::button_release_event
#define panel_press_event    CGAL::panel_press_event
#define configure_event      CGAL::configure_event
#define exposure_event       CGAL::exposure_event 
#define motion_event         CGAL::motion_event
#define destroy_event        CGAL::destroy_event  
#define timer_event          CGAL::timer_event
#define no_event             CGAL::no_event


#endif
