// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel       <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_DEBUG_H
#define CGAL_SURFACE_SWEEP_2_DEBUG_H

#include <CGAL/license/Surface_sweep_2.h>

#include <CGAL/No_intersection_surface_sweep_2.h>

namespace CGAL {
namespace Surface_sweep_2 {

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//                         DEBUG UTILITIES                                //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_text(const char* text,
                                                      bool do_eol)
{
  if (m_need_indent)
    for (uint8_t i = m_indent_size; i != 0; --i) std::cout << " ";
  std::cout << text;
  m_need_indent = false;
  if (do_eol) print_eol();
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_eol()
{
  std::cout << std::endl;
  m_need_indent = true;
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::increase_indent()
{ m_indent_size += 2; }

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::decrease_indent()
{ m_indent_size -= 2; }

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_start(const char* name,
                                                       bool do_eol)
{
  print_text("Start ");
  print_text(name);
  if (do_eol) print_eol();
  increase_indent();
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_end(const char* name,
                                                     bool do_eol)
{
  decrease_indent();
  print_text("End ");
  print_text(name);
  if (do_eol) print_eol();
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_curve(const Subcurve* sc)
{
  if (m_need_indent)
    for (uint8_t i = m_indent_size; i != 0; --i) std::cout << " ";
  sc->Print();
  m_need_indent = false;
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::PrintEventQueue()
{
  print_text("Event queue: ", true);
  Event_queue_iterator iter = m_queue->begin();
  while (iter != m_queue->end()) {
    Event* e = *iter++;
    e->Print();
  }
  CGAL_SS_DEBUG(std::cout << "--------------------------------" << std::endl;)
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::PrintSubCurves()
{
  print_text("Sub curves: ", true);
  for (size_t i = 0; i < m_num_of_subCurves; ++i) m_subCurves[i].Print();
  print_eol();
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::PrintStatusLine()
{
  if (m_statusLine.size() == 0) {
    print_text("Status line: empty", true);
    return;
  }
  print_text("Status line: ");
  if (m_currentEvent->is_closed())
    std::cout << "(" << m_currentEvent->point() << ")";
  else {
    Arr_parameter_space x = m_currentEvent->parameter_space_in_x();
    Arr_parameter_space y = m_currentEvent->parameter_space_in_y();
    PrintOpenBoundaryType(x, y);
  }
  print_eol();
  increase_indent();
  Status_line_iterator iter = m_statusLine.begin();
  while (iter != m_statusLine.end()) {
    print_curve(*iter);
    print_eol();
    ++iter;
  }
  decrease_indent();
  print_text("Status line end");
  print_eol();
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::
PrintOpenBoundaryType(Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  switch (ps_x) {
   case ARR_LEFT_BOUNDARY:  std::cout << "left boundary"; return;
   case ARR_RIGHT_BOUNDARY: std::cout << "right boundary"; return;
   case ARR_INTERIOR:
   default: break;
  }

  switch (ps_y) {
   case ARR_BOTTOM_BOUNDARY: std::cout << "bottom boundary"; return;
   case ARR_TOP_BOUNDARY:    std::cout << "top boundary"; return;
   case ARR_INTERIOR:
   default: CGAL_error();
  }
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::PrintEvent(const Event* e)
{
  if (e->is_closed()) std::cout << e->point();
  else {
    Arr_parameter_space x = e->parameter_space_in_x();
    Arr_parameter_space y = e->parameter_space_in_y();
    PrintOpenBoundaryType(x, y);
    std::cout << " with open curve: " << e->curve();
  }
}

template <typename Vis>
void No_intersection_surface_sweep_2<Vis>::print_event_info(const Event* e)
{
  print_text("Event Info: ");
  PrintEvent(e);
  print_eol();
  increase_indent();
  print_text("Left curves:", true);
  increase_indent();
  Event_subcurve_const_iterator iter;
  for (iter = e->left_curves_begin(); iter != e->left_curves_end(); ++iter) {
    print_curve(*iter);
    print_eol();
  }
  decrease_indent();
  print_text("Right curves:", true);
  increase_indent();
  for (iter = e->right_curves_begin(); iter != e->right_curves_end(); ++iter) {
    print_curve(*iter);
    print_eol();
  }
  decrease_indent();
  decrease_indent();
  print_text("End Event Info", true);
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
