// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_2_DEBUG_H
#define CGAL_SWEEP_LINE_2_DEBUG_H

#include <CGAL/Basic_sweep_line_2.h>

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//                         DEBUG UTILITIES                                //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
print_text(const char* text, bool do_eol)
{
  if (m_need_indent)
    for (uint8_t i = m_indent_size; i != 0; --i) std::cout << " ";
  std::cout << text;
  m_need_indent = false;
  if (do_eol) print_eol();
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::print_eol()
{
  std::cout << std::endl;
  m_need_indent = true;
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::increase_indent()
{ m_indent_size += 2; }

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::decrease_indent()
{ m_indent_size -= 2; }

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
print_start(const char* name, bool do_eol)
{
  print_text("Start ");
  print_text(name);
  if (do_eol) print_eol();
  increase_indent();
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
print_end(const char* name, bool do_eol)
{
  decrease_indent();
  print_text("End ");
  print_text(name);
  if (do_eol) print_eol();
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
print_curve(const Base_subcurve* sc)
{
  if (m_need_indent)
    for (uint8_t i = m_indent_size; i != 0; --i) std::cout << " ";
  sc->Print();
  m_need_indent = false;
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::PrintEventQueue()
{
  CGAL_SL_DEBUG(std::cout << std::endl << "Event queue: " << std::endl;)
  Event_queue_iterator iter = m_queue->begin();
  while (iter != m_queue->end()) {
    Event *e = *iter;
     e->Print();
    ++iter;
  }
  CGAL_SL_DEBUG(std::cout << "--------------------------------" << std::endl;)
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::PrintSubCurves()
{
  CGAL_SL_DEBUG(std::cout << std::endl << "Sub curves: " << std::endl;)
  for (unsigned int i = 0 ; i < m_num_of_subCurves ; ++i)
    m_subCurves[i].Print();
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::PrintStatusLine()
{
  if ( m_statusLine.size() == 0) {
    std::cout << std::endl << "Status line: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Status line: (" ;
  if(m_currentEvent->is_closed())
    std::cout << m_currentEvent->point() << ")" << std::endl;
  else {
    Arr_parameter_space x = m_currentEvent->parameter_space_in_x(),
                  y = m_currentEvent->parameter_space_in_y();

    PrintOpenBoundaryType(x, y);
  }
  Status_line_iterator iter = m_statusLine.begin();
  while ( iter != m_statusLine.end()) {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Status line - end" << std::endl;
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
PrintOpenBoundaryType (Arr_parameter_space ps_x, Arr_parameter_space ps_y)
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

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
PrintEvent(const Event* e)
{
  if (e->is_closed()) std::cout << e->point();
  else {
    Arr_parameter_space x = e->parameter_space_in_x();
    Arr_parameter_space y = e->parameter_space_in_y();
    PrintOpenBoundaryType(x, y);
    std::cout << " with open curve: " << e->curve();
  }
}

template <typename Tr, typename Visit, typename Crv, typename Evnt,
          typename Alloc>
void Basic_sweep_line_2<Tr, Visit, Crv, Evnt, Alloc>::
print_event_info(const Event* e)
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

#endif
