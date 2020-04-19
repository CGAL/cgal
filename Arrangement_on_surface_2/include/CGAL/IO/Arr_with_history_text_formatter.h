// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>

#ifndef CGAL_ARR_WITH_HISTORY_TEXT_FORMATTER_H
#define CGAL_ARR_WITH_HISTORY_TEXT_FORMATTER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The header file for the text-formatter classes.
 */

#include <CGAL/basic.h>
#include <iostream>

namespace CGAL {

/*! \class
 * A class defining a textual (ASCII) input/output format for arrangements
 * with history and supports reading and writing an arrangement from or to
 * input/output streams.
 */
template <class ArrFormatter_>
class Arr_with_history_text_formatter : public ArrFormatter_
{
public:

  typedef ArrFormatter_                                   Base;
  typedef Arr_with_history_text_formatter<Base>           Self;

  typedef typename Base::Arrangement_2                    Arr_with_history_2;
  typedef typename Arr_with_history_2::Size                             Size;
  typedef typename Arr_with_history_2::Dcel                             Dcel;
  typedef typename Arr_with_history_2::Curve_2            Curve_2;
  typedef typename Arr_with_history_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Arr_with_history_2::Point_2            Point_2;

  typedef typename Arr_with_history_2::Vertex_handle      Vertex_handle;
  typedef typename Arr_with_history_2::Halfedge_handle    Halfedge_handle;
  typedef typename Arr_with_history_2::Face_handle        Face_handle;

  typedef typename Arr_with_history_2::Vertex_const_handle
                                                      Vertex_const_handle;
  typedef typename Arr_with_history_2::Halfedge_const_handle
                                                      Halfedge_const_handle;
  typedef typename Arr_with_history_2::Face_const_handle
                                                      Face_const_handle;

  /*! Default constructor.*/
  Arr_with_history_text_formatter ():
    Base ()
  {}

  /*! Construct an output formatter. */
  Arr_with_history_text_formatter (std::ostream& os) :
    Base (os)
  {}

  /*! Construct an input formatter. */
  Arr_with_history_text_formatter (std::istream& is) :
    Base (is)
  {}

  /// \name Functions for writing curves.
  //@{

  /*! Write a begin-curves comment. */
  void write_curves_begin ()
  {
    __write_comment ("BEGIN CURVES");
    return;
  }

  /*! Write an end-curves comment. */
  void write_curves_end ()
  {
    __write_comment ("END CURVES");
    return;
  }

  /*! Write a specific curve. */
  void write_curve_begin ()
  {}

  void write_curve_end ()
  {}

  void write_curve (const Curve_2& c)
  {
    this->out() << c << std::endl;
    return;
  }

  void write_induced_edges_begin ()
  {}

  void write_induced_edges_end ()
  {
    this->out() << std::endl;
  }
  //@}

  /// \name Functions for reading curves.
  //@{

  /*! Start reading the curves. */
  void read_curves_begin ()
  {
    __skip_comments();
  }

  /*! Read the end-curves message. */
  void read_curves_end()
  {
    __skip_comments();
  }

  /*! Read a specific curve. */
  void read_curve_begin ()
  {}

  void read_curve_end ()
  {}

  void read_curve (Curve_2& c)
  {
    this->in() >> c;
    __skip_until_EOL();

    return;
  }

  void read_induced_edges_begin ()
  {}

  void read_induced_edges_end ()
  {
     __skip_until_EOL();
  }
  //@}

private:

  /*! Write a comment line. */
  void __write_comment (const char *str)
  {
    this->out() << "# " << str << std::endl;
    return;
  }

  /*! Skip until end of line. */
  void __skip_until_EOL ()
  {
    int     c;
    while ((c = this->in().get()) != EOF && c != '\n') {};
    return;
  }

  /*! Skip comment lines. */
  void __skip_comments ()
  {
    int     c;
    while ((c = this->in().get()) != EOF && c == '#')
      __skip_until_EOL();
    this->in().putback (c);

    return;
  }

};

} //namespace CGAL

#endif
