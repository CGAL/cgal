// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef CGAL_CCB_CURVE_ITERATOR_H
#define CGAL_CCB_CURVE_ITERATOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


namespace CGAL {

 /*! \class
    * A circulator adapter that iterates over all curves in a CCB.
    */
    template <class Arrangement_>
   class Ccb_curve_iterator
   {
   public:

     typedef Arrangement_                            Arrangement;
     typedef typename Arrangement::Geometry_traits_2 Traits;
     typedef Ccb_curve_iterator<Arrangement>         Self;
     typedef typename Arrangement::Ccb_halfedge_const_circulator
                                                     Ccb_halfedge_const_circulator;
     typedef typename Traits::X_monotone_curve_2     value_type;
     typedef std::forward_iterator_tag               iterator_category;
     typedef const value_type&                       reference;
     typedef const value_type*                       pointer;
     typedef std::size_t                             difference_type;

   private:

     Ccb_halfedge_const_circulator    _first;    // The first halfedge.
     Ccb_halfedge_const_circulator    _circ;     // The current circulator.
     bool                       _done;     // Indicates whether we completed
                                           // a full traversal of the CCB.

   public:

     /*! Default constructor. */
     Ccb_curve_iterator () :
       _done (true)
     {}

     /*!
      * Constructor from a circulator.
      * \param circ A circulator for the first halfedge in the CCB.
      * \param done (true) in order to create a past-the-end iterator.
      */
     Ccb_curve_iterator (Ccb_halfedge_const_circulator circ,
                         bool done = false) :
       _first (circ),
       _circ (circ),
       _done (done)
     {}

     /*! Dereference operators. */
     reference operator* () const
     {
       return (_circ->curve());
     }

     pointer operator-> () const
     {
       return (&(_circ->curve()));
     }

     /*! Equality operators.*/
     bool operator== (const Self& it) const
     {
       return (_done == it._done && _circ == it._circ);
     }

     bool operator!= (const Self& it) const
     {
       return (_done != it._done || _circ != it._circ);
     }

     /*! Increment operators. */
     Self& operator++ ()
     {
       CGAL_assertion(!_done);
      
      --_circ;

      if (_circ == _first)
        _done = true;
      

       return (*this);
     }

     Self operator++ (int )
     {
       Self   temp = *this;
       ++(*this);
       return (temp);
     }

   };


} //namespace CGAL

#endif
