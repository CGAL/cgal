// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Base_Node.h
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
//
// ======================================================================

#ifndef CGAL_BASE_NODE_H
#define CGAL_BASE_NODE_H
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/IO/PS_stream.h>

namespace CGAL {

	template < class Traits > // = Kd_tree_traits_point >
	class Base_node {
	public:
		typedef Traits::Item Item;
		typedef Base_node<Traits> Node;
		typedef typename Item::FT NT;
		typedef Traits::Separator Separator;

		Base_node() {};

                // is_leaf() should be defined for all derived classes
		// virtual const bool is_leaf() const =0;
                virtual const bool is_leaf() const {assert(false); return 0;}

                // methods below are not always defined for all classes
                // since it is not allowed to call them if they are
                // not defined, calling the virtual function below is erroneous
                // assert(false) is used to detect this type of error
                // an alternative would be to not define the virtual
                // functions below and to apply static_cast.
                // This does not work for generic classes
                // See the example at the start of Binary_search_tree.h

		virtual const int size() const {assert(false); return 0;}
		virtual Node* lower() const {assert(false); return 0;}
		virtual Node* upper() const {assert(false); return 0;}
		virtual Item** const begin() const {assert(false); return 0;}
		virtual Item** const end() const {assert(false); return 0;}
                virtual const Separator* separator() const {assert(false); return 0;}
                virtual const NT low_value() const {assert(false); return 0;}
                virtual const NT high_value() const {assert(false); return 0;}

           virtual void data_to_postscript(PS_Stream& PS,
            const int i, const int j,
            const NT mini, const NT maxi,
	    const NT minj, const NT maxj) {assert(false);};

        virtual ~Base_node() {};
   };


} // namespace CGAL
#endif // CGAL_BASE_NODE_H
