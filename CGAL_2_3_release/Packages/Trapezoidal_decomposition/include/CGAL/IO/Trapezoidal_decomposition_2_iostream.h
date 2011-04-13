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
// release_date  : 1999, October 13
//
// file          : include/CGAL/Trapezoidal_decomposition_2/IO/Trapezoidal_decomposition_2_iostream.h
// package       : Trapezoidal decomposition 2 (1.04)
// author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator	 : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter		 : 
// ======================================================================
#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#include <CGAL/Trapezoidal_decomposition_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class Traits>
std::ostream& operator<<(
  std::ostream &out,
  const typename Trapezoidal_decomposition_2<Traits>::In_face_iterator& it)
{
	out << "In_face_iterator(";
	if (it.operator->()) out << *it; else out <<"end";
	return out << ")" << std::flush;
}

template < class Traits>
std::ostream& write(
   std::ostream &out,const Trapezoidal_decomposition_2<Traits>& td)
{
	return out << td;
}

template < class Traits>
std::ostream& operator<<(
   std::ostream &out,const Trapezoidal_decomposition_2<Traits>& td)
{
	return write(out,td.get_data_structure(),td.get_traits());
}

template < class Traits>
std::ostream& write(std::ostream &out,const Td_X_trapezoid<Traits>& t,
					const Traits& traits,bool validate=true)
{
	typedef Trapezoidal_decomposition_2<Traits> TD;
	typedef Td_X_trapezoid<Traits> X_trapezoid;
	typedef X_trapezoid* pointer;

	out << "(";
	if (!t.is_left_unbounded()) out << t.left(); else out << "-oo";
	out << ",";
	if (!t.is_right_unbounded()) out << t.right(); else out << "+oo";
	out << "," << std::flush;
	if (!t.is_bottom_unbounded()) out << t.bottom(); else out << "-oo";
	out << ",";
	if (!t.is_top_unbounded()) out << t.top(); else out << "+oo";
	out << ",neighbours(" << std::flush;
	
	// debug neighbours equivalence relation
	int max_size=4+1;
	int null_size=2,size=null_size,i,j;
	
	X_trapezoid* value[] =
	{
		0,
			(X_trapezoid*)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE,
			t.left_bottom_neighbour(),
			t.left_top_neighbour(),
			t.right_bottom_neighbour(),
			t.right_top_neighbour()
	};
	typedef char debug_string[256];
	debug_string name[]=
	{
		"none",
			"deleted",
			"lb",
			"lt",
			"rb",
			"rt"
	};
	for (j=null_size;j<=max_size;j++)
	{
		for (i=0;i<size;i++)
			if (value[j]==value[i])
			{
				CGAL_CLIB_STD::strcat(name[i],"=");
				CGAL_CLIB_STD::strcat(name[i],name[j]);
				break;
			}
            if (i==size)
			{
                value[size]=value[j];
                CGAL_CLIB_STD::strcpy(name[size++],name[j]);
			}
	}
	if (size==null_size) {out << "none";}
	for(j=null_size;j<size;j++)
	{
		out << name[j];
		// identify neighbours
		if (traits.is_degenerate_point(t) && value[j])
			out << "=" << value[j]->top();
		out << " ";
	}
	out << ")" << std::flush;
	
	if (t.is_active())
	{
		if (!traits.is_degenerate(t))
		{
            if (t.is_unbounded())
			{
                out << ",U";
                if (t.is_left_unbounded())
					out << (char) 174;
                if (t.is_bottom_unbounded())
					out << (char) 25;
                if (t.is_top_unbounded())
					out << (char) 24;
                if (t.is_right_unbounded())
					out << (char) 175;
			}
            else
				out << ",T";
		}
		else if (traits.is_degenerate_curve(t))
			out << ",C";
		else // if (t.is_degenerate_point())
			out << ",P";
	}
	else
		out << ",D";
	
/*	Calling t.is_valid requires the traits to be initialized with the proper 
bounding box */
	if (validate)
	{
		if (t.is_valid(&traits))
			out << ",+";
		else
			out << ",-";
	}
	
	out << ")" << std::flush;
	
	return out;
}

template < class Traits>
std::ostream& operator<<(std::ostream &out,const Td_X_trapezoid<Traits>& t)
{
	typedef Trapezoidal_decomposition_2<Traits> TD;
	typedef Td_X_trapezoid<Traits> X_trapezoid;
	typedef X_trapezoid* pointer;
	Traits traits;

	out << "(";
	if (!t.is_left_unbounded()) out << t.left(); else out << "-oo";
	out << ",";
	if (!t.is_right_unbounded()) out << t.right(); else out << "+oo";
	out << "," << std::flush;
	if (!t.is_bottom_unbounded()) out << t.bottom(); else out << "-oo";
	out << ",";
	if (!t.is_top_unbounded()) out << t.top(); else out << "+oo";
	out << ",neighbours(" << std::flush;
	
	// debug neighbours equivalence relation
	int max_size=4+1;
	int null_size=2,size=null_size,i,j;
	
	X_trapezoid* value[] =
	{
		0,
			(X_trapezoid*)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE,
			t.left_bottom_neighbour(),
			t.left_top_neighbour(),
			t.right_bottom_neighbour(),
			t.right_top_neighbour()
	};
	typedef char debug_string[256];
	debug_string name[]=
	{
		"none",
			"deleted",
			"lb",
			"lt",
			"rb",
			"rt"
	};
	for (j=null_size;j<=max_size;j++)
	{
		for (i=0;i<size;i++)
			if (value[j]==value[i])
			{
				CGAL_CLIB_STD::strcat(name[i],"=");
				CGAL_CLIB_STD::strcat(name[i],name[j]);
				break;
			}
            if (i==size)
			{
                value[size]=value[j];
                CGAL_CLIB_STD::strcpy(name[size++],name[j]);
			}
	}
	if (size==null_size) {out << "none";}
	for(j=null_size;j<size;j++)
	{
		out << name[j];
		// identify neighbours
		if (traits.is_degenerate_point(t) && value[j])
			out << "=" << value[j]->top();
		out << " ";
	}
	out << ")" << std::flush;
	
	if (t.is_active())
	{
		if (!traits.is_degenerate(t))
		{
            if (t.is_unbounded())
			{
                out << ",U";
                if (t.is_left_unbounded())
					out << (char) 174;
                if (t.is_bottom_unbounded())
					out << (char) 25;
                if (t.is_top_unbounded())
					out << (char) 24;
                if (t.is_right_unbounded())
					out << (char) 175;
			}
            else
				out << ",T";
		}
		else if (traits.is_degenerate_curve(t))
			out << ",C";
		else // if (t.is_degenerate_point())
			out << ",P";
	}
	else
		out << ",D";
	
/*	Calling t.is_valid requires the traits to be initialized with the proper 
bounding box
	if (traits)
	{
		if (t.is_valid(traits))
			out << ",+";
		else
			out << ",-";
	}
*/	
	out << ")" << std::flush;
	
	return out;
}

CGAL_END_NAMESPACE

#endif //CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H
