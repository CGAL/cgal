// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_OFF_TO_NEF_3_H
#define CGAL_OFF_TO_NEF_3_H

#include <CGAL/Nef_polyhedron_3.h>

// --- begin preliminary number type converter -----------------

#ifndef CGAL_NUMBER_TYPE_CONVERTER_NEF_3_H
#define CGAL_NUMBER_TYPE_CONVERTER_NEF_3_H

#include<sstream>
#include<CGAL/Gmpz.h>
#include<CGAL/Gmpq.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Nef_S2/Normalizing.h>
#include<CGAL/Nef_nary_union_3.h>

#ifdef CGAL_USE_LEDA
#include<CGAL/leda_integer.h>
#include<CGAL/leda_rational.h>
#endif // CGAL_USE_LEDA

CGAL_BEGIN_NAMESPACE

class Homogeneous_tag;
class Cartesian_tag;
template<typename Tag, typename Kernel> class number_type_converter_nef_3;

template<class Kernel>
class number_type_converter_nef_3<Homogeneous_tag, Kernel> {
 public:
  static typename Kernel::Point_3
    convert (const CGAL::Cartesian<double>::Point_3 &d)
  {
    typedef typename Kernel::Point_3   Point_3;
    typedef typename Kernel::RT        RT;
    
    CGAL::Gmpq x(d.x()), y(d.y()), z(d.z());
    
      CGAL::Homogeneous<CGAL::Gmpz>::Point_3 b =
	normalized ( CGAL::Homogeneous<CGAL::Gmpz>::Point_3 (
	 x.numerator()   * y.denominator() * z.denominator(),
         x.denominator() * y.numerator()   * z.denominator(),
         x.denominator() * y.denominator() * z.numerator(),
         x.denominator() * y.denominator() * z.denominator() ) );

      std::ostringstream outx, outy, outz, outw;
      outx << b.hx();
      outy << b.hy();
      outz << b.hz();
      outw << b.hw();
      RT rx (outx.str().c_str());
      RT ry (outy.str().c_str());
      RT rz (outz.str().c_str());
      RT rw (outw.str().c_str());
      
      return Point_3 (rx, ry, rz, rw);
  }
};

template<class Kernel>
class number_type_converter_nef_3<Cartesian_tag, Kernel> {
 public:
  static typename Kernel::Point_3
    convert (const CGAL::Cartesian<double>::Point_3 &d)
  {  return typename Kernel::Point_3(d.x(), d.y(), d.z());
  }
};

CGAL_END_NAMESPACE
#endif // NUMBER_TYPE_CONVERTER_NEF_3_H

// --- end preliminary number type converter -------------------

#include <iostream>
#include <map>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/Nef_3/Mark_bounded_volumes.h>
#include <CGAL/IO/Scanner_OFF.h>
#include <CGAL/normal_vector_newell_3.h>

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
#include <CGAL/Timer.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Nef_3>
std::size_t
OFF_to_nef_3 (std::istream &i_st, Nef_3 &nef_union, bool verb=false)
{
   // Nef
   typedef typename Nef_3::Kernel          Kernel;
   typedef typename Nef_3::Point_3         Point_3;
   typedef typename std::vector<Point_3>   Point_set;
   CGAL::Nef_nary_union_3<Nef_3>            nary_union;
   // input data structure
   typedef double	                   Scan_NT;
   typedef CGAL::Cartesian<Scan_NT>        Scan_kernel;
   typedef Scan_kernel::Point_3            Scan_point;
   typedef std::vector<Scan_point>         Scan_point_set;
   typedef Scan_kernel::Vector_3           Scan_vector;

   typedef CGAL::Scanner_OFF<Scan_kernel>  Scan_OFF;
   typedef Scan_OFF::Vertex_iterator       Scan_vertex_it;
   typedef Scan_OFF::Facet_iterator        Scan_facet_it;
   typedef Scan_OFF::Index_iterator        Scan_index_it;

   typedef typename Kernel::Kernel_tag Kernel_tag;
   typedef typename CGAL::number_type_converter_nef_3<Kernel_tag,Kernel> ntc;

   // declarations and defaults
   std::size_t discarded_facets=0;
   long idx;

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
   CGAL::Timer t_convert, t_union;
#endif

   // input: description of polyhedron in object file format (OFF)
   // with cartesian double coordinates
   Scan_OFF scan (i_st);
   int NOV = scan.size_of_vertices();

   // read and store vertices
   Scan_point_set V_scan;
   V_scan.reserve (NOV);
   Point_set V;
   V.reserve (NOV);

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
   t_convert.start();
#endif

   Scan_vertex_it v_it = scan.vertices_begin();
   for (idx=0; v_it != scan.vertices_end(); ++v_it, ++idx)
   {  V_scan.push_back (*v_it);
     V.push_back (ntc::convert(*v_it));
   }
   CGAL_warning ( idx==NOV );
   NOV = idx;

   // for each facet
   Scan_facet_it f_it = scan.facets_begin();
   for (idx=0; f_it != scan.facets_end(); ++f_it, ++idx)
   {  // read facet
      Scan_facet_it::indices_size_type NOI=f_it.size_of_indices(), jdx;
      Scan_point_set V_f_scan;
      V_f_scan.reserve(NOI);
      Point_set V_f;
      V_f.reserve(NOI);

      Scan_index_it ind_it = f_it->begin();
      for (jdx=0; ind_it != f_it->end(); ++ind_it, ++jdx)
      {  // assertion: index out of range?
         CGAL_assertion ( 0 <= *ind_it && *ind_it < NOV );
         V_f_scan.push_back (V_scan[*ind_it]);
         V_f.push_back (V[*ind_it]);
      }
      CGAL_warning ( jdx==NOI );
      NOI = jdx;

      bool is_nef = false;
      CGAL_assertion_msg( V_f.size() >= 1 || !verb, "empty vertex cycle");
      if ( V_f.size() >= 1 )
      {  // compute Newell vector <double>
         Scan_vector normal;
         normal_vector_newell_3(V_f_scan.begin(),V_f_scan.end(),normal);

         // construct and enqueue Nef_polyhedron_3 <Kernel>
         Nef_3 nef (V_f.begin(), V_f.end(), normal, verb);
         if ( !nef.is_empty() )
	   {nary_union.add_polyhedron(nef);
            is_nef = true;
         }
      }

      if ( !is_nef )
      {  ++discarded_facets;
         if (verb)
	 {  std::cerr << "Hence, discard input facet " << (idx+1)
	       << " (enumerated beginning with 1)."
	       << " Check semantics!\n" << std::endl;
	 }
      }
   }

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
   t_convert.stop();
   std::cout << "time (conversion): " << t_convert.time()<< std::endl;
#endif

   // union of queue entries

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
   t_union.start();
#endif

#ifdef CGAL_NEF_OFF_TO_NEF_TIMER
   t_union.stop();
   std::cout << "time (union): " << t_union.time() << "\n" << std::endl;
#endif

   // return values
   //   if ( nef_map.size() == 0 ) nef_union = Nef_3 ();
   //   else
   nef_union = nary_union.get_union();
   CGAL::Mark_bounded_volumes<Nef_3> mbv (true);
   nef_union.delegate (mbv);

   return discarded_facets;
}

CGAL_END_NAMESPACE
#endif // CGAL_OFF_TO_NEF_3_H
