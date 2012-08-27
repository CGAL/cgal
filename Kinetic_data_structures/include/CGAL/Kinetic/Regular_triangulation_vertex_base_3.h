// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_KINETIC_REGULAR_VERTEX_BASE_3_H
#define CGAL_KINETIC_KINETIC_REGULAR_VERTEX_BASE_3_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

namespace CGAL { namespace Kinetic {
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits,
class Vertex_base= CGAL::Triangulation_vertex_base_3<typename SimulationTraits::Instantaneous_kernel> >
class Regular_triangulation_vertex_base_3:
public CGAL::Triangulation_vertex_base_with_info_3<typename SimulationTraits::Simulator::Event_key,
						   typename SimulationTraits::Instantaneous_kernel,
						   Vertex_base>
{
    private:
        typedef CGAL::Triangulation_vertex_base_with_info_3<typename SimulationTraits::Simulator::Event_key,
            typename SimulationTraits::Instantaneous_kernel,
            Vertex_base> P;
        typedef typename Vertex_base::Triangulation_data_structure   TDS;
    public:
        typedef TDS                            Triangulation_data_structure;
        typedef typename TDS::Cell_handle     Cell_handle;
        typedef typename TDS::Vertex_handle   Vertex_handle;
        typedef typename Vertex_base::Geom_traits Traits;

        typedef typename SimulationTraits::Simulator::Event_key Label;

        Regular_triangulation_vertex_base_3(): P() {
        }

        Regular_triangulation_vertex_base_3(Cell_handle f): P(f) {
        }

        template < typename TDS3 >
            struct Rebind_TDS
        {
            typedef typename Vertex_base::template Rebind_TDS<TDS3>::Other Cb3;
            typedef Regular_triangulation_vertex_base_3<SimulationTraits, Cb3>  Other;
        };
};

} } //namespace CGAL::Kinetic
#endif
