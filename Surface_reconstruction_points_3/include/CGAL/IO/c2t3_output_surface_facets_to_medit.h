// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Steve Oudot, Laurent Rineau, Nader Salman

#ifndef CGAL_C2T3_OUTPUT_SURFACE_FACETS_TO_MEDIT_H
#define CGAL_C2T3_OUTPUT_SURFACE_FACETS_TO_MEDIT_H

//#include <CGAL/Complex_2_in_triangulation_3.h>

#include <iomanip>
#include <stack>

namespace CGAL { namespace Surface_mesher {
    template < class Tr>
    typename Tr::size_type number_of_facets_on_surface(const Tr& T);
}
}

namespace Surface_mesher_io {

    template <class C2t3>
    void
        output_surface_facets_to_medit (std::ostream& os, const C2t3& c2t3)
    {
        using CGAL::Surface_mesher::number_of_facets_on_surface;

        typedef typename C2t3::Triangulation Tr;
        typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
        typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename Tr::Facet Facet;
        typedef typename Tr::Edge Edge;
        typedef typename Tr::Vertex_handle Vertex_handle;
        typedef typename Tr::Point Point;
        typedef typename Tr::Geom_traits Gt;

        // Header.
        const Tr& tr = c2t3.triangulation();

        os << "MeshVersionFormatted 1 \n"
           << "Dimension \n"
           << "3\n\n";

        CGAL_assertion(c2t3.number_of_facets() == number_of_facets_on_surface(tr));

        //os << std::setprecision(20);

        // Finite vertices coordinates.
        os << "Vertices\n"
           << tr.number_of_vertices() << " \n";

        std::map<Vertex_handle, int> V;
        int inum = 0;
        for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
            vit != tr.finite_vertices_end();
            ++vit)
        {
            V[vit] = inum++;
            Point p = static_cast<Point>(vit->point());
            os << p.x() << " " << p.y() << " " << p.z() << " 0 \n";
        }


        // Finite facets indices.
        os << "\nTriangles\n"
           << c2t3.number_of_facets() << " \n";

        for( Finite_facets_iterator fit = tr.finite_facets_begin();
            fit != tr.finite_facets_end(); ++fit)
            if ((*fit).first->is_facet_on_surface((*fit).second)==true)
            {
                for (int i=0; i<4; i++)
                    if (i != (*fit).second)
                        os << V[(*fit).first->vertex(i)]+1 << " ";

                os << "0 \n"; // without color.
            }

        // Footer
        os << "\nEnd\n";
    }

} // end namespace Surface_mesher_io

// backward compatibility: the CGAL namespace was first forgotten
// TODO: fix this
namespace CGAL {
    using namespace Surface_mesher_io;
}
using namespace Surface_mesher_io;

#endif  // CGAL_C2T3_OUTPUT_SURFACE_FACETS_TO_MEDIT_H
