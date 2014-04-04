//Different types of shapes with the same API.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#ifndef CGAL_INTERNAL_SHAPE_TYPE_H
#define CGAL_INTERNAL_SHAPE_TYPE_H

#include <CGAL/internal/Auto_count.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>

#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>


namespace CGAL {

namespace internal {

/// The standard type for the shape of a set of points.
/** \tparam Kernel the geometric traits class. It should be a model of
 *  DelaunayTriangulationTraits_3.
 *  \tparam Fixed_shape whether the shape is constructed for a fixed scale.
 *  It should be a model of Boolean_tags.
 */
template < class Kernel, class Fixed_shape >
class Shape_type {
	typedef typename Kernel::FT                                                     Scalar;

    typedef CGAL::Alpha_shape_vertex_base_3< Kernel,
            CGAL::Triangulation_vertex_base_with_info_3< unsigned int, Kernel > >   Vb;
    typedef CGAL::Alpha_shape_cell_base_3< Kernel,
            CGAL::Triangulation_cell_base_with_info_3< unsigned int, Kernel > >     Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb,Cb>                             Tds;

public:
    typedef CGAL::Delaunay_triangulation_3< Kernel, Tds >                           Structure;  ///< The triangulation that spatially orders the point set.
    typedef CGAL::Alpha_shape_3< Structure >                                        Shape;      ///< The structure that identifies the triangles in the surface.
    
    /// Default constructor.
    Shape_type() {}

    /// Construct a new shape.
    /** \param shape the shape to base the new shape on.
     *  If shape is NULL, the new shape will not contain any vertices.
     *  Otherwise, the new shape will clone the vertices.
     *  \param squared_radius the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     *
     *  Note: this class does not take responsibility for destroying
     *  the object after use.
     */
    Shape* construct( Shape* shape, const Scalar& squared_radius ) {
        if( shape ) return new Shape( *shape, squared_radius, Shape::GENERAL );
        else return new Shape( squared_radius, Shape::GENERAL );
    }
    
    /// Construct a new shape.
    /** \tparam InputIterator an interator over the points of the shape.
     *  It should point to a CGAL::Point_3<Kernel>.
     *  \param start an iterator to the first point of the shape.
     *  \param end a past-the-end iterator for the points of the shape.
     *  \param squared_radius the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     *
     *  Note: this class does not take responsibility for destroying
     *  the object after use.
     */
    template < class InputIterator >
    Shape* construct( InputIterator start, InputIterator end, const Scalar& squared_radius ) {
        return new Shape( boost::make_transform_iterator( start, Auto_count<Point>() ),
                          boost::make_transform_iterator( end, Auto_count<Point>() ),
                          squared_radius, Shape::GENERAL );
    }
}; // class Shape_type

/// The yype for the shape of a set of points with fixed scale.
/** \tparam Kernel the geometric traits class. It should be a model of
 *  DelaunayTriangulationTraits_3.
 */
template < class Kernel >
class Shape_type < Kernel, CGAL::Tag_true > {
	typedef typename Kernel::FT                                                     Scalar;

    typedef CGAL::Fixed_alpha_shape_vertex_base_3< Kernel,
            CGAL::Triangulation_vertex_base_with_info_3< unsigned int, Kernel > >   Vb;
    typedef CGAL::Fixed_alpha_shape_cell_base_3< Kernel,
            CGAL::Triangulation_cell_base_with_info_3< unsigned int, Kernel > >     Cb;

    typedef CGAL::Triangulation_data_structure_3<Vb,Cb>                             Tds;

public:
    typedef CGAL::Delaunay_triangulation_3< Kernel, Tds >                           Structure;  ///< The triangulation that spatially orders the point set.
    typedef CGAL::Fixed_alpha_shape_3< Structure >                                  Shape;      ///< The structure that identifies the triangles in the surface.
       
    /// Default constructor.
    Shape_type() {}

    /// Construct a new shape.
    /** \param shape the shape to base the new shape on.
     *  If shape is NULL, the new shape will not contain any vertices.
     *  Otherwise, the new shape will clone the vertices.
     *  \param squared_radius the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     *
     *  Note: this class does not take responsibility for destroying
     *  the object after use.
     */
    Shape* construct( Shape* shape, const Scalar& squared_radius ) {
        if( shape ) return new Shape( *shape, squared_radius );
        else return new Shape( squared_radius );
    }
    
    /// Construct a new shape.
    /** \tparam InputIterator an interator over the points of the shape.
     *  It should point to a CGAL::Point_3<Kernel>.
     *  \param start an iterator to the first point of the shape.
     *  \param end a past-the-end iterator for the points of the shape.
     *  \param squared_radius the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     *
     *  Note: this class does not take responsibility for destroying
     *  the object after use.
     */
    template < class InputIterator >
    Shape* construct( InputIterator start, InputIterator end, const Scalar& squared_radius ) {
        return new Shape( boost::make_transform_iterator( start, Auto_count<Point>() ),
                          boost::make_transform_iterator( end, Auto_count<Point>() ),
                          squared_radius );
    }
}; // class Shape_type

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_SHAPE_TYPE_H