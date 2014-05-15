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

#include <CGAL/Scale_space_reconstruction_3/internal/Auto_count.h>

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

/// The standard type for the shape of a set of points.
/** The shape of a set of points is ill-defined. Specifically,
 *  because a set of points has no inherent notion of connectivity,
 *  the contour or outline of a set of points is no more descriptive than
 *  the set itself.
 *
 *  For this reason, a meaningful shape of a set of points can only be
 *  defined in the presence of some sort of indication of scale. This
 *  scale factor describes at which level of detail we wish to examine the shape
 *  of the set of points. A shape of a point set at a large scale will
 *  generally have less details than a shape of the same point set at
 *  a smaller scale.
 *
 *  The shape can be constructed either at a fixed predefined scale,
 *  or at a dynamic scale. The first option is faster when constructing
 *  a single shape. It is undefined whether a shape with fixed scale may
 *  have its scale changed, but if so, this will likely require more time
 *  than changing a dynamic scale. In either case, it is possible to change
 *  the point set while maintaining the same scale.
 *
 *  A shape is generally stored as a subset of the elements of a triangulation.
 *  \tparam Kernel the geometric traits class. It should be a model of
 *  DelaunayTriangulationTraits_3.
 *  \tparam Fixed_scale whether the shape is constructed for a fixed scale.
 *  It should be a model of Boolean_tags.
 */
template < class Kernel, class Fixed_scale >
class Shape_of_points_3 {
    typedef Triangulation_vertex_base_with_info_3< unsigned int, Kernel >   Vb;
    typedef Alpha_shape_vertex_base_3< Kernel, Vb >                         aVb;
    typedef Triangulation_cell_base_with_info_3< unsigned int, Kernel >     Cb;
    typedef Alpha_shape_cell_base_3< Kernel, Cb >                           aCb;
    typedef Triangulation_data_structure_3<aVb,aCb>                         Tds;

public:
/// \name Types
/// \{
	typedef typename Kernel::FT                                             FT;                 ///< The number field type.
	typedef typename Kernel::Point_3                                        Point;              ///< The point type.
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                  Triangulation_data_structure;                     ///< The triangulation data structure type.
    typedef Delaunay_triangulation_3< Kernel, Triangulation_data_structure >    Triangulation;  ///< The triangulation type.
#else
    typedef Tds                                                             Triangulation_data_structure;
    typedef Delaunay_triangulation_3< Kernel, Tds >                         Triangulation;
#endif // DOXYGEN_RUNNING
    typedef Alpha_shape_3< Triangulation >                                  Shape;              ///< The shape type.

/// \}

private:
    typedef internal::Auto_count<Point>                                     PointIndex;

public:
/// \name Constructors
/// \{
    /// Default constructor.
    Shape_of_points_3() {}
/// \}

/// \name Operations
/// \{
    /// Construct a new shape.
    /** Important note: Shape_of_points_3 does not take responsibility for destroying
     *  the object after use.
     *
     *  \param shape the shape to base the new shape on.
     *  If `shape` is NULL, the new shape will not contain any vertices.
     *  Otherwise, the new shape will clone the vertices.
     *  \param squared_radius the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     */
    Shape* construct( Shape* shape, const FT& squared_radius ) const {
        if( shape ) return new Shape( *shape, squared_radius, Shape::GENERAL );
        else return new Shape( squared_radius, Shape::GENERAL );
    }
    
    /// Construct a new shape.
    /** Important note: Shape_of_points_3 does not take responsibility for destroying
     *  the object after use.
     *
     *  \tparam InputIterator an interator over the points.
     *  The iterator should point to a model of Point.
     *  \param begin an iterator to the first point of the shape.
     *  \param end a past-the-end iterator for the points.
     *  \param squared_radius the squared scale parameter.
     *  \return a pointer to the new shape.
     */
    template < class InputIterator >
    Shape* construct( InputIterator begin, InputIterator end, const FT& squared_radius ) const {
        return new Shape( boost::make_transform_iterator( begin, PointIndex() ),
                          boost::make_transform_iterator( end, PointIndex() ),
                          squared_radius, Shape::GENERAL );
    }

    /// Set the scale of a shape.
    /** Important note: Shape_of_points_3 may destroy the shape and
     *  replace it by a new shape.
     *
     *  \param shape the shape to adjust.
     *  \param squared_radius the new squared scale parameter of the shape.
     *  \pre `shape` is not NULL.
     */
    void set_scale( Shape* shape, const FT& squared_radius ) const {
        CGAL_assertion( shape != NULL );
        shape->set_alpha( squared_radius );
    }
/// \}
}; // class Shape_of_points_3

// The type for the shape of a set of points with fixed scale.
template < class Kernel >
class Shape_of_points_3 < Kernel, Tag_true > {
    
    typedef Triangulation_vertex_base_with_info_3< unsigned int, Kernel >   Vb;
    typedef Fixed_alpha_shape_vertex_base_3< Kernel, Vb >                   aVb;
    typedef Triangulation_cell_base_with_info_3< unsigned int, Kernel >     Cb;
    typedef Fixed_alpha_shape_cell_base_3< Kernel, Cb >                     aCb;

    typedef Triangulation_data_structure_3<aVb,aCb>                         Tds;

public:
    typedef Tds                                                             Triangulation_data_structure;
    typedef Delaunay_triangulation_3< Kernel, Tds >                         Triangulation;
    typedef Fixed_alpha_shape_3< Triangulation >                            Shape;

	typedef typename Kernel::FT                                             FT;
	typedef typename Kernel::Point_3                                        Point;
private:
    typedef internal::Auto_count<Point>                                     PointIndex;
       
public:
    Shape_of_points_3() {}

    //  Construct a new shape, possibly cloning an existing shape.
    /*  Note: Shape_of_points_3 does not take responsibility for destroying
     *  the object after use.
     */
    Shape* construct( Shape* shape, const FT& squared_radius ) const {
        if( shape ) return new Shape( *shape, squared_radius );
        else return new Shape( squared_radius );
    }
    
    //  Construct a new shape.
    /*  Note: Shape_of_points_3 does not take responsibility for destroying
     *  the object after use.
     */
    template < class InputIterator >
    Shape* construct( InputIterator begin, InputIterator end, const FT& squared_radius ) const {
        return new Shape( boost::make_transform_iterator( begin, PointIndex() ),
                          boost::make_transform_iterator( end, PointIndex() ),
                          squared_radius );
    }

    //  Set the scale of a shape.
    /*  Important note: Shape_of_points_3 may destroy the shape and
     *  replace it by a new shape.
     */
    void set_scale( Shape* shape, const FT& squared_radius ) const {
        CGAL_assertion( shape != NULL );
        Shape* tmp = construct( shape, squared_radius );
        delete shape;
        shape = tmp;
    }
}; // class Shape_of_points_3

} // namespace CGAL

#endif // CGAL_INTERNAL_SHAPE_TYPE_H