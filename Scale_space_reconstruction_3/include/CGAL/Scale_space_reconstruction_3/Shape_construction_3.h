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


#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_SHAPE_CONSTRUCTION_3_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_SHAPE_CONSTRUCTION_3_H

#include <CGAL/license/Scale_space_reconstruction_3.h>


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

// provides a generalized constructor for the shape of a set of points. 
/* \ingroup PkgScaleSpaceReconstruction3Classes
 *  The shape of a set of points is ill-defined. Specifically,
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
 *  \tparam GeomTraits is the geometric traits class. It must be a model of
 *  `DelaunayTriangulationTraits_3`. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.
 *  \tparam FixedSurface determines whether the shape is constructed for a fixed
 *  scale. It must be a Boolean_tag type.
 */
template < class GeomTraits, class FixedSurface >
class Shape_construction_3 {
    typedef Triangulation_vertex_base_with_info_3< unsigned int, GeomTraits >   Vb;
    typedef Alpha_shape_vertex_base_3< GeomTraits, Vb >                     aVb;
    typedef Triangulation_cell_base_with_info_3< unsigned int, GeomTraits >     Cb;
    typedef Alpha_shape_cell_base_3< GeomTraits, Cb >                       aCb;
    typedef Triangulation_data_structure_3<aVb,aCb>                         Tds;

public:
/// \name Types
/// \{
	typedef typename GeomTraits::FT                                         FT;                 ///< defines the number field type.
	typedef typename GeomTraits::Point_3                                    Point;              ///< defines the point type.
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                  Triangulation_data_structure;                     ///< defines the triangulation data structure type.
    typedef Delaunay_triangulation_3< GeomTraits, Triangulation_data_structure >    Triangulation;  ///< defines the triangulation type.
#else
    typedef Tds                                                             Triangulation_data_structure;
    typedef Delaunay_triangulation_3< GeomTraits, Tds >                         Triangulation;
#endif // DOXYGEN_RUNNING
    typedef Alpha_shape_3< Triangulation >                                  Shape;              ///< defines the shape type.

/// \}

private:
    typedef internal::Auto_count<Point>                                     PointIndex;

public:
/// \name Constructors
/// \{
    /// constructs a default shape constructor.
    Shape_construction_3() {}
/// \}

/// \name Operations
/// \{
    /// constructs a new shape.
    /** Important note: Shape_construction_3 does not take responsibility for destroying
     *  the object after use.
     *
     *  \param shape points to the shape to base the new shape on.
     *  If `shape` is NULL, the new shape will not contain any vertices.
     *  Otherwise, the new shape will clone the vertices.
     *  \param squared_radius is the squared scale parameter of the shape.
     *  \return a pointer to the new shape.
     */
    Shape* construct( Shape* shape, const FT& squared_radius ) const {
        if( shape ) return new Shape( *shape, squared_radius, Shape::GENERAL );
        else return new Shape( squared_radius, Shape::GENERAL );
    }
    
    /// constructs a new shape.
    /** Important note: Shape_construction_3 does not take responsibility for destroying
     *  the object after use.
     *
     *  \tparam InputIterator an interator over the points.
     *  The iterator should point to a model of Point.
     *  \param begin is an iterator to the first point of the shape.
     *  \param end is a past-the-end iterator for the points.
     *  \param squared_radius is the squared scale parameter.
     *  \return a pointer to the new shape.
     */
    template < class InputIterator >
    Shape* construct( InputIterator begin, InputIterator end, const FT& squared_radius ) const {
        return new Shape( boost::make_transform_iterator( begin, PointIndex() ),
                          boost::make_transform_iterator( end, PointIndex() ),
                          squared_radius, Shape::GENERAL );
    }

    /// changes the scale of a shape.
    /** Important note: Shape_construction_3 may destroy the shape object and
     *  replace it by a new shape.
     *
     *  \param shape points to the shape to adjust.
     *  \param squared_radius is the new squared scale parameter of the shape.
     *  \pre `shape` is not NULL.
     */
    void change_scale( Shape* shape, const FT& squared_radius ) const {
        CGAL_assertion( shape != NULL );
        shape->set_alpha( squared_radius );
    }
/// \}
}; // class Shape_construction_3

// The type for the shape of a set of points with fixed scale.
template < class GeomTraits >
class Shape_construction_3 < GeomTraits, Tag_true > {
    
    typedef Triangulation_vertex_base_with_info_3< unsigned int, GeomTraits >   Vb;
    typedef Fixed_alpha_shape_vertex_base_3< GeomTraits, Vb >               aVb;
    typedef Triangulation_cell_base_with_info_3< unsigned int, GeomTraits >     Cb;
    typedef Fixed_alpha_shape_cell_base_3< GeomTraits, Cb >                 aCb;

    typedef Triangulation_data_structure_3<aVb,aCb>                         Tds;

public:
    typedef Tds                                                             Triangulation_data_structure;
    typedef Delaunay_triangulation_3< GeomTraits, Tds >                     Triangulation;
    typedef Fixed_alpha_shape_3< Triangulation >                            Shape;

	typedef typename GeomTraits::FT                                         FT;
	typedef typename GeomTraits::Point_3                                    Point;
private:
    typedef internal::Auto_count<Point>                                     PointIndex;
       
public:
    Shape_construction_3() {}

    //  Construct a new shape, possibly cloning an existing shape.
    /*  Note: Shape_construction_3 does not take responsibility for destroying
     *  the object after use.
     */
    Shape* construct( Shape* shape, const FT& squared_radius ) const {
        if( shape ) return new Shape( *shape, squared_radius );
        else return new Shape( squared_radius );
    }
    
    //  Construct a new shape.
    /*  Note: Shape_construction_3 does not take responsibility for destroying
     *  the object after use.
     */
    template < class InputIterator >
    Shape* construct( InputIterator begin, InputIterator end, const FT& squared_radius ) const {
        return new Shape( boost::make_transform_iterator( begin, PointIndex() ),
                          boost::make_transform_iterator( end, PointIndex() ),
                          squared_radius );
    }

    //  Set the scale of a shape.
    /*  Important note: Shape_construction_3 may destroy the shape and
     *  replace it by a new shape.
     */
    void change_scale( Shape*& shape, const FT& squared_radius ) const {
        CGAL_assertion( shape != NULL );
        Shape* tmp = construct( shape, squared_radius );
        delete shape;
        shape = tmp;
    }
}; // class Shape_construction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_SHAPE_CONSTRUCTION_3_H
