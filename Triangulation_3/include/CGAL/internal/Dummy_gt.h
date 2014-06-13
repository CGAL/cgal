// A dummy geometry traits class.
// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s): Thijs van Lankveld

#ifndef CGAL_INTERNAL_TRIANGULATION_DUMMY_GT_H
#define CGAL_INTERNAL_TRIANGULATION_DUMMY_GT_H

namespace CGAL { namespace internal {

// Dummy triangulation which provides all types that a triangulation can use.
struct Dummy_gt {
    struct Point_3 {};
    struct Segment_3 {};
    struct Triangle_3 {};
    struct Tetrahedron_3 {};
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_DUMMY_GT_H
