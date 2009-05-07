// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

double random_in(const double a,
                 const double b)
{
    double r = rand() / (double)RAND_MAX;
    return a + (b - a) * r;
}

template <class K>
typename K::Point_3 random_point_in(const CGAL::Bbox_3& bbox)
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(bbox.xmin(),bbox.xmax());
    FT y = (FT)random_in(bbox.ymin(),bbox.ymax());
    FT z = (FT)random_in(bbox.zmin(),bbox.zmax());
    return typename K::Point_3(x,y,z);
}

template <class K>
typename K::Vector_3 random_vector()
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(0.0,1.0);
    FT y = (FT)random_in(0.0,1.0);
    FT z = (FT)random_in(0.0,1.0);
    return typename K::Vector_3(x,y,z);
}
