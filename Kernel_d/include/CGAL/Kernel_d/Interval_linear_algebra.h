// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus, Olivier Devillers

#ifndef CGAL_INTERVAL_LINEAR_ALGEBRA_H
#define CGAL_INTERVAL_LINEAR_ALGEBRA_H

#include <CGAL/Interval_nt.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {

namespace internal
{

template<class Matrix>
Sign // The matrix is row major: M[i] represents row i.
/*
FIXME : the function DOES MODIFY the matrix M, but calling functions
assume M is not modified --> BUGS. (e.g. Side_of_oriented_subsphere)
*/
sign_of_determinantDxD_with_interval_arithmetic(Matrix & M)
// attempts to compute the determinant using interval arithmetic
{
    typedef Interval_nt<false> NT; // Field number type
    int sign = 1;
    int curdim = M.dimension().first;
    for(int col = 0; col < curdim; ++col)
    {
        int pivot_line = -1;
        NT pivot = 0.0;
        for(int l = col; l < curdim; ++l)
        {
            NT candidate = CGAL::abs(M[l][col]);
            if( CGAL::certainly(0.0 < candidate) )
                if( pivot.inf() < candidate.inf() )
                {
                    pivot_line = l;
                    pivot = candidate;
                }
        }
        if( -1 == pivot_line )
        {
            throw CGAL::Uncertain_conversion_exception("undecidable interval computation of determinant");
        }
        // if the pivot value is negative, invert the sign
        if( M[pivot_line][col] < 0.0 )
            sign = - sign;
        // put the pivot cell on the diagonal
        if( pivot_line > col )
        {
            std::swap(M[col], M[pivot_line]); // this takes constant time with std::vector<>
            sign = - sign;
        }
        // using matrix-line operations, set zero in all cells below the pivot cell.
        // This costs roughly (curdim-col-1)^2 mults and adds, because we don't actually
        // compute the  zeros below the pivot cell
        NT denom = NT(1.0) / M[col][col];
        for(int line = col + 1; line < curdim; ++line)
        {
            NT numer = M[line][col] * denom;
            for (int c = col + 1; c < curdim; ++c)
                M[line][c] -= numer * M[col][c];
        }
    }
    return Sign(sign);
}

} // end of namespace internal

template<>
inline
Sign
Linear_algebraCd<Interval_nt_advanced>::sign_of_determinant(const Matrix & M)
{
    switch( M.dimension().first )
    {
        case 2:
            return CGAL::sign_of_determinant( M(0,0), M(0,1),
                    M(1,0), M(1,1));
            break;
        case 3:
            return CGAL::sign_of_determinant( M(0,0), M(0,1), M(0,2),
                    M(1,0), M(1,1), M(1,2),
                    M(2,0), M(2,1), M(2,2));
            break;
        case 4:
            return CGAL::sign_of_determinant( M(0,0), M(0,1), M(0,2), M(0,3),
                    M(1,0), M(1,1), M(1,2), M(1,3),
                    M(2,0), M(2,1), M(2,2), M(2,3),
                    M(3,0), M(3,1), M(3,2), M(3,3));
            break;
        default:
            return internal::sign_of_determinantDxD_with_interval_arithmetic(const_cast<Matrix &>(M));
            break;
    }
}

} //namespace CGAL

#endif // CGAL_INTERVAL_LINEAR_ALGEBRA_H
