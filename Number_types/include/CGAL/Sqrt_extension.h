// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//                 Ron Wein         <wein@post.tau.ac.il>


#ifndef CGAL_SQRT_EXTENSION_H
#define CGAL_SQRT_EXTENSION_H

// COMMENTS FROM EXACUS
/*! \ingroup NiX_Sqrt_extension
\brief represents an extension of a number type by one square root. 

 An instance of this class
represents  an extension of the type NT by a square root of the 
type ROOT. In case NT and ROOT do not coincide, 
NT must be constructible from ROOT.  The number type NT 
must be at least a model of the IntegralDomainWithoutDiv concept. 

An Sqrt_extension is a model of RealComparable if NT is RealComparable.\n  
The <B>algebraic type</B> of NiX::Sqrt_extension depends on the algebraic type 
of NT: 
- IntegralDomainWithoutDiv -> IntegralDomainWithoutDiv
- IntegralDomain           -> IntegralDomain
- UFDomain                 -> IntegralDomain
- EuclideanRing            -> IntegralDomain
- Field                    -> Field
- FieldWithSqrt            -> Field (not recommended)


Note that NT and ROOT can themselves be an instance of
NiX::Sqrt_extension, yielding a nested extension.\n
Note that the extension of an UFDomain or EuclideanRing is just an 
IntegralDomain, since the extension in general destroys the unique 
factorization property. 
*/

#include <CGAL/number_type_basic.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>
#include <CGAL/Sqrt_extension/Algebraic_structure_traits.h>
#include <CGAL/Sqrt_extension/Real_embeddable_traits.h>
#include <CGAL/Sqrt_extension/Fraction_traits.h>
#include <CGAL/Sqrt_extension/Coercion_traits.h>
#include <CGAL/Sqrt_extension/Modular_traits.h>
#include <CGAL/Sqrt_extension/Scalar_factor_traits.h>
#include <CGAL/Sqrt_extension/Algebraic_extension_traits.h>
#include <CGAL/Sqrt_extension/Chinese_remainder_traits.h>
#include <CGAL/Sqrt_extension/io.h>
#include <CGAL/Sqrt_extension/Get_arithmetic_kernel.h>
#include <CGAL/Sqrt_extension/convert_to_bfi.h>
#include <CGAL/Sqrt_extension/Wang_traits.h>
#include <CGAL/Sqrt_extension/Eigen_NumTraits.h>


#endif  // CGAL_SQRT_EXTENSION_H

// EOF
