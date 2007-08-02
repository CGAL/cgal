// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_MODULAR_TRAITS_H
#define CGAL_SQRT_EXTENSION_MODULAR_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Modular.h>
#include <CGAL/Modular_traits.h>

CGAL_BEGIN_NAMESPACE

/////////// MODULAR_TRAITS BEGIN

template< class COEFF, class ROOT>
class Modular_traits< Sqrt_extension<COEFF,ROOT> > {
    
private:
    typedef Sqrt_extension<COEFF,ROOT> EXT; 
    typedef Modular_traits<COEFF> MT_COEFF;
    typedef Modular_traits<ROOT>  MT_ROOT;
    typedef typename MT_COEFF::Modular_NT Modular_NT_coeff;
    typedef typename MT_ROOT::Modular_NT  Modular_NT_root;
public:
    typedef Sqrt_extension<COEFF, ROOT > NT;
    typedef typename MT_COEFF::Is_modularizable Is_modularizable;
    typedef Sqrt_extension<Modular_NT_coeff, Modular_NT_root> Modular_NT;
    
    struct Modular_image{
        Modular_NT operator()(const EXT& a){
            typename MT_ROOT::Modular_image mod_image_root;
            typename MT_COEFF::Modular_image mod_image_coeff;
            Modular_NT_root  root_mod = mod_image_root(a.root());
            if(root_mod != Modular_NT_root(0)){
                return Modular_NT(mod_image_coeff(a.a0()),
                                  mod_image_coeff(a.a1()),
                                  root_mod);
            }else{
                return Modular_NT(mod_image_coeff(a.a0()));
            }
        }
    };

    struct Modular_image_inv{
        NT operator()(const Modular_NT& a){
            typename MT_ROOT::Modular_image_inv mod_image_inv_root;
            typename MT_COEFF::Modular_image_inv mod_image_inv_coeff;
            
            if(a.is_extended()){
                return NT(
                        mod_image_inv_coeff(a.a0()),
                        mod_image_inv_coeff(a.a1()),
                        mod_image_inv_root(a.root()));
            }else{
                return NT(mod_image_inv_coeff(a.a0()));
            }
        }
    };
};

/////////// MODULAR_TRAITS END

CGAL_END_NAMESPACE

#endif
