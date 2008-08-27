// Copyright (c) 2007  INRIA Sophia-Antipolis (France).  All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kernel_23/include/CGAL/basic.h $
// $Id: basic.h 28567 2006-02-16 14:30:13Z lsaboret $
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_DEFAULT_ARGUMENT_H
#define CGAL_DEFAULT_ARGUMENT_H

// Default_argument is a tag that can be used to shrink mangled names and
// error messages in place of the default value of template arguments.
// It could also be used by users to specify default values to arguments which
// are not at the end of the argument list.
// It can also be useful to easily break cyclic dependencies in templates.
// Maybe we should provide a macro to disable it so as to show the arguments?
// Maybe we could document it?

// If_default_argument is a helper class which helps using this scheme.

namespace CGAL {

struct Default_argument;

template < typename Argument, typename Default >
struct If_default_argument {
    typedef Argument type;
};

template < typename Default >
struct If_default_argument <Default_argument, Default> {
    typedef Default type;
};

}

#endif // CGAL_DEFAULT_ARGUMENT_H
