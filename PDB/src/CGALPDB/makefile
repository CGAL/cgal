# Copyright (c) 2007 INRIA Sophia-Antipolis (France).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org); you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; version 2.1 of the License.
# See the file LICENSE.LGPL distributed with CGAL.
#
# Licensees holding a valid commercial license may use this file in
# accordance with the commercial license agreement provided with the software.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# $URL: svn+ssh://lrineau@scm.gforge.inria.fr/svn/cgal/trunk/Installation/src/CGAL/makefile $
# $Id: makefile 36165 2007-02-09 23:33:22Z lrineau $
#
# Author(s)     : Laurent Rineau


#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

OBJECTS = all_files$(OBJ_EXT)

OLD_OBJECTS= Atom$(OBJ_EXT)\
Error_logger$(OBJ_EXT)\
Monomer_data$(OBJ_EXT)\
Chain$(OBJ_EXT)\
Model$(OBJ_EXT)\
PDB$(OBJ_EXT)\
Heterogen$(OBJ_EXT)\
Chain_pdb$(OBJ_EXT)\
Monomer$(OBJ_EXT)\
Quaternion$(OBJ_EXT)\
pdb_utils$(OBJ_EXT)


# out for now align$(OBJ_EXT)
#---------------------------------------------------------------------#
#                    libCGALPDB settings
#---------------------------------------------------------------------#

SOVERSION=2.0.0
SOMAJOR=2

LIB_NAME_BASE = CGALPDB


SHARED_LIB_ADDITIONAL_LDFLAGS =  $(CGAL_LIBPATHFLAGS)  $(CGAL_CGAL_LDFLAGS)  


#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

