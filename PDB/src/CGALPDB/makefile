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

OBJECTS = \
PDB_Atom$(OBJ_EXT) \
PDB_Protein_pdb$(OBJ_EXT) \
PDB_align$(OBJ_EXT) \
PDB_Error_logger$(OBJ_EXT) \
PDB_Quaternion$(OBJ_EXT) \
PDB_distance$(OBJ_EXT) \
PDB_Model$(OBJ_EXT) \
PDB_Residue$(OBJ_EXT) \
PDB_pdb_utils$(OBJ_EXT) \
PDB_PDB$(OBJ_EXT) \
PDB_Residue_data$(OBJ_EXT) \
PDB_transforms$(OBJ_EXT) \
PDB_Protein$(OBJ_EXT) \
PDB_Transform$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    libCGALPDB settings
#---------------------------------------------------------------------#

SOVERSION=2.0.1
SOMAJOR=2

LIB_NAME_BASE = CGALPDB


SHARED_LIB_ADDITIONAL_LDFLAGS =  $(CGAL_LIBPATHFLAGS)  $(CGAL_CGAL_LDFLAGS)  


#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

