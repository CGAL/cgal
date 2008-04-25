# Copyright (c) 2005,2006  Utrecht University (The Netherlands),
# ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
# INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
# (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
# and Tel-Aviv University (Israel).  All rights reserved.
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
# $URL$
# $Id$
#
# Author(s)     : Daniel Russel

#---------------------------------------------------------------------#
#                    specific rules for kds
#---------------------------------------------------------------------#

vpath %.h ../../../Kinetic_data_structures/include/CGAL/Kinetic/IO/internal ../../include/CGAL/Kinetic/IO/internal

#---------------------------------------------------------------------#
#                    moc rules
#---------------------------------------------------------------------#

#KDS_Qt_examiner_viewer$(OBJ_EXT) : KDS_Qt_examiner_viewer.moc

Kinetic_Qt_widget_2_core$(OBJ_EXT): Kinetic_Qt_widget_2_core.moc

Kinetic_Qt_core$(OBJ_EXT): Kinetic_Qt_core.moc

Kinetic_Qt_timer$(OBJ_EXT): Kinetic_Qt_timer.moc

Kinetic_Qt_window_2$(OBJ_EXT): Kinetic_Qt_window_2.moc

Kinetic_pixmaps$(OBJ_EXT): Kinetic_*.xpm

# It would be nice the *.moc file base name is identical to the *.h.
# In this case the following rules can be removed all together.

Kinetic_Qt_core.moc: Qt_core.h
	$(QT_MOC) $< -o $@

Kinetic_Qt_timer.moc: Qt_timer.h
	$(QT_MOC) $< -o $@

Kinetic_Qt_widget_2_core.moc: Qt_widget_2_core.h
	$(QT_MOC) $< -o $@

Kinetic_Qt_window_2.moc: Qt_window_2.h
	$(QT_MOC) $< -o $@
