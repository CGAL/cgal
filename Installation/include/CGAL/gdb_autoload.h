// Copyright (c) 2025  GeometryFactory Sarl (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau
//
//
// GDB Pretty Printer Auto-load Header
//
// This header embeds a .debug_gdb_scripts section that tells GDB to
// automatically load CGAL pretty printers when debugging.
//
// Include this file in any CGAL executable or library to enable
// auto-loading of the gdb pretty printers.

#ifndef CGAL_GDB_AUTOLOAD_H
#define CGAL_GDB_AUTOLOAD_H

#if defined(__ELF__) && (defined(__GNUC__) || defined(__clang__))
// ELF format (Linux, BSD) with GCC or Clang

// Use inline assembly to embed the .debug_gdb_scripts section
__asm__(
    ".pushsection \".debug_gdb_scripts\", \"MS\",@progbits,1\n"
    ".byte 1\n"  // Python script marker
    ".asciz \"share/gdb/auto-load/cgal_printers.py\"\n"
    ".popsection\n"
);

#endif // __ELF__

#endif // CGAL_GDB_AUTOLOAD_H
