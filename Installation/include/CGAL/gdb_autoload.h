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
    ".ascii \"\\4gdb.inlined-script.CGAL\\n\"\n"
    ".ascii \"import gdb\\n\"\n"
    ".ascii \"path=[line.strip() for line in gdb.execute('info sources', to_string=True).split(',') if line.count('include/CGAL/gdb_autoload.h')][0].replace('/include/CGAL/gdb_autoload.h', '')\\n\"\n"
    ".ascii \"gdb.execute('source ' + path + '/share/gdb/auto-load/cgal_printers.py')\\n\"\n"
    ".popsection\n"
);

namespace CGAL {

  // This is only here to have a CGAL symbol in the object file,
  // so the gdb 'info sources' command can find this file's path.
  [[maybe_unused]] inline static int gdb_autoload_dummy_instance = 0;

}

#endif // __ELF__

#endif // CGAL_GDB_AUTOLOAD_H
