"""
CGAL GDB Pretty Printers

This module provides pretty printers for common CGAL types to make
debugging more convenient.

Auto-loaded via .debug_gdb_scripts section.
"""

import gdb
import re


class PointC3Printer:
    """Pretty printer for CGAL::PointC3<*> types"""

    def __init__(self, val):
        self.val = val

    def to_string(self):
        try:
            # Navigate through the nested structure:
            # ((((pa).base).base)._M_elems)[0]
            base1 = self.val['base']
            base2 = base1['base']
            elems = base2['_M_elems']

            # Extract coordinates
            x = float(elems[0])
            y = float(elems[1])
            z = float(elems[2])

            return f"PointC3({x}, {y}, {z})"
        except (AttributeError, gdb.MemoryError, RuntimeError) as e:
            return f"PointC3(<error: {e}>)"

    def children(self):
        """Show individual coordinates"""
        try:
            base1 = self.val['base']
            base2 = base1['base']
            elems = base2['_M_elems']

            yield 'x', elems[0]
            yield 'y', elems[1]
            yield 'z', elems[2]
        except (AttributeError, gdb.MemoryError, RuntimeError):
            pass

    def display_hint(self):
        return 'array'


def cgal_lookup_function(val):
    """Lookup function to register CGAL pretty printers"""
    typename = str(val.type.strip_typedefs())

    # Match CGAL::PointC3<anything>
    if re.match(r'^CGAL::PointC3<.*>$', typename):
        return PointC3Printer(val)

    return None


# Register the pretty printer globally
# When auto-loaded via .debug_gdb_scripts, use global registration
objfile = gdb.current_objfile()
if objfile:
    # Loaded for specific objfile
    gdb.printing.register_pretty_printer(
        objfile,
        cgal_lookup_function,
        replace=True
    )
else:
    # Loaded globally (e.g., via source command)
    gdb.printing.register_pretty_printer(
        gdb,
        cgal_lookup_function,
        replace=True
    )

print(f"CGAL pretty printers loaded (objfile: {objfile})")
