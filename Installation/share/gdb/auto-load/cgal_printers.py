# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: 2025 GeometryFactory Sarl
"""
CGAL GDB Pretty Printers

This module provides pretty printers for common CGAL types to make
debugging more convenient.

Auto-loaded via .debug_gdb_scripts section.
"""

import re
import os
import gdb

def base(val):
    """Get the base of a possibly derived CGAL type"""
    v = val
    while True:
        # print([f.name for f in v.type.fields()])
        if 'base' in [f.name for f in v.type.fields()]:
            v = v['base']
            continue
        if len(v.type.fields()) == 1 and v.type.fields()[0].is_base_class:
            v = v[v.type.fields()[0]]
            continue
        break
    # print(v.type)
    return v

class CGALPointPrinterBase:
    """Base class for CGAL 3D point pretty printers"""

    def __init__(self, val, label):
        self.val = val
        self.label = label

    def extract_coords(self):
        try:
            v = base(self.val)
            elems = v['_M_elems']
            x = float(elems[0])
            y = float(elems[1])
            z = float(elems[2])
            return x, y, z
        except (KeyError, IndexError) as e:
            return (f"<error: {e}>", None, None)

    def to_string(self):
        x, y, z = self.extract_coords()
        if x is not None and y is not None and z is not None:
            return f"{self.label}({x}, {y}, {z})"
        else:
            # Print raw value for debugging
            return f"{self.label} ERROR {x}"

    def children(self):
        # yield 'raw', str(self.val.format_string(raw=True, pretty=False))
        x, y, z = self.extract_coords()
        if x is not None and y is not None and z is not None:
            yield 'x', x
            yield 'y', y
            yield 'z', z

    def display_hint(self):
        return 'array'


class PointC3Printer(CGALPointPrinterBase):
    """Pretty printer for CGAL::PointC3<any>"""
    def __init__(self, val):
        super().__init__(val, "PointC3")


class Point3Printer(CGALPointPrinterBase):
    """Pretty printer for CGAL::Point_3<any>"""
    def __init__(self, val):
        super().__init__(val, "Point_3")


class UniqueHashMapPrinter:
    """Pretty printer for CGAL::Unique_hash_map<Key, Data, ...>"""

    class _iter:
        def __init__(self, m_map):
            self._table = m_map['table']
            self._freelist = m_map['freelist']
            self._table_end = int(m_map['table_end'])
            try:
                self._nullkey = (1 << (gdb.lookup_type('size_t').sizeof * 8)) - 1
            except gdb.error:
                self._nullkey = (1 << (gdb.lookup_type('unsigned long').sizeof * 8)) - 1
            self._index = 0
            self._item = self._table

        def __iter__(self):
            return self

        def __next__(self):
            index = self._index
            self._index += 1
            if(self._item == self._freelist):
                raise StopIteration
            key = self._item['k']
            if int(key) == self._nullkey:
                self._item = self._item + 1
                return self.__next__()
            elt = self._item.dereference()
            self._item = self._item + 1
            return ('[%d]' % index, elt)

    def __init__(self, val):
        self.val = val

    def to_string(self):
        try:
            # Get the internal chained_map
            m_map = self.val['m_map']
            
            # Get the size and default value
            if m_map['table']:
                table_size = int(m_map['table_size'])
                default_val = m_map['def']
                return f"Unique_hash_map<...> with {table_size} buckets (default={default_val})"
            else:
                return "Unique_hash_map<...> (empty)"
        except (KeyError, ValueError, gdb.error) as e:
            return f"Unique_hash_map<...> <error: {e}>"

    def children(self):
        try:
            m_map = self.val['m_map']
            
            if not m_map['table']:
                return iter(())

            return self._iter(m_map)

        except (KeyError, ValueError, gdb.error) as e:
            return f"Unique_hash_map<...> <error: {e}>"

    def display_hint(self):
        return 'array'


def cgal_lookup_function(val):
    """Lookup function to register CGAL pretty printers"""
    typename = str(val.type.strip_typedefs())

    # Match CGAL::PointC3<anything>
    if re.match(r'^CGAL::PointC3<.*>$', typename):
        return PointC3Printer(val)

    # Match CGAL::Point_3<R_>
    if re.match(r'^CGAL::Point_3<.*>$', typename):
        return Point3Printer(val)

    # Match CGAL::Unique_hash_map<Key, Data, ...>
    if re.match(r'^CGAL::Unique_hash_map<.*>$', typename):
        return UniqueHashMapPrinter(val)

    return None

# De=register existing CGAL pretty printers to avoid duplicates
def deregister_existing_cgal_printers():
    try:
        objfile = gdb.current_objfile()
    except (AttributeError, RuntimeError, gdb.error):
        objfile = None

    if objfile and hasattr(objfile, 'pretty_printers'):
        printers = objfile.pretty_printers
    elif hasattr(gdb, 'pretty_printers'):
        printers = gdb.pretty_printersq
    else:
        printers = []
    to_remove = []
    for printer in printers:
        if printer.__name__ == 'cgal_lookup_function':
            to_remove.append(printer)

    for printer in to_remove:
        printers.remove(printer)


# Register the pretty printer for GDB (manual fallback for older GDB)
def register_cgal_printer(printer):
    objfile = None
    try:
        objfile = gdb.current_objfile()
    except (AttributeError, RuntimeError, gdb.error):
        pass
    if objfile and hasattr(objfile, 'pretty_printers'):
        objfile.pretty_printers.append(printer)
    elif hasattr(gdb, 'pretty_printers'):
        gdb.pretty_printers.append(printer)
    else:
        print(
            "[CGAL] Could not register pretty printer: no pretty_printers attribute found.")


deregister_existing_cgal_printers()
register_cgal_printer(cgal_lookup_function)
# Get the full path of this script file
current_file = os.path.abspath(__file__)
print(f"CGAL pretty printers auto-loaded from {current_file}")
