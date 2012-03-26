# Copyright (c) 2011  GeometryFactory Sarl (France)
#
# This file is part of CGAL (www.cgal.org); you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
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
# Author(s)     : Laurent Rineau

import gdb
import itertools
import re

def lookup_function (val):
    "Look-up and return a pretty-printer that can print val."

    # Get the type.
    type = val.type

    # If it points to a reference, get the reference.
    if type.code == gdb.TYPE_CODE_REF:
        type = type.target ()

    # Get the unqualified type, stripped of typedefs.
    type = type.unqualified ().strip_typedefs ()

    # Get the type name.    
    typename = type.tag
    if typename == None:
        return None

    # Iterate over local dictionary of types to determine
    # if a printer is registered for that type.  Return an
    # instantiation of the printer if found.
    for function in CGAL_pretty_printers_dict:
        if function.search (typename):
            return CGAL_pretty_printers_dict[function] (val)
        
    # Cannot find a pretty printer.  Return None.
    return None



class CGAL_Handle_for:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        node = self.val['ptr_'].dereference()
        return 'Handle_for(%s , refcount=%d)' % (node['t'],node['count'])

class CGAL_Point_2:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        node = self.val;
        type = self.val.type
        for field in type.fields():
            if field.is_base_class:
                node = node.cast(field.type)
                break
        return 'CGAL::Point_2(%s)' % node['base']['base']

class CGAL_Tdsvb3:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        node = self.val;
        return 'CGAL::Tvb_3(%s)' % node['_p']

class CGAL_Point_3:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        node = self.val;
        type = self.val.type
        for field in type.fields():
            if field.is_base_class:
                node = node.cast(field.type)
                break
        return 'CGAL::Point_3(%s)' % node['base']['base']

class CGAL_Vector_2:
    def __init__(self, val, name):
        self.val = val
        self.name = name
        
    def to_string (self):
        node = self.val['base']
        return 'CGAL::%s(%s)' % (self.name, node)

class CGAL_Array:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        node = self.val['_M_instance']
        return node

class CGAL_Boost_tuples:
    def __init__(self, val):
        self.val = val
        
    def to_string (self):
        return '{%s}' % self.display_head_and_continue(self.val)

    def display_head_and_continue(self, val):
        has_tail = False;
        for field in val.type.fields():
            if(field.name != 'head'):
                has_tail = 1
                break
        if has_tail:
            return '%s, %s' % (val['head'], self.display_head_and_continue(val['tail']))
        else:
            return val['head']



gdb.pretty_printers.append(lookup_function)

print "Hello from CGAL_pretty_printers"
CGAL_pretty_printers_dict = {}
CGAL_pretty_printers_dict[re.compile('^CGAL::Handle_for<.*>$')] = lambda val: CGAL_Handle_for(val)
CGAL_pretty_printers_dict[re.compile('^CGAL::Point_2<.*>$')] = lambda val: CGAL_Point_2(val)
CGAL_pretty_printers_dict[re.compile('^CGAL::Point_3<.*>$')] = lambda val: CGAL_Point_3(val)
CGAL_pretty_printers_dict[re.compile('^CGAL::Vector_2<.*>$')] = lambda val: CGAL_Vector_2(val, 'Vector_2')
CGAL_pretty_printers_dict[re.compile('^CGAL::Circle_2<.*>$')] = lambda val: CGAL_Vector_2(val, 'Circle_2')

#CGAL_pretty_printers_dict[re.compile('^CGAL::Triangulation_ds_vertex_base_3<.*>$')] = lambda val: CGAL_Tdsvb3(val)

CGAL_pretty_printers_dict[re.compile('^(std|boost)(::tr1)?::array<.*>')] = lambda val: CGAL_Array(val)

CGAL_pretty_printers_dict[re.compile('^(std|boost)(::tr1)?(::tuples)?::tuple<.*>')] = lambda val: CGAL_Boost_tuples(val)
#p2 = gdb.selected_frame().read_var('p2')

