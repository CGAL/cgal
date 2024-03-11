#!/bin/env python
# count_events.py

__date__='2015-09-03'
__author__='Gernot Walzl'

import sys
import os
import sqlite3
import re


def read_polyhedrons(database):
    result = [['PolyhedronID', 'description',
            'Vertices', 'Edges', 'Facets']]
    connection = sqlite3.connect(database)
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM Polyhedrons ORDER BY PolyhedronID ASC')
    polyhedrons = cursor.fetchall()
    for polyhedron in polyhedrons:
        poly_id = polyhedron[0]
        poly_description = polyhedron[1]

        cursor.execute('SELECT COUNT(VID) FROM Vertices WHERE PolyhedronID=?',
                (poly_id,))
        num_vertices = cursor.fetchone()[0]
        cursor.execute('SELECT COUNT(EID) FROM Edges WHERE PolyhedronID=?',
                (poly_id,))
        num_edges = cursor.fetchone()[0]
        cursor.execute('SELECT COUNT(FID) FROM Facets WHERE PolyhedronID=?',
                (poly_id,))
        num_facets = cursor.fetchone()[0]

        result.append([poly_id, poly_description,
                num_vertices, num_edges, num_facets])
    return result


def read_skeletons(database):
    result = [['SkelID', 'PolyhedronID', 'config', 'description',
            'Nodes', 'Arcs', 'Sheets',
            'E1', 'E2', 'E3', 'E4', 'E5', 'E6',
            'V-E', 'V-V I', 'V-V II', 'V-V-E I', 'V-V-E II', 'E-E', 'V-F']]
    connection = sqlite3.connect(database)
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM StraightSkeletons ORDER BY SkelID ASC')
    skeletons = cursor.fetchall()
    for skeleton in skeletons:
        skel_id = skeleton[0]
        poly_id = skeleton[1]
        skel_config = skeleton[2]
        skel_description = skeleton[3]

        cursor.execute('SELECT COUNT(NID) FROM Nodes WHERE SkelID=?',
                (skel_id,))
        num_nodes = cursor.fetchone()[0]
        cursor.execute('SELECT COUNT(AID) FROM Arcs WHERE SkelID=?',
                (skel_id,))
        num_arcs = cursor.fetchone()[0]
        cursor.execute('SELECT COUNT(SID) FROM Sheets WHERE SkelID=?',
                (skel_id,))
        num_sheets = cursor.fetchone()[0]

        cursor.execute('SELECT SkelID, etype, COUNT(etype) FROM Events ' +
                'WHERE SkelID=? ' +
                'GROUP BY SkelID, etype ORDER BY SkelID, etype ASC',
                (skel_id,))
        events = cursor.fetchall()
        num_events = [0] * 15
        for event in events:
            num_events[event[1]] = event[2]

        result.append([skel_id, poly_id, skel_config, skel_description,
                num_nodes, num_arcs, num_sheets,
                num_events[2], num_events[3], num_events[4], num_events[5], num_events[6], num_events[7],
                num_events[10], num_events[8], num_events[9], num_events[11], num_events[12], num_events[13], num_events[14]])
    return result


def regex_replace_column(mytable, column, pattern, repl):
    result = []
    for row in mytable:
        row_c = list(row)
        row_c[column] = re.sub(pattern, repl, row_c[column])
        result.append(row_c)
    return result


def select_columns(mytable, columns):
    result = []
    for row in mytable:
        row_c = []
        for index in columns:
            row_c.append(row[index])
        result.append(row_c)
    return result


def print_csv(mytable):
    for row in mytable:
        str_row = ''
        for cell in row:
            if len(str_row) > 0:
                str_row += ' | '
            str_row += str(cell)
        print(str_row)


def print_latex_table(mytable):
    num_cols = len(mytable[0])
    print('\\begin{tabular}{'+'|c'*num_cols+'|}')
    print('\\hline')
    for row in mytable:
        str_row = ''
        for cell in row:
            if len(str_row) > 0:
                str_row += ' & '
            str_row += str(cell)
        str_row = str_row.replace('_', '\_')
        str_row += ' \\\\'
        print(str_row)
    print('\\hline')
    print('\\end{tabular}')


if __name__ == '__main__':
    database = 'skeldata3d.db3'
    if len(sys.argv) > 1:
        database = sys.argv[1]
    if not os.path.isfile(database):
        print('Error: Database '+database+' does not exist.')
        exit(1)
    polyhedrons = read_polyhedrons(database)
    polyhedrons = regex_replace_column(polyhedrons, 1,
            ".*filename='../res/polyhedrons/([^;]*)';.*", r"\1")
    print_latex_table(polyhedrons)
    skeletons = read_skeletons(database)
    skeletons = regex_replace_column(skeletons, 3,
            ".*time=([^;]*);.*", r"\1")
#    print_latex_table(skeletons)
    print_latex_table(select_columns(skeletons, [1, 4, 5, 6, 3]))
    print_latex_table(select_columns(skeletons, [1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]))
