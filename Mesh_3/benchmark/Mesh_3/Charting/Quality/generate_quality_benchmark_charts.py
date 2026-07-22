# Copyright (c) 2019-2023 Google LLC (USA).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org).
#
# $URL$
# $Id$
# SPDX-License-Identifier: GPL-3.0-or-later
#
#
# Author(s)     : Pierre Alliez
#                 Michael Hemmer
#                 Cedric Portaneri
#
#!/usr/bin/python

import os, sys, subprocess, datetime, time, signal, getopt
import numpy as np
import matplotlib.pyplot as plt

def main(argv):

  inputdir=""
  outputdir=""
  commit_hash=""
  do_diff=False
  diffdir=""
  diff_hash=""
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:a:o:c:d:p:')
  except getopt.GetoptError:
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-i":
      inputdir = arg
    elif opt == "-o":
      outputdir = arg
    elif opt == "-c":
      commit_hash = arg
    elif opt == "-d":
      diff_hash = arg
      do_diff = True
    elif opt == "-p":
      diffdir = arg

  print("Generating quality charts for inputdir =", inputdir)
  print("Outputdir =", outputdir)
  print("Commit hash =", commit_hash)
  print("Do diff =", do_diff)
  print("Diff hash =", diff_hash)
  print("Diffdir =", diffdir)

  all_metric = {
      "Complexity_(#_of_Vertices)" : {},
      "Complexity_(#_of_Facets)" : {},
      "Complexity_(#_of_Cells)" : {},
      "Minimum_Edge_Length_" : {},
      "Mean_Edge_Length_" : {},
      "Maximum_Edge_Length_" : {},
      "Minimum_Facet_Area_" : {},
      "Mean_Facet_Area_" : {},
      "Maximum_Facet_Area_" : {},
      "Total_Area_" : {},
      "Minimum_Facet_Angle_(degree)" : {},
      "Maximum_Facet_Angle_(degree)" : {},
      "Minimum_Cell_Volume_" : {},
      "Mean_Cell_Volume_" : {},
      "Maximum_Cell_Volume_" : {},
      "Total_Volume_" : {},
      "Minimum_Cell_Angle_(degree)" : {},
      "Mean_Cell_Angle_(degree)" : {},
      "Maximum_Cell_Angle_(degree)" : {},
      "Smallest_edge_radius_ratio" : {},
      "Smallest_radius_radius_ratio" : {},
      "Bigget_V_SMA" : {}}

  num_input = 0
  for filename in os.listdir(inputdir) :
    mesh_id = str(filename.split('.')[0])
    new_path = os.path.join(inputdir,filename)
    print("new_path?", new_path)
    new_file = open(new_path)
    is_empty_new = os.path.getsize(new_path) <= 1
    print("is_empty_new?", is_empty_new)
    if do_diff :
      old_path = os.path.join(diffdir,filename)
      old_file = open(old_path)
      is_empty_old = os.path.getsize(old_path) <= 1
      print("is_empty_old?", is_empty_old)
      for key in all_metric :
        if is_empty_new or is_empty_old :
          new_val = 0.
          old_val = 0.
        else :
          new_entry = new_file.readline().strip()
          old_entry = old_file.readline().strip()
          new_val = float(new_entry) if new_entry else 0.
          old_val = float(old_entry) if old_entry else 0.
          print("new entry, val", new_entry, new_val)
          print("old entry, val", new_entry, new_val)
          all_metric[key][mesh_id] = [new_val, old_val]
    else :
      for key in all_metric :
        if is_empty_new :
          new_val = 0.
        else :
          new_entry = new_file.readline().strip()
          new_val = float(new_entry) if new_entry else 0.
          all_metric[key][mesh_id] = [new_val, new_val]
    num_input = num_input+1

  if num_input == 0 :
    sys.exit(0)

  # update .pdf chart
  date_now = datetime.datetime.now()
  date_for_filename = str(date_now.year) +"_"+ str(date_now.month) +"_"+ str(date_now.day) +"_"+ str(date_now.hour) +"h"+ str(date_now.minute) +"mn"
  for key in all_metric:
    goal = 0
    if key == "Minimum_Facet_Angle_(degree)" or key == "Maximum_Facet_Angle_(degree)" :
      goal = 60.
    elif key == "Minimum_Cell_Angle_(degree)" or key == "Mean_Cell_Angle_(degree)" or key == "Maximum_Cell_Angle_(degree)" :
      goal = 70.5 # for a regular tetrahedron

    num_el = range(len(all_metric[key]))
    avg_diff_to_goal = 0.
    avg = 0.
    x1 = []
    x2 = []

    for value in all_metric[key].values() :
      avg += value[0]
      diff_to_goal = abs(value[1]-goal) - abs(value[0]-goal)
      avg_diff_to_goal += diff_to_goal
      x1.append(value[0])
      x2.append(value[1])
    avg_diff_to_goal /= float(len(all_metric[key]))
    avg /= float(len(all_metric[key]))

    plt.figure(figsize=(8,8))
    if do_diff :
      plt.hist(x2, bins=100, color='tab:green', alpha=0.5)
    plt.hist(x1, bins=100, color='tab:blue', alpha=0.5)
    plt.vlines(x = goal, ymin=plt.ylim()[0], ymax=plt.ylim()[1], linestyles='dashed')

    title = ""
    if do_diff :
      title += "Diff between " + commit_hash + " and " + diff_hash + " on " + str(num_input) + " meshes"
    else :
      title += "Benchmarking on " + str(num_input) + " meshes"

    avg_str = str(format(abs(avg), '.4f'))

    if key == "Complexity_(#_of_Vertices)":
      title += "\nIn average we have " +  avg_str + " vertices"
    elif key == "Complexity_(#_of_Facets)":
      title += "\nIn average we have " +  avg_str + " facets"
    elif key == "Complexity_(#_of_Cells)":
      title += "\nIn average we have " +  avg_str + " cells"
    elif key == "Minimum_Edge_Length_":
      title += "\nIn average we have a minimum edge length of " +  avg_str
    elif key == "Mean_Edge_Length_":
      title += "\nIn average we have an average edge length of " +  avg_str
    elif key == "Maximum_Edge_Length_":
      title += "\nIn average we have a maximum edge of length " +  avg_str
    elif key == "Minimum_Facet_Area_":
      title += "\nIn average we have a minimum facet area of " +  avg_str
    elif key == "Mean_Facet_Area_":
      title += "\nIn average we have an average facet area of " +  avg_str
    elif key == "Maximum_Facet_Area_":
      title += "\nIn average we have a maximum facet area of " +  avg_str
    elif key == "Total_Area_":
      title += "\nIn average we have a total area of " +  avg_str
    elif key == "Minimum_Facet_Angle_(degree)":
      title += "\nIn average we have a minimum facet angle of " + avg_str + "°"
    elif key == "Maximum_Facet_Angle_(degree)":
      title += "\nIn average we have a maximum facet angle of " + avg_str + "°"
    elif key == "Minimum_Cell_Volume_":
      title += "\nIn average we have a minimum cell volume of " +  avg_str
    elif key == "Mean_Cell_Volume_":
      title += "\nIn average we have an average cell volume of " +  avg_str
    elif key == "Maximum_Cell_Volume_":
      title += "\nIn average we have a maximum cell volume " +  avg_str
    elif key == "Total_Volume_":
      title += "\nIn average we have a total volume of " +  avg_str
    elif key == "Minimum_Cell_Angle_(degree)":
      title += "\nIn average we have a minimum dihedral angle of " + avg_str + "°"
    elif key == "Mean_Cell_Angle_(degree)":
      title += "\nIn average we have an average dihedral angle of " + avg_str + "°"
    elif key == "Maximum_Cell_Angle_(degree)":
      title += "\nIn average we have a maximum dihedral angle of " + avg_str + "°"
    elif key == "Smallest_edge_radius_ratio":
      title += "\nIn average we have a minimum edge radius ratio of " +  avg_str
    elif key == "Smallest_radius_radius_ratio":
      title += "\nIn average we have a minimum radius radius ratio of " +  avg_str
    elif key == "Bigget_V_SMA":
      title += "\nIn average we have a maximum V_SMA of " +  avg_str

    if do_diff and avg_diff_to_goal == 0. :
      title += "\nNo change between the two commits"
    elif do_diff :
      avg_diff_str = str(format(abs(avg_diff_to_goal), '.4f'))
      if key == "Minimum_Facet_Angle_(degree)" or key == "Maximum_Facet_Angle_(degree)":
        if avg_diff_to_goal < 0 :
          title += "\nIn average we lose "
        else :
          title += "\nIn average we gain "
        title += avg_diff_str + "° toward 60°"
      elif key == "Minimum_Cell_Angle_(degree)" or key == "Mean_Cell_Angle_(degree)" or key == "Maximum_Cell_Angle_(degree)":
        if avg_diff_to_goal < 0 :
          title += "\nIn average we lose "
        else :
          title += "\nIn average we gain "
        title += avg_diff_str + "° toward 70.5°"
      else :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we get " + avg_diff_str + " more"
        else :
          title += "\nIn average we get " +  avg_diff_str + " less"
        title += " " + key.replace("_"," ")

    plt.title(title, fontsize=15)
    plt.xlabel(key.replace("_"," "), fontsize=14)
    plt.ylabel("# of meshes", fontsize=14)
    plt.tick_params(axis="x", labelsize=9)
    plt.tick_params(axis="y", labelsize=9)

    chart_filename = ""
    if do_diff :
      chart_filename += "diff_"+commit_hash+"_"+diff_hash+"_"+key+"_"+date_for_filename+".pdf"
    else :
      chart_filename += "results_"+commit_hash+"_"+key+"_"+date_for_filename+".pdf"
    chart_path = os.path.join(outputdir+"/charts",chart_filename)
    if os.path.isfile(chart_path) :
      os.remove(chart_path)
    plt.savefig(chart_path, bbox_inches="tight")
    plt.close()

  print("Quality charts have been generated")

  sys.exit()

if __name__ == "__main__":
  main(sys.argv[1:])
