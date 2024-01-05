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
  alpha=""
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
    elif opt == "-a":
      alpha = arg
    elif opt == "-o":
      outputdir = arg
    elif opt == "-c":
      commit_hash = arg
    elif opt == "-d":
      diff_hash = arg
      do_diff = True
    elif opt == "-p":
      diffdir = arg

  all_metric = {
      "Mean_Min_Angle_(degree)" : {},
      "Mean_Max_Angle_(degree)" : {},
      "Mean_Radius_Ratio" : {},
      "Mean_Edge_Ratio" : {},
      "Mean_Aspect_Ratio" : {},
      "Complexity_(#_of_triangle)" : {},
      "#_of_almost_degenerate_triangle" : {},
      "Hausdorff_distance_output_to_input_(%_of_bbox_diag)" : {}}
  num_input = 0
  print("inputdir = ", inputdir)
  for filename in os.listdir(inputdir) :
    new_path = os.path.join(inputdir,filename)
    new_file = open(new_path)
    if do_diff :
      old_path = os.path.join(diffdir,filename)
      old_file = open(old_path)
      is_empty_old = os.path.getsize(old_path) <= 1
      for key in all_metric :
        try :
          new_val = float(new_file.readline().rstrip('\n'))
          old_val = float(old_file.readline().rstrip('\n'))
          mesh_id = str(filename.split('.')[0])
          all_metric[key][mesh_id] = [new_val, old_val]
        except ValueError:
          pass
    else :
      for key in all_metric :
        try :
          new_val = float(new_file.readline().rstrip('\n'))
          mesh_id = str(filename.split('.')[0])
          all_metric[key][mesh_id] = [new_val, new_val]
        except ValueError:
          pass
    num_input = num_input+1

  # update .pdf chart
  date_now = datetime.datetime.now()
  date_for_filename = str(date_now.year) +"_"+ str(date_now.month) +"_"+ str(date_now.day) +"_"+ str(date_now.hour) +"h"+ str(date_now.minute) +"mn"
  for key in all_metric:
    goal = 0
    if key == "Mean_Min_Angle_(degree)" or key == "Mean_Max_Angle_(degree)":
      goal = 60
    elif key == "Mean_Radius_Ratio" or key == "Mean_Edge_Ratio" or key == "Mean_Aspect_Ratio" :
      goal = 1

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
      title += "Diff between " + commit_hash + " and " + diff_hash + " on " + str(num_input) + " meshes from Thingi10K\nAlpha = Bbox diag length / " + alpha
    else :
      title += "Benchmarking on " + str(num_input) + " meshes from Thingi10K\nAlpha = Bbox diag length / " + alpha

    avg_str = str(format(abs(avg), '.2f'))
    if key == "Mean_Min_Angle_(degree)" or key == "Mean_Max_Angle_(degree)":
      title += "\nIn average we have " + avg_str + "°"
    elif key == "Mean_Radius_Ratio" or key == "Mean_Edge_Ratio" or key == "Mean_Aspect_Ratio" :
      title += "\nIn average we have a ratio of " +  avg_str
    elif key == "Hausdorff_distance_output_to_input_(%_of_bbox_diag)" :
      title += "\nIn average we have a distance of " +  avg_str + "% of bbox diag"
    elif key == "Complexity_(#_of_triangle)" or key == "#_of_almost_degenerate_triangle" :
      title += "\nIn average we have " +  avg_str + " triangles"

    if do_diff and avg_diff_to_goal == 0. :
      title += "\nNo change between the two commits"
    elif do_diff :
      avg_diff_str = str(format(abs(avg_diff_to_goal), '.2f'))
      if key == "Mean_Min_Angle_(degree)" or key == "Mean_Max_Angle_(degree)":
        if avg_diff_to_goal < 0 :
          title += "\nIn average we loose "
        else :
          title += "\nIn average we gain "
        title += avg_diff_str + "° toward 60°"
      elif key == "Mean_Radius_Ratio" or key == "Mean_Edge_Ratio" or key == "Mean_Aspect_Ratio" :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we loose "
        else :
          title += "\nIn average we gain "
        title += avg_diff_str + " of ratio toward 1"
      elif key == "Hausdorff_distance_output_to_input_(%_of_bbox_diag)" :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we increase by "
        else :
          title += "\nIn average we reduce by "
        title += avg_diff_str + " the bbox ratio"
      elif key == "Complexity_(#_of_triangle)" or key == "#_of_almost_degenerate_triangle" :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we get " + avg_diff_str + " more"
        else :
          title += "\nIn average we get " +  avg_diff_str + " less"
        title += " triangles"

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

  print("pdf updated")

  sys.exit()

if __name__ == "__main__":
  main(sys.argv[1:])
