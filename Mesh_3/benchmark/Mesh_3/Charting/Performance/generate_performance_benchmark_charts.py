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

  print("Generating performance charts for inputdir =", inputdir)
  print("Outputdir =", outputdir)
  print("Commit hash =", commit_hash)
  print("Do diff =", do_diff)
  print("Diff hash =", diff_hash)
  print("Diffdir =", diffdir)

  all_metric = {
      "Facet_Scan_Time_(second)" : {},
      "Facet_Refine_Time_(second)" : {},
      "Cell_Scan_Time_(second)" : {},
      "Cell_Refine_Time_(second)" : {},
      "Lloyd_Time_(second)" : {},
      "ODT_Time_(second)" : {},
      "Perturber_Time_(second)" : {},
      "Exuder_Time_(second)" : {},
      "Total_Time_(second)" : {},
      "Memory_Peak_(mbytes)" : {}}

  num_input = 0
  for filename in os.listdir(inputdir) :
    mesh_id = str(filename.split('.')[0])
    print("perf charting, filename", filename)
    new_path = os.path.join(inputdir,filename)
    new_file = open(new_path)
    is_empty_new = os.path.getsize(new_path) <= 1
    if do_diff :
      old_path = os.path.join(diffdir,filename)
      old_file = open(old_path)
      is_empty_old = os.path.getsize(old_path) <= 1
      for key in all_metric:
        if is_empty_new or is_empty_old :
          new_val = 0.
          old_val = 0.
        else :
          new_entry = new_file.readline().strip()
          old_entry = old_file.readline().strip()
          new_val = float(new_entry) if new_entry else 0.
          old_val = float(old_entry) if old_entry else 0.
        all_metric[key][mesh_id] = [new_val, old_val]
    else :
      for key in all_metric:
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
    if key == "Memory_Peak_(mbytes)" :
      title += "\nIn average we use up to " +  avg_str + " mbytes"
    else :
      title += "\nIn average we spend " + avg_str + " seconds"

    if do_diff and avg_diff_to_goal == 0. :
      title += "\nNo change between the two commits"
    elif do_diff :
      avg_diff_str = str(format(abs(avg_diff_to_goal), '.4f'))
      if key == "Memory_Peak_(mbytes)" :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we use " +  avg_diff_str + " more"
        else :
          title += "\nIn average we use " +  avg_diff_str + " less"
        title += " mbytes"
      else :
        if avg_diff_to_goal < 0 :
          title += "\nIn average we get slower by "
        else :
          title += "\nIn average we get faster "
        title += avg_diff_str + " seconds"

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

  print("Performance charts have been generated")

  sys.exit()

if __name__ == "__main__":
  main(sys.argv[1:])
