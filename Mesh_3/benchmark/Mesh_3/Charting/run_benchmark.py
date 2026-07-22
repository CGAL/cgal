#!/usr/bin/python

import os, sys, subprocess, datetime, time, signal, getopt

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import XMLParser

####################################################################################################
def signal_handler(signum, frame):
  raise Exception("Timed out!")

####################################################################################################
def run_benchmark(execname, filename, facet_size, facet_approx, facet_angle, cell_size, cell_shape, max_time, test_ID, outpath):

  exit_codes = {
   0 : "VALID_SOLID_OUTPUT",
   1 : "INPUT_IS_INVALID",
   2 : "OUTPUT_IS_INVALID",
   3 : "SIGSEGV",
   4 : "SIGABRT",
   5 : "SIGFPE",
   6 : "TIMEOUT"
  }

  raw_filename = os.path.splitext(os.path.basename(filename))[0]
  xml_filename = outpath + "/logs/" + test_ID + "/" + raw_filename + ".xml"
  print("xml_filename = ", xml_filename)

  exit_code = 0
  output = ""
  cmd = ("/usr/bin/time", "-v", execname, xml_filename, filename, facet_size, facet_approx, facet_angle, cell_size, cell_shape)

  print("cmd = ", cmd)

  proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)

  try:
    outs, errs = proc.communicate(timeout=int(max_time))
    exit_code = proc.returncode
    output = outs.decode("utf-8") + errs.decode("utf-8")

    print("output = ", output)

    for output_line in output.split("\n"):
      if output_line == "Command terminated by signal 11":
        exit_code = 3
        continue
      elif output_line == "Command terminated by signal 6":
        exit_code = 4
        continue
      elif output_line == "Command terminated by signal 8":
        exit_code = 5
        continue

  except subprocess.TimeoutExpired:
    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    exit_code = 6
    output = "process ran too long"

  print("EXIT CODE: ", exit_codes[exit_code])

  outfile = outpath + "/Robustness/results/" + test_ID + "/" + raw_filename + ".txt"
  print("writing to", outfile)
  file = open(outfile, "w")
  file.write(exit_codes[exit_code])
  file.close()

  return exit_code

####################################################################################################
def parse_xml_file(filename, tag):
  tree = ET.parse(filename)
  root = tree.getroot()

  elems = root.findall(f'.//{tag}')
  if len(elems) == 0:
    print("Error: no elements found for ", tag)
    return sys.exit(1)

  elem = elems[0] # take the first
  return elem.text.strip()

####################################################################################################
def main(argv):
  execname=""
  filename=""
  facet_size=""
  facet_approx=""
  facet_angle=""
  cell_size=""
  cell_shape=""
  max_time=""
  test_ID=""
  outpath=""
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:t:', ["exec=", "facet_size=", "facet_approx=", "facet_angle=", "cell_size=", "cell_shape=", "test_ID=", "out="])
  except getopt.GetoptError:
    sys.exit(2)
  for opt, arg in opts:
    if opt == "--exec": # executable
      execname = arg
    elif opt == "-i": # filename
      filename = arg
    elif opt == "--facet_size":
      facet_size = arg
    elif opt == "--facet_approx":
      facet_approx = arg
    elif opt == "--facet_angle":
      facet_angle = arg
    elif opt == "--cell_size":
      cell_size = arg
    elif opt == "--cell_shape":
      cell_shape = arg
    elif opt == "-t":
      max_time = arg
    elif opt == "--test_ID":
      test_ID = arg
    elif opt == "--out": # output file
      outpath = arg

  print("execname = ", execname)
  print("filename = ", filename)
  print("facet_size = ", facet_size)
  print("facet_approx = ", facet_approx)
  print("facet_angle = ", facet_angle)
  print("cell_size = ", cell_size)
  print("cell_shape = ", cell_shape)
  print("max_time = ", max_time)
  print("test_ID = ", test_ID)
  print("outpath = ", outpath)

  exit_code = run_benchmark(execname, filename, facet_size, facet_approx, facet_angle, cell_size, cell_shape, max_time, test_ID, outpath)

  # if the exit code is different from 0, then there is nothing to analyze
  if exit_code != 0:
    sys.exit(exit_code)

  raw_filename = os.path.splitext(os.path.basename(filename))[0]
  xml_filename = outpath + "/logs/" + test_ID + "/" + raw_filename + ".xml"

  # Parse the XML output to extract performance and quality metrics
  print("parsing", xml_filename)

  parser = XMLParser(encoding="utf-8")
  try:
    ET.parse(xml_filename, parser)
    print("XML is valid")
  except Exception as e:
    print("XML is invalid -", e)

  # --- Performance
  perf_results_filename = outpath + "/Performance/results/" + test_ID + "/" + raw_filename + ".txt"
  perf_results = open(perf_results_filename, "w")

  # Refinement
  facet_scan_time = parse_xml_file(xml_filename, "Facets_scan_time")
  facet_refine_time = parse_xml_file(xml_filename, "Facets_refine_time")
  cell_scan_time = parse_xml_file(xml_filename, "Cells_scan_time")
  cell_refine_time = parse_xml_file(xml_filename, "Cells_refine_time")

  # Optimization
  lloyd_optim_time = parse_xml_file(xml_filename, "Lloyd_optim_time")
  odt_optim_time = parse_xml_file(xml_filename, "Odt_optim_time")
  perturber_optim_time = parse_xml_file(xml_filename, "Perturber_optim_time")
  exuder_optim_time = parse_xml_file(xml_filename, "Exuder_optim_time")

  # Total
  total_time = parse_xml_file(xml_filename, "Total_time")

  # Memory
  memory = parse_xml_file(xml_filename, "Mem")

  perf_results.write(facet_scan_time + "\n")
  perf_results.write(facet_refine_time + "\n")
  perf_results.write(cell_scan_time + "\n")
  perf_results.write(cell_refine_time + "\n")
  perf_results.write(lloyd_optim_time + "\n")
  perf_results.write(odt_optim_time + "\n")
  perf_results.write(perturber_optim_time + "\n")
  perf_results.write(exuder_optim_time + "\n")
  perf_results.write(total_time + "\n")
  perf_results.write(memory + "\n")

  perf_results.close()

  # --- Quality
  qual_results_filename = outpath + "/Quality/results/" + test_ID + "/" + raw_filename + ".txt"
  qual_results = open(qual_results_filename, "w")

  number_of_vertices = parse_xml_file(xml_filename, "V")
  number_of_facets = parse_xml_file(xml_filename, "F")
  number_of_cells = parse_xml_file(xml_filename, "C")

  min_edge_size = parse_xml_file(xml_filename, "Minimum_edge_length")
  mean_edge_size = parse_xml_file(xml_filename, "Mean_edge_length")
  max_edge_size = parse_xml_file(xml_filename, "Maximum_edge_length")
  min_facet_size = parse_xml_file(xml_filename, "Minimum_facet_area")
  mean_facet_size = parse_xml_file(xml_filename, "Mean_facet_area")
  max_facet_size = parse_xml_file(xml_filename, "Maximum_facet_area")
  total_area = parse_xml_file(xml_filename, "Total_area")
  # min_facet_distance = parse_xml_file(xml_filename, "Min_facet_distance")
  # mean_facet_distance = parse_xml_file(xml_filename, "Mean_facet_distance")
  # max_facet_distance = parse_xml_file(xml_filename, "Mean_facet_distance")
  min_facet_angle = parse_xml_file(xml_filename, "Minimum_facet_angle")
  max_facet_angle = parse_xml_file(xml_filename, "Maximum_facet_angle")
  min_cell_size = parse_xml_file(xml_filename, "Minimum_cell_volume")
  mean_cell_size = parse_xml_file(xml_filename, "Mean_cell_volume")
  max_cell_size = parse_xml_file(xml_filename, "Maximum_cell_volume")
  total_volume = parse_xml_file(xml_filename, "Total_volume")
  min_cell_shape = parse_xml_file(xml_filename, "Minimum_dihedral_angle")
  mean_cell_shape = parse_xml_file(xml_filename, "Mean_dihedral_angle")
  max_cell_shape = parse_xml_file(xml_filename, "Maximum_dihedral_angle")
  smallest_edge_radius_ratio = parse_xml_file(xml_filename, "Smallest_edge_radius_ratio")
  smallest_radius_radius_ratio = parse_xml_file(xml_filename, "Smallest_radius_radius_ratio")
  biggest_v_sma = parse_xml_file(xml_filename, "Biggest_V_SMA")

  print("number_of_vertices = ", number_of_vertices)
  print("number_of_facets = ", number_of_facets)
  print("number_of_cells = ", number_of_cells)
  print("min_edge_size = ", min_edge_size)
  print("mean_edge_size = ", mean_edge_size)
  print("max_edge_size = ", max_edge_size)
  print("min_facet_size = ", min_facet_size)
  print("mean_facet_size = ", mean_facet_size)
  print("max_facet_size = ", max_facet_size)
  print("total_area = ", total_area)
  # print("min_facet_distance = ", min_facet_distance)
  # print("mean_facet_distance = ", mean_facet_distance)
  # print("max_facet_distance = ", max_facet_distance)
  print("min_facet_angle = ", min_facet_angle)
  print("max_facet_angle = ", max_facet_angle)
  print("min_cell_size = ", min_cell_size)
  print("mean_cell_size = ", mean_cell_size)
  print("max_cell_size = ", max_cell_size)
  print("total_volume = ", total_volume)
  print("min_cell_shape = ", min_cell_shape)
  print("mean_cell_shape = ", mean_cell_shape)
  print("max_cell_shape = ", max_cell_shape)
  print("smallest_edge_radius_ratio = ", smallest_edge_radius_ratio)
  print("smallest_radius_radius_ratio = ", smallest_radius_radius_ratio)
  print("biggest_v_sma = ", biggest_v_sma)

  qual_results.write(number_of_vertices + "\n")
  qual_results.write(number_of_facets + "\n")
  qual_results.write(number_of_cells + "\n")
  qual_results.write(min_edge_size + "\n")
  qual_results.write(mean_edge_size + "\n")
  qual_results.write(max_edge_size + "\n")
  qual_results.write(min_facet_size + "\n")
  qual_results.write(mean_facet_size + "\n")
  qual_results.write(max_facet_size + "\n")
  qual_results.write(total_area + "\n")
  qual_results.write(min_facet_angle + "\n")
  qual_results.write(max_facet_angle + "\n")
  qual_results.write(min_cell_size + "\n")
  qual_results.write(mean_cell_size + "\n")
  qual_results.write(max_cell_size + "\n")
  qual_results.write(total_volume + "\n")
  qual_results.write(min_cell_shape + "\n")
  qual_results.write(mean_cell_shape + "\n")
  qual_results.write(max_cell_shape + "\n")
  qual_results.write(smallest_edge_radius_ratio + "\n")
  qual_results.write(smallest_radius_radius_ratio + "\n")
  qual_results.write(biggest_v_sma + "\n")

  qual_results.close()

if __name__ == "__main__":
  main(sys.argv[1:])
