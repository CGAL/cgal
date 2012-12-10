#! /usr/bin/ruby

# == Synopsis
#
# Prints a list of all files corresponding to a given package.
# Each pathname in the list is relative to the package directory.
# Unless specified otherwise, internal files (such as those listed
# in a 'dont_submit' text file) are excluded.
#
# == Usage
#
# list_package_files [OPTIONS] package_folder
#
# OPTIONS:
#
# -h, --help:
#    show help
#
# -i, --include_internal
#    Include internal files (excluded by default).
#


require 'getoptlong'
require 'rdoc/usage'

require 'list_package_files_impl.rb'

# -- TEST -- 
# ARGV = [TEST_PKG_DIR]
# ARGV = ['-i',TEST_PKG_DIR]
# -- TEST -- 

opts = GetoptLong.new( [ '--help'            , '-h', GetoptLong::NO_ARGUMENT ],
                       [ '--include_internal', '-i', GetoptLong::NO_ARGUMENT ]
                     )


include_internal = false

opts.each do |opt, arg|
  case opt
   when '--help'
      RDoc::usage
   when '--include_internal'
      include_internal = true           
  end
end

if ARGV.length > 0
  package_dir = ARGV.shift
  puts list_package_files(package_dir, include_internal)
else
  RDoc::usage
end
