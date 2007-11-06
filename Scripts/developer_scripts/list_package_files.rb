#! /bin/ruby

# == Synopsis
#
# Prints a list of all files corresponding to a given package.
# Each pathname in the list is relative to the package directory.
# Unless specified otherwise, internal files (such as those listed
# in a 'dont_submit' text file) are excluded.
#
# == Usage
#
# list_package_files [OPTIONS]
#
# OPTIONS:
#
# -h, --help:
#    show help
#
# -d, -package_dir DIR:
#    Directory where the package exist. Default is the current directory.
#
# -a, --all_files
#    Do not exclude internal 'dont submit' files.
#

require 'getoptlong'
require 'rdoc/usage'

load 'list_pkg_files_impl.rb'

opts = GetoptLong.new( [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
                       [ '--package_dir', '-d', GetoptLong::OPTIONAL_ARGUMENT ],
                       [ '--all_files', '-a', GetoptLong::OPTIONAL_ARGUMENT ]
                     )

package_dir = '.'

exclude_internal = true

opts.each do |opt, arg|
  case opt
   when '--help'
      RDoc::usage
    when '--package_dir'
      package_dir = arg 
    when '--all_files'
       exclude_internal = false           
  end
end

dont_submit = exclude_internal ? pkg_dont_submit_list(package_dir) : ExcludedFiles.new()

puts list_pkg_files(dont_submit,package_dir)





     