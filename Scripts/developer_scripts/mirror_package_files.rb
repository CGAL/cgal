#! /usr/bin/ruby

# == Synopsis
#
# "Mirrors" package files into a build folder for compilation.
#
# The package files are taken from a source package directory
# (typically a folder under an svn working copy)
#
# The directory structure under the package folder is replicated under
# a target build directory (creating subdirs as neccesary)
#
# "Mirroring" a package file consist on creating within the build directory a "view"
# to the file located into the package directory.
# This allows compilation within a "CGAL build directory structure" while
# keeping the source files in its original location within an svn working copy.
#
# In true POSIX platforms this is done via a symlink.
#
# In Windows platforms hardlinks are used instead. To ensure that broken hardlinks
# are properly syncronized, these are re-created whenever the source file is newer.
# 
# The mirroring operation is automatically selected according to the platform,
# but can be overrided if needed via the '--mirror_operation' option.
#
# Unless specified otherwise, via the '-all' option, internal "dont_submit" files
# are excluded.
#
# == Usage
#
# mirror_package_files [OPTIONS] PACKAGE_DIR BUILD_DIR
#
# PACKAGE_DIR  the source package directory
#
# BUILD_DIR  the destination build directory
#
# OPTIONS:
#
# -h, --help:
#    show help
#
# -o, -mirror_operation [symlink|hardlink]
#    Mirroring operation
#       symlink  is the default in platforms supporting it (*nix)
#       hardlink is the default in platforms not supporting symlinks (windows)
#    
# -a, all
#    Do not exclude internal 'dont submit' files.
#

require 'getoptlong'
require 'rdoc/usage'

load 'list_pkg_files_impl.rb'
load 'mirror_pkg_files_impl.rb'

opts = GetoptLong.new( [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
                       [ '--mirror_operation', '-o', GetoptLong::OPTIONAL_ARGUMENT ],
                       [ '--all', '-a', GetoptLong::OPTIONAL_ARGUMENT ]
                     )

exclude_internal = true

mirror_operation = default_mirror_op

opts.each do |opt, arg|
  case opt
    when '--help'
      RDoc::usage
      
    when '--mirror_operation'
      case arg
        when 'symlink'
          mirror_operation = :symlink
        when 'hardlink'
          mirror_operation = :hardlink
        else
          $stderr << "Invalid mirror operation: " << arg
      end      
      
    when '--all'
      exclude_internal = false           
       
  end
end

src_package_dir = ARGV.length >= 1 ? ARGV.shift : ''
tgt_build_dir   = ARGV.length >= 1 ? ARGV.shift : ''

unless src_package_dir.empty? || tgt_build_dir.empty? then

  excluded = exclude_internal ? pkg_dont_submit_list(src_package_dir) : ExcludedFiles.new
  
  files = list_pkg_files(excluded,src_package_dir)

  mirror_pkg_files(files,src_package_dir,tgt_build_dir,mirror_operation)
else
  RDoc::usage
end
