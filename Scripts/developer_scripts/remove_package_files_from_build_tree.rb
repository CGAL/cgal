#! /usr/bin/ruby

# == Synopsis
#
# Remove package files from a build folder.
#
# The list of files to remove are taken from a source package directory
# (typically a folder under an svn working copy) 
# (see list_package_files.rb script)
#
# == Usage
#
# remove_package_files_from_buid_tree [OPTIONS] PACKAGE_SUBDIR BUILD_ROOT
#
# PACKAGE_SUBDIR  the source package sub directory used to generate the file list
#
# BUILD_ROOT  the destination build root directory where the files will be removed
#
# OPTIONS:
#
# -h, --help:
#    show help
#
# -r, -rename
#    Rename files instead, adding a ".removed" suffix
#    

require 'getoptlong'
require 'rdoc/usage'

require 'list_package_files_impl.rb'
require 'remove_package_files_from_build_tree_impl.rb'

# -- TEST -- 
 ARGV = [TEST_PKG_DIR, TEST_BUILD_ROOT]
# ARGV = ['-r','-i',TEST_PKG_DIR, TEST_BUILD_ROOT]
# -- TEST -- 

opts = GetoptLong.new( [ '--help'            , '-h', GetoptLong::NO_ARGUMENT ],
                       [ '--rename'          , '-r', GetoptLong::NO_ARGUMENT ],
                       [ '--include_internal', '-i', GetoptLong::NO_ARGUMENT ]
                     )

rename_instead   = false
include_internal = false

opts.each do |opt, arg|
  case opt
    when '--help'
      RDoc::usage
      
    when '--rename'
      rename_instead = true
      
    when '--include_internal'
       include_internal = true           
  end
end

src_package_subdir = ARGV.length >= 1 ? ARGV.shift : ''
tgt_build_root     = ARGV.length >= 1 ? ARGV.shift : ''

unless src_package_subdir.empty? || tgt_build_root.empty? then

  files = list_package_files(src_package_subdir,include_internal)

  remove_package_files_from_build_tree(files,src_package_subdir,tgt_build_root,:rename_instead => rename_instead )
else
  RDoc::usage
end
