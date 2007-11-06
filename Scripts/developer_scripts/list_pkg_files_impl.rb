require 'Pathname'

load 'common_impl.rb'

#
# excluded files filter
#
# uses fnmatch to filter out entries that should not be included in the package list
# this allows the 'don_submit' listing to use wildcards
#
class ExcludedFiles < Array
  def exclude? ( entry )
    self.each do |to_exclude|
      return true if File.fnmatch(to_exclude.chomp,entry)
    end
    return false
  end
end


#
# Fills the 'list' array with relative pathnames to the files in the 'package_dir' folder.
# Subdirectories are recursed.
# Directory entries mathching any string in the 'excluded' array are skipped.
#
def list_pkg_files(filter, package_dir, list = [] )
  __aux_list_pkg_files(filter, as_pathname(package_dir), as_pathname(''), list)  
end

def __aux_list_pkg_files(filter, package_dir, local_subdir, list = [] )

  assert_exist!(package_dir, 'package directory' )
  
  filter = ExcludedFiles.new unless filter
  
  subdir_from_cwd = package_dir + local_subdir ;
 
  # Process each entry in the package's folder
  Dir.foreach(subdir_from_cwd) do |entry|
 
    unless entry == '.' || entry == '..'  
    
      unless filter.exclude?(entry)
      
        path_from_cwd     = subdir_from_cwd + entry 
        path_from_package = local_subdir + entry
        
        # Recurse if this entry is a sub folder
        if FileTest.directory?( path_from_cwd )
          __aux_list_pkg_files(filter, package_dir, path_from_package, list)
        else
          list << path_from_package
        end
        
      end
    end  
  end
   
  return list 
end


#
# Returns an array listing the entries that should be excluded from the package file list.
# This includes fixed entries hard-coded here and those read from a 'dont_submit' text file, if exits.
#
def pkg_dont_submit_list(package_dir)
 
  list = ['TODO',
         'dont_submit',
         'maintainer',
         'description.txt',
         'changes.txt',
         'doc_tex',
         '.svn'
        ]
                        
  dont_submit_file = package_dir + '/dont_submit' 
 
  list += IO.readlines(dont_submit_file) if FileTest.exist?(dont_submit_file) 
   
  return ExcludedFiles.new(list)
 
end

#puts pkg_dont_submit_list($test_pkg_dir)
#puts list_pkg_files(pkg_dont_submit_list($test_pkg_dir), $test_pkg_dir)
#puts list_pkg_files(nil, $test_pkg_dir)
#puts list_pkg_files(nil, 'unexitent_folder')

