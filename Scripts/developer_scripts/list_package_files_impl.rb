require 'Pathname'

require 'common_impl.rb'

# The list can contain wildcards
def matches_any_filename_in_list(filename,list)
  list.each do |item|
    return true if File.fnmatch(item.chomp,filename)
  end
  return false
end

#
# Returns an array listing the entries that should be excluded from the package file list.
# This includes fixed entries hard-coded here and those read from a 'dont_submit' text file, if exits.
#
def get_internal_package_files_list(package_dir)
 
  list = ['TODO',
         'dont_submit',
         'maintainer',
         'description.txt',
        ]
                        
  dont_submit_file = package_dir + '/dont_submit' 
 
  list += IO.readlines(dont_submit_file) if File.exist?(dont_submit_file) 
   
  return list ;
 
end



#
# Fills the 'list' array with relative pathnames to the files in the 'package_dir' folder.
# Subdirectories are recursed.
# Directory entries mathching any filename in 'internal_files' are skipped.
#
def list_package_files(package_dir, include_internal = false )

  internal_files = include_internal ? [] : get_internal_package_files_list(package_dir)  
  
  __aux_list_package_files(package_dir, '', [], internal_files )  
  
end

def __aux_list_package_files(package_dir, local_subdir, list, internal_files )

  assert_exist!(package_dir, 'package directory' )
  
  subdir_from_cwd = package_dir + '/' + local_subdir ;
 
  # Process each entry in the package's folder
  Dir.foreach(subdir_from_cwd) do |entry|
 
    entry.chomp!
    
    unless entry == '.' || entry == '..'  
    
      unless matches_any_filename_in_list(entry, internal_files)
      
        path_from_cwd     = subdir_from_cwd + '/' + entry 
        path_from_package = local_subdir    + '/' + entry
        
        # Recurse if this entry is a sub folder
        if File.directory?( path_from_cwd )
          __aux_list_package_files(package_dir, path_from_package, list, internal_files )
        else
          list << path_from_package
        end
        
      end
    end  
  end
   
  return list 
end


# -- TEST -- 
#puts get_internal_package_files_list(TEST_PKG_DIR)
#puts list_package_files(TEST_PKG_DIR, false)
#puts list_package_files(TEST_PKG_DIR)
#puts list_package_files('inexistent_folder')
# -- TEST -- 

