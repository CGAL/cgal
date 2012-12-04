require 'Pathname'
require 'FileUtils'

require 'common_impl.rb'
require 'list_package_files_impl.rb'

  
#
# Removes (or renames) a list of 'files' from a given 'build_root'
# Entries in the file list are pathnames relative to the 'package_subdir'
#
def remove_package_files_from_build_tree(files,
										 package_subdir,
										 build_root,
										 options = { :rename_instead => false }
										)

  assert_exist!(package_subdir, 'package sub directory' )
  assert_exist!(build_root    , 'build root directory'   )
  
  $report << "Removing package [#{package_subdir}] from [#{build_root}]\n"

  files.each do |file|
  
    begin # exception block

      dir_name, file_name = File.split(file)
      
      src_file = package_subdir + '/' + dir_name + '/' + file_name
      dst_file = build_root     + '/' + dir_name + '/' + file_name
      
      assert_exist!(src_file, 'source file')

      if ( File.exist?(dst_file) ) then
        if ( options[:rename_instead] ) then
          $report << "Renaming [#{dst_file}] to [#{dst_file}.removed].\n"
          File.rename(dst_file, dst_file + ".removed" )
        else  
          $report << "Removing [#{dst_file}].\n"
          File.delete(dst_file)
        end
      else
        $report << "[#{dst_file}] does not exist in the build tree.\n"
      end
      
    rescue Exception => e
      $report << "Error removig #{file}: #{e}" << ENDL
    end
  end
end


# -- TEST -- 
#remove_package_files_from_build_tree(["include/CGAL/Straight_skeleton_2.h"],TEST_PKG_DIR, TEST_BUILD_ROOT, :rename_instead => false )
#remove_package_files_from_build_tree(list_package_files(TEST_PKG_DIR,false),TEST_PKG_DIR, TEST_BUILD_ROOT, :rename_instead => false )
# -- TEST -- 


