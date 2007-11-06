require 'Pathname'
require 'FileUtils'

load 'common_impl.rb'

def default_mirror_op
  return RUBY_PLATFORM =~ /mswin32|cygwin|mingw|bccwin/ ? :hardlink : :symlink
end
  
#
# "Mirror" the listed files from a given package folder into a given target folder.
# Entries in the file list are pathnames relative to the package folder.
# The directory structure under the package folder is replicated in the build folder,
# creating new folders as neccesary
#
# "Mirroring" a file consist on creating within the build folder a view to a file located into the package folder.
# In true POSIX platforms this is done via a symlink.
# In Windows platforms hardlinks are used instead. Hardlinks are recrated whenever the source file is newer, this
# ensures that broken hardlinks are syncronized again.
# 
# The mirroring operation is automatically selected according to the platform, but can be overrided if needed.
#
# By setting mirror_op = :copy, this function can also be used to create
# a clean copy (instead of a mirror) of the package files inside a build folder.
# In this case, if the destination file must not exist.
#
def mirror_pkg_files(files,
                     package_dir,
                     build_dir, 
                     mirror_op = default_mirror_op
                    )

  package_dir = as_pathname(package_dir)
  build_dir   = as_pathname(build_dir)

  assert_exist!(package_dir, 'package directory' )
  assert_exist!(build_dir  , 'build directory'   )
  
  $report << "Mirroring package [#{package_dir}] into [#{build_dir}]\n"

  # Keep a local hash of subfolders to avoid accessing the filesystem redudantly
  subdir_exist = {}

  files.each do |file|
  
    begin # exception block

      dir_name, file_name = File.split(file)
      
      src_file = package_dir + dir_name + file_name
      dst_file = build_dir   + dir_name + file_name
      
      assert_exist!(src_file, 'source file')

      # Remove existing hardlink if outdated
      if ( mirror_op == :hardlink ) then
        if ( File.exist?(dst_file) ) then
          if ( File.mtime(src_file) != File.mtime(dst_file) ) then
            $report << "Hardlink [#{dst_file}] is outdated. Removing it.\n"
            File.delete(dst_file)
          end
        end
      end

      # Already mirror files are never mirrored agaiin.
      unless File.exist?(dst_file) then
      
        #
        # Replicate directory structure as needed
        #
        Pathname.new(dir_name).descend do |local_subdir|
        
          dst_subdir = build_dir + local_subdir 
          
          unless subdir_exist[dst_subdir]
            unless FileTest.exist?(dst_subdir)
            
              $report << "Creating build subdir [#{dst_subdir}]\n"
              
              Dir.mkdir(dst_subdir)
              
            end  
            subdir_exist[dst_subdir] = true
          end  
        end
      
        #
        # Mirror file
        #
        case ( mirror_op )
          when :symlink
            $report << "Creating symlink from [#{src_file}] to [#{dst_file}]\n" 
            File.symlink(src_file, dst_file)
          when :hardlink
            $report << "Creating hardlink from [#{src_file}] to [#{dst_file}]\n"
            File.link(src_file, dst_file)
          when :copy
            $report << "Creating a clean copy of [#{src_file}] as [#{dst_file}]\n"
            File.copy(src_file, dst_file)
        end
        
      end
      
    rescue Exception => e
      $report << "Error installing #{file}: #{e}" << ENDL
    end
  end
end


#mirror_pkg_files(["include/CGAL/Straight_skeleton_2.h"],$test_pkg_dir, $test_build_dir)
#mirror_pkg_files(list_pkg_files(pkg_dont_submit_list($test_pkg_dir), $test_pkg_dir),$test_pkg_dir, $test_build_dir)


