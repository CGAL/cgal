require 'Pathname'

class Logger

  attr_accessor(:device)
  
  def initialize(dev)
    @device = dev
  end
  
  def <<(msg)
    @device << msg
  end
end

ENDL = "\n" unless defined?(ENDL)

$debout = Logger.new($stdout) 
$report = Logger.new($stdout)

def as_pathname(path)
  return path.kind_of?(Pathname) ? path : Pathname.new(path)
end

def assert_exist!(path, what)
  raise "CGAL ERROR: Could not found #{what}: #{path}" unless FileTest.exist?(path)  
end  

$test_pkg_dir   = '../../Straight_skeleton_2'
$test_build_dir = '../../Build'


