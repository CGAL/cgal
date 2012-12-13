class Logger

  attr_accessor(:device)
  
  def initialize(dev)
    @device = dev
  end
  
  def <<(msg)
    @device << msg
  end
end

ENDL = "\n"

$debout = Logger.new($stdout) 
$report = Logger.new($stdout)

def assert_exist!(path, what)
  raise "CGAL ERROR: Could not found #{what}: #{path}" unless File.exist?(path)  
end  

TEST_PKG_ROOT   = '../..'
TEST_BUILD_ROOT = '../../Build'
TEST_PKG_DIR    = '../../Straight_skeleton_2'

