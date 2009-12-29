----------------------------------------------------------------------
-- CGAL Skeleton and offset ipelet description
----------------------------------------------------------------------

label = "Skeleton and offset"

about = [[
This ipelet is part of the CGAL_ipelet package. See www.cgal.org.
]]

-- this variable will store the C++ ipelet when it has been loaded
ipelet = false

function run(model, num)
  if not ipelet then ipelet = assert(ipe.Ipelet(dllname)) end
  model:runIpelet(methods[num].label, ipelet, num)
end

methods = {
  { label="Interior skeleton" },
  { label="Exterior skeleton" },
  { label="Interior offset" },
  { label="Exterior offset" },
  { label="Interior offsets" },
  { label="Exterior offsets" },
  { label="Help" },
}

----------------------------------------------------------------------
