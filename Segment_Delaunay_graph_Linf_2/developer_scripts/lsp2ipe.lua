-- lsp2ipe.lua
-- convert segments in lsp format from stdin
-- to an ipe file that is given in the command line;

function printusage()
  io.stderr:write("Usage: ipescript lsp2ipe <file_to_save>")
  io.stderr:write(" [[offsetx offsety] scale]\n")
end

-- print('#argv = ', #argv)
if #argv < 1 then
  printusage()
  return
end
fname = argv[1]

offsetx = 0
offsety = 0
scale   = 1

if #argv > 1 then
  if #argv == 3 or #argv == 4 then
    offsetx = tonumber(argv[2])
    offsety = tonumber(argv[3])
    if #argv == 4 then
      scale = tonumber(argv[4])
    end
  else
    printusage()
  end
end

if offsetx == nil or offsety == nil or scale == nil then
  io.stderr:write("Error: bad number input in command line\n")
  return
end


-- find basic.isy ipe style sheet

-- list of possible locations
basiclist = {
  "/opt/local/share/ipe/7.1.4/styles/basic.isy",
  "/usr/share/ipe/7.1.2/styles/basic.isy",
  "/usr/share/ipe/7.1.1/styles/basic.isy",
  "/opt/local/share/ipe/7.1.2/styles/basic.isy"
}

function file_exists(name)
   local f = io.open(name,"r")
   if f ~= nil then io.close(f) return true else return false end
end

basicfilename = nil

for i, filename in ipairs(basiclist) do
  --print("trying ", filename)
  if file_exists(filename) then
    basicfilename = filename
    break
  end
end

if basicfilename == nil then
  io.stderr:write("Error: missing basic.isy\n")
  return
end

doc = ipe.Document()
sheet = ipe.Sheet(basicfilename)
sheets = ipe.Sheets()
sheets:insert(1, sheet)
doc:replaceSheets(sheets)


attr = {}
for l in io.lines() do
  print(l)
  t = {}
  for w in l:gmatch("[^\t ]+") do
    table.insert(t, w)
  end
  print("#t=", #t)
  if #t ~= 0 then
    print(t[1])
    if t[1] == 'p' then
      if #t ~= 3 then
        io.stderr:write("Error: bad point line, t=", #t, "\n")
        return
      end
      --print("point found")
      px = tonumber(t[2])
      py = tonumber(t[3])
      if px == nil or py == nil then
        io.stderr:write("Error: bad number input\n")
        return
      end
      fpx = offsetx + (px * scale)
      fpy = offsety + (py * scale)
      v = ipe.Vector(fpx, fpy)
      obj = ipe.Reference(attr, "mark/disk(sx)", v)
    elseif t[1] == 's' then
      if #t ~= 5 then
        io.stderr:write("Error: bad segment line, t=", #t, "\n")
        return
      end
      px = tonumber(t[2])
      py = tonumber(t[3])
      qx = tonumber(t[4])
      qy = tonumber(t[5])
      if px == nil or py == nil or qx == nil or qy == nil then
        io.stderr:write("Error: bad number input\n")
        return
      end
      fpx = offsetx + (px * scale)
      fpy = offsety + (py * scale)
      fqx = offsetx + (qx * scale)
      fqy = offsety + (qy * scale)

      seg = {}
      seg["type"] = "segment"
      seg[1] = ipe.Vector(fpx, fpy)
      seg[2] = ipe.Vector(fqx, fqy)
      subpath = {}
      subpath["type"] = "curve"
      subpath["closed"] = false
      subpath[1] = seg
      shape = { subpath }
      obj = ipe.Path(attr, shape, false)
    else
      -- here we have neither 'p' nor 's'
      io.stderr:write("Error: bad line type (neither p nor s):\n")
      io.stderr:write(l)
      return
    end
    doc[1]:insert(nil, obj, nil, "alpha")
  end
end

doc:save(fname, xml, nil)
