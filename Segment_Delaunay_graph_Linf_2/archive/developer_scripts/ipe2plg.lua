-- ipe2cin.plg
-- convert points, segments, and closed polygons in ipe file to plg format,
-- only for first page and its first view of ipe file
-- and only for visible layers in the first view

if #argv ~= 1 then
  io.stderr:write("Usage: ipescript ipe2cin <file>\n")
  return
end
fname = argv[1]

doc = assert(ipe.Document(fname))

-- io.write("# number of pages = ", #doc, "\n")
-- print("# of pages =", #doc)

p = assert(doc[1])

-- print("# of views=", p:countViews())

-- print("# of objects in p1 =", #p)

function process_object(obj)
  -- print(obj:matrix())
  m = obj:matrix()
  is_id = m:isIdentity()
  if obj:type() == "reference" then
    if (obj:get("markshape")):find("mark") then
      if is_id then
        vec = obj:position()
      else
        vec = m * obj:position()
      end
      print("1")
      print(vec.x, vec.y)
    end
  end
  if obj:type() == "path" then
    shape = obj:shape()
    for i, subpath in pairs(shape) do
      if subpath["type"] == "curve" then
        opencurve = false
        allsegments = true
        manycomponents = false
        num_components = #subpath
        --print("#subpath components =", num_components)
        --print("#subpath is closed =", subpath["closed"])
        pointfirst = nil
        pointlast = nil
        points = {}
        for i, seg in ipairs(subpath) do
          atleastoneseg = false
          if seg["type"] == "segment" then
            atleastoneseg = true
            if is_id then
              endp1 = seg[1]
              endp2 = seg[2]
            else
              endp1 = m * seg[1]
              endp2 = m * seg[2]
            end
          else
            -- ignore non-segment component and set false variable
            allsegments = false
            break
          end
          table.insert(points, endp1)
          if i == 1 then
            pointfirst = endp1
            --print("pointfirst =', pointfirst)
          else
            -- i > 1
            manycomponents = true
          end
          if i == num_components then
            pointlast = endp2
            if subpath["closed"] then
              if pointfirst ~= pointlast then
                table.insert(points, pointlast)
              end
            else
              -- here, subpath is not closed
              if pointfirst ~= pointlast then
                -- problematic curve, ignore it
                opencurve = true
                break;
              end
            end
          end
        end -- of for i, seg in ipairs(subpath)
        if allsegments then
          if subpath["closed"] or
             ((not subpath["closed"]) and (not opencurve)) then
            print(#points)
            for i, p in ipairs(points) do
              print(p.x, p.y)
            end
          end
        end -- of all segments case
      end -- of curve case
    end -- of for subpath
    --print(obj:xml())
  end
end

for i, obj, sel, layer in p:objects() do
  --print(i, obj, sel, layer)
  if p:visible(1, layer) then
    process_object(obj)
  end
end
