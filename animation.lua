local gasket = require 'gasket'

local reg = 2*math.sqrt(3) - 3

local anim = {
	gasket.radii_phase(0.5, 0.5, 0),
	gasket.radii_phase(reg, reg, 1/12),
	gasket.new(0,0,-0.5)
}

function set_frame(t)
	local k = math.floor(t)
	if k+1 == #anim then
		k = k - 1
	end
	local a = anim[k+1]
	local b = anim[k+2]
	local c = t - k
	local g = gasket.new((1-c)*a.fu+c*b.fu, (1-c)*a.ha+c*b.ha, (1-c)*a.fv+c*b.fv)
	g.apply_to_genome(frame:get_genome())
end

frame:get_genome():center(0,0)
frame:get_genome():pixels_per_unit(190)
frame:get_genome():width(400)
frame:get_genome():height(400)
frame:get_genome():background(0,0,0)

for k=1,120 do
	i = k - 1
	filename = 'gasket/frame'..string.format('%03d',i)..'.png'
	set_frame(2.0*i/120)
	frame:render(filename)
end