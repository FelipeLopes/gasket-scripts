local complex = require 'complex'
local mobius = require 'mobius'
local util = require 'util'

local EPS = 1e-10
local ii = complex { 0, 1 }
local sr3 = complex(math.sqrt(3))
local pi = 2*math.acos(0)

local gasket = {}

gasket.new = function(fu_, ha_, fv_)
    local self = {
        fu = fu_,
        ha = ha_,
        fv = fv_
    }

    local u = math.cosh(self.ha)*complex.exp(2*self.fu*pi*ii)
    local v = math.sinh(self.ha)*complex.exp(2*self.fv*pi*ii)

    local tr = mobius.translation(2)
    local rot = mobius.scaling(0.5*(-1+sr3*ii)).conjugate(mobius.new(1,sr3,-1,sr3))
    local s = mobius.new(0,ii,1,0).compose(mobius.new(u,v,v:conjugate(),u:conjugate()))

    self.mobius_list = {}
    
    table.insert(self.mobius_list, tr.conjugate(s))
    table.insert(self.mobius_list, rot.compose(tr.inverse()).conjugate(s))

    local function clear_genome(g, sz)
        for i = g:num_xforms(),2,-1 do
            g:del_xform(i)
        end
        local xf = g:get_xform(1)
        for i, _ in ipairs(VARIATIONS) do
            xf:var(i, 0)
        end
        xf:density(0.5)
        for i=2,sz do
            local xf = g:add_xform()
            xf:var(LINEAR, 0)
        end
    end

    local function set_xform(xf, data)
        xf:var(data.var, 1)
        xf:o(data.pre.o[1], data.pre.o[2])
        xf:x(data.pre.x[1], data.pre.x[2])
        xf:y(data.pre.y[1], data.pre.y[2])
        xf:op(data.post.o[1], data.post.o[2])
        xf:xp(data.post.x[1], data.post.x[2])
        xf:yp(data.post.y[1], data.post.y[2])
    end

    function self.apply_to_genome(g)
        clear_genome(g, #self.mobius_list)
        for i=1,#self.mobius_list do
          set_xform(g:get_xform(i), self.mobius_list[i].xform())
        end
    end

    return self
end

gasket.radii_phase = function(r1, r2, f)
	if r2 > r1 and math.abs(r1-r2)>EPS then
		error('first radius parameter should be greater than or equal to the second')
	elseif r1 + r2 > 1 and math.abs(r1+r2-1)>EPS then
		error('radii sum cannot be larger than 1')
	elseif r1<0 or math.abs(r1)<EPS or r2<0 or math.abs(r2)<EPS then
		error('all radii should be positive')
	end
	local a = -1
	local b = 1/r1
	local c = 1/r2
	local s1 = a + b + c
	local s2 = math.max(a*b + b*c + c*a, 0)
	local d = s1 -2*math.sqrt(s2)
	if d<0 or math.abs(d)<EPS then
		d = s1 + 2*math.sqrt(s2)
	end
	if d < c and math.abs(d-c)>EPS then
		error('radii given are not the largest circles')
	end

	local l1 = r1 + r2
	local l2 = 1 - r1
	local l3 = 1 - r2
	local cosx = (l2*l2+l3*l3-l1*l1)/(2*l2*l3)
	local sinx = math.sqrt(math.max(1-cosx*cosx, 0))

	local p1 = complex.exp(f*2*pi*ii)
	local p2 = p1*(cosx+sinx*ii)

	local v1 = (1-r1)*p1
	local v2 = (1-r2)*p2
	local v3 = (r2*v1+r1*v2)/(r1+r2)

	local m = mobius.pts_to_pts(p1,p2,v3,complex(1),-1,0)

	local _, fu = complex.get(complex.ln(m.a)/(2*pi))
	local ha = util.asinh(complex.abs(m.b))
	local fv = 0
	if (math.abs(ha) > EPS) then
		_, fv = complex.get(complex.ln(m.b)/(2*pi))
	end

	return gasket.new(fu, ha, fv)
end

return gasket
