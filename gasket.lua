local complex = require 'complex'
local mobius = require 'mobius'
local sdf = require 'sdf'
local util = require 'util'

local EPS = 1e-10
local ii = complex { 0, 1 }
local sr3 = complex(math.sqrt(3))
local pi = 2*math.acos(0)

local gasket = {}

gasket.new = function(fu_, ha_, fv_, ar_, pc_, sc_, minimal_, corr_)
    local self = {
        fu = fu_,
        ha = ha_,
        fv = fv_,
        ar = (ar_ and ar_ or 16/9),
        pc = (pc_ and pc_ or complex(0)),
        sc = (sc_ and sc_ or 1),
        minimal = (minimal_ == true),
        corr = (corr_ and corr_ or 5)
    }

    local function adapt(p1, p2, p3, dive, rot)
        local ans = mobius.scaling(1)
        local height = 2/self.sc
        local width = height*self.ar
        local arr = {dive, dive.conjugate(rot), dive.conjugate(rot.inverse())}
        ans = mobius.scaling(1)
        local search = true
        while search do
            search = false
            for i=1,3 do
                local q1 = ans.compose(arr[i]).apply(p1)
                local q2 = ans.compose(arr[i]).apply(p2)
                local q3 = ans.compose(arr[i]).apply(p3)
                local shape = sdf.from_pts(q1,q2,q3)
                if shape.rect_inside(self.pc, width, height) then
                    search = true
                    ans = ans.compose(arr[i])
                    break
                end
            end
        end
        return ans
    end

    local u = math.cosh(self.ha)*complex.exp(2*self.fu*pi*ii)
    local v = math.sinh(self.ha)*complex.exp(2*self.fv*pi*ii)

    local s = mobius.new(u,v,v:conjugate(),u:conjugate())
    local st = mobius.scaling(self.sc).compose(mobius.translation(-self.pc)).inverse()

    self.mobius_list = {}

    local tr = mobius.new(1,0,-2*ii,1).conjugate(s)
    local rot = mobius.pts_to_pts(complex(0),1,-1,complex(1),-1,0).conjugate(s)

    local p1 = s.inverse().apply(-1)
    local p2 = s.inverse().apply(0)
    local p3 = s.inverse().apply(1)

    local shape = sdf.from_pts(p1,p2,p3)
    local height = 2/self.sc
    local width = height*self.ar

    if self.minimal then
        table.insert(self.mobius_list, tr.conjugate(st))
        table.insert(self.mobius_list, rot.compose(tr.inverse()).conjugate(st))
    elseif shape.rect_inside(self.pc, width, height) then
        local m = adapt(p1,p2,p3,tr,rot)
        tr = tr.conjugate(m.inverse())
        rot = rot.conjugate(m.inverse())
        table.insert(self.mobius_list, tr.conjugate(st))
        table.insert(self.mobius_list, tr.conjugate(rot).conjugate(st))
        table.insert(self.mobius_list, tr.conjugate(rot.inverse()).conjugate(st))
    elseif shape.flip().rect_inside(self.pc, width, height) then
        local m = adapt(p1,p3,p2,tr.inverse(),rot)
        tr = tr.conjugate(m.inverse())
        rot = rot.conjugate(m.inverse())
        table.insert(self.mobius_list, tr.inverse().conjugate(st))
        table.insert(self.mobius_list, tr.inverse().conjugate(rot).conjugate(st))
        table.insert(self.mobius_list, tr.inverse().conjugate(rot.inverse()).conjugate(st))
    else
        table.insert(self.mobius_list, tr.conjugate(st))
        table.insert(self.mobius_list, tr.conjugate(rot).conjugate(st))
        table.insert(self.mobius_list, tr.conjugate(rot.inverse()).conjugate(st))
        table.insert(self.mobius_list, tr.inverse().conjugate(st))
        table.insert(self.mobius_list, tr.inverse().conjugate(rot).conjugate(st))
        table.insert(self.mobius_list, tr.inverse().conjugate(rot.inverse()).conjugate(st))
    end

    local function clear_genome(g, sz)
        for i = g:num_xforms(),2,-1 do
            g:del_xform(i)
        end
        local xf = g:get_xform(1)
        for i, _ in ipairs(VARIATIONS) do
            xf:var(i, 0)
        end
        xf:density(0.5)
        g:chaos(1,{1})
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
        if not self.minimal and #self.mobius_list == 6 then
            local arr = {}
            for i=1,3 do
                arr[i] = {1,1,1,0,0,0}
            end
            for i=4,6 do
                arr[i] = {0,0,0,1,1,1}
            end
            for i=1,6 do
                arr[i][i] = self.corr
            end
            for i=1,6 do
                g:chaos(i,arr[i])
            end
        end
        if not self.minimal and #self.mobius_list == 3 then
            local arr = {}
            for i=1,3 do
                arr[i] = {1,1,1}
            end
            for i=1,3 do
                arr[i][i] = self.corr
            end
            for i=1,3 do
                g:chaos(i,arr[i])
            end
        end
    end

    return self
end

gasket.radii_phase = function(r1, r2, f, flip, ar, pc, sc, minimal, corr)
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

    if flip then
        return gasket.new(-fu, ha, -fv, ar, pc, sc, minimal, corr)
    else
	    return gasket.new(fu, ha, fv, ar, pc, sc, minimal, corr)
    end
end

gasket.params = function(a, b, f, flip, ar, pc, sc, minimal, corr)
    local reg = 2*math.sqrt(3) - 3
    local r1 = reg*(1-a)+a
    local top = math.min(r1,1-r1)
    local bot = 4*(1-r1)*r1/((1+r1)*(1+r1))
    local r2 = top*(1-b)+bot*b
    return gasket.radii_phase(r1, r2, f, flip, ar, pc, sc, minimal, corr)
end

return gasket
