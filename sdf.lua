local complex = require 'complex'
local ii = complex { 0, 1 }

sdf = {}

sdf.new = function(a_, b_, c_, d_)
    local self = {
        a = a_,
        b = b_,
        c = c_,
        d = d_
    }

    function self.inside(z)
        local x, y = complex.get(z)
        return self.a*(x*x+y*y)+self.b*x+self.c*y+self.d < 0
    end

    function self.circ_rect_coll(p, w, h)
        local b = self.b/self.a
        local c = self.c/self.a
        local d = self.d/self.a
        local xc = -b/2
        local yc = -c/2
        local r = math.sqrt(xc*xc+yc*yc-d)
        local xp, yp = complex.get(p)
        local dx = math.abs(xc-xp)
        local dy = math.abs(yc-yp)
        if (dx > w/2 + r) or (dy > h/2 + r) then
            return false
        elseif (dx < w/2) or (dy < h/2) then
            return true
        else
            return (dx-w/2)*(dx-w/2)+(dy-h/2)*(dy-h/2) < r*r
        end
    end

    function self.rect_inside(p, w, h)
        local a = p - w/2 - ii*h/2
        local b = p + w/2 - ii*h/2
        local c = p - w/2 + ii*h/2
        local d = p + w/2 + ii*h/2
        if math.abs(self.a)>0 and self.a<0 and self.circ_rect_coll(p,w,h) then
            return false
        else
            return self.inside(a) and self.inside(b) and self.inside(c) and self.inside(d)
        end
    end

    function self.flip()
        return sdf.new(-self.a, -self.b, -self.c, -self.d)
    end

    return self
end

local function det(m)
    return m[1][1]*m[2][2]*m[3][3]+m[1][2]*m[2][3]*m[3][1]+m[1][3]*m[2][1]*m[3][2]-
        m[1][1]*m[2][3]*m[3][2]-m[1][2]*m[2][1]*m[3][3]-m[1][3]*m[2][2]*m[3][1]
end

sdf.from_pts = function(p, q, r)
    local xp, yp = complex.get(p)
    local np = xp*xp+yp*yp
    local xq, yq = complex.get(q)
    local nq = xq*xq+yq*yq
    local xr, yr = complex.get(r)
    local nr = xr*xr+yr*yr

    local a = det({
        {xp, yp, 1},
        {xq, yq, 1},
        {xr, yr, 1}
    })

    local b = -det({
        {np, yp, 1},
        {nq, yq, 1},
        {nr, yr, 1}
    })

    local c = det({
        {np, xp, 1},
        {nq, xq, 1},
        {nr, xr, 1}
    })

    local d = -det({
        {np, xp, yp},
        {nq, xq, yq},
        {nr, xr, yr}
    })

    return sdf.new(a,b,c,d)
end

return sdf
