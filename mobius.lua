local complex = require 'complex'
local mobius = {}
local EPS = 1e-10
local ii = complex { 0, 1 }

mobius.new = function(a_, b_, c_, d_)
    local sdet = complex.sqrt(a_*d_ - b_*c_)
    local self = {
        a = a_/sdet,
        b = b_/sdet,
        c = c_/sdet,
        d = d_/sdet
    }

    function self.apply(z)
        return (self.a*z+self.b)/(self.c*z+self.d)
    end

    function self.inverse()
        return mobius.new(self.d, -self.b, -self.c, self.a)
    end

    function self.compose(m)
        local a = self.a
        local b = self.b
        local c = self.c
        local d = self.d
        return mobius.new(a*m.a+b*m.c, a*m.b+b*m.d, c*m.a+d*m.c, c*m.b+d*m.d)
    end

    function self.conjugate(s)
        return s.inverse().compose(self).compose(s)
    end

    function self.xform()
        local a = self.a
        local b = self.b
        local c = self.c
        local d = self.d
        if (complex.abs(c) > EPS) then
            local pre = mobius.translation(d/c)
            local post = mobius.translation(a/c).compose(mobius.scaling(-1/(c*c)))
            return {
                var = 3,
                pre = {
                    o = {complex.get(pre.apply(0):conjugate())},
                    x = {complex.get(pre.apply(1):conjugate())},
                    y = {complex.get(pre.apply(ii):conjugate())}
                },
                post = {
                    o = {complex.get(post.apply(0))},
                    x = {complex.get(post.apply(1))},
                    y = {complex.get(post.apply(ii))}
                }
            }
        else
            return {
                var = 1,
                pre = {
                    o = {complex.get(self.apply(0))},
                    x = {complex.get(self.apply(1))},
                    y = {complex.get(self.apply(ii))}
                },
                post = {
                    o = {0, 0},
                    x = {1, 0},
                    y = {0, 1}
                }
            }
        end
    end

    return self
end

mobius.scaling = function(a)
    return mobius.new(complex(a),0,0,1)
end

mobius.translation = function(b)
    return mobius.new(1,complex(b),0,1)
end

mobius.pts = function(p, q, r)
    return mobius.new(q-r,-p*(q-r),q-p,-r*(q-p))
end

mobius.pts_to_pts = function(p1, q1, r1, p2, q2, r2)
    return mobius.pts(p2, q2, r2).inverse().compose(mobius.pts(p1, q1, r1))
end

return mobius
