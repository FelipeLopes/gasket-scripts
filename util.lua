local util = {}

util.asinh = function(x)
    return math.log(x + math.sqrt(x*x+1))
end

util.acosh = function(x)
    return math.log(x + math.sqrt(x*x-1))
end

util.tprint = function(tbl, indent)
    if not indent then indent = 0 end
    for k, v in pairs(tbl) do
        formatting = string.rep("  ", indent) .. k .. ": "
        if type(v) == "table" then
            print(formatting)
            util.tprint(v, indent+1)
        elseif type(v) == 'boolean' then
            print(formatting .. tostring(v))
        else
            print(formatting .. v)
        end
    end
end

return util
