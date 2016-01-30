macro any(comp)
    if !isa(comp, Expr) || (comp.head != :comprehension)
        error("@any may only precede a comprehension")
    end
    loop = Expr(:for, esc(comp.args[2]), 
                    quote
                        if $(esc(comp.args[1]))
                            result = true
                            break
                        end
                    end)
    quote
        result = false
        $loop
        result
    end
end

macro all(comp)
    if !isa(comp, Expr) || (comp.head != :comprehension)
        error("@all may only precede a comprehension")
    end
    loop = Expr(:for, esc(comp.args[2]), 
                    quote
                        if !$(esc(comp.args[1]))
                            result = false
                            break
                        end
                    end)
    quote
        result = true
        $loop
        result
    end
end

include("vec2Dutils.jl")
include("FSAutils.jl")