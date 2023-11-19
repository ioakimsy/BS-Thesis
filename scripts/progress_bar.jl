using Pkg

#= Pkg.add("ProgressMeter")
Pkg.add("Alert") =#
using ProgressMeter
using Alert

@alert begin
    x,n = 1,1000
    p = Progress(n)
    for iter = 1:n
        x *= 2
        sleep(0.001)
        ProgressMeter.next!(p; showvalues = [("iter",iter), ("x", x)])
    end
end
