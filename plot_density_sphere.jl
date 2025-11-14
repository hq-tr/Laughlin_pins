include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/Density.jl")
using .FQH_states
using .ParticleDensity
using BenchmarkTools
using ArgMacros
using Plots; plotlyjs()

ENV["GKSwstype"] = "100" # Disconnect from Display

function main()
@inlinearguments begin
    @argumentrequired String fname "-f" "--filename"
end

state = readwf(fname)
No = length(state.basis[1])
R_max = sqrt(2*No)+0.2

println("Check norm = $(wfnorm(state))")

θ = LinRange(0,π,40)
ϕ = LinRange(0,2π, 80)

N₁ = length(θ)
N₂ = length(ϕ)

Θ = repeat(θ, inner=(1,N₂))
Φ = repeat(ϕ', inner=(N₁,1))

#den = zeros(size(Z))
@time den = get_density_sphere(state, Θ, Φ )

println("-----")
#den = @time get_density_disk(state, X, Y)

open("$(fname)_density_sphere.dat", "w") do f
    for i in 1:N₁, j in 1:N₂
        write(f, "$(Θ[i,j])\t$(Φ[i,j])\t$(den[i,j])\n")
    end
end

#p = plot(heatmap(z=density, aspect_ratio=:equal)) 

xx = sin.(Θ) .* cos.(Φ)
yy = sin.(Θ) .* sin.(Φ)
zz = cos.(Θ)

p = surface(xx, yy, zz, fill_z = den, size=(600,600))
savefig(p, "$(fname)_density.html")
end

main()