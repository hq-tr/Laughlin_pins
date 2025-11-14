include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("/home/trung/_qhe-julia/Potentials.jl")
include("/home/trung/two_body_potential/v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using .HilbertSpaceGenerator
using .Potentials
using LinearAlgebra
using Arpack
using KrylovKit
using SparseArrays
using ArgMacros
using BenchmarkTools
using Printf

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @argumentdefault String "" fname "-f" "--filename"
        @argumentdefault Int 5 k "-n" "--nev"
        @argumentdefault Int 1 pin_size "-k" "--pin_size"
        @argumentdefault String "" intname "-i" "--interaction-file"
        #@argumentdefault String "" pinname "-p" "--pin-potential-file"
        @argumentoptional Int n_el "-e" "--n_el"
        @argumentoptional Int n_orb "-o" "--n_orb"
        @argumentflag special "-s" "--special"
        @argumentflag positive_pin "-p" "--positive"
        @argumentdefault String "Arpack" EDmethod "--method"
        @argumentdefault Float64 0. λ_start "--lambda_start"
        @argumentdefault Float64 1.5 λ_stop  "--lambda_stop"
        @argumentdefault Float64 0.1 λ_step "--lambda_step"
        #@argumentoptional Float64 L_z "-Lz" "--Lz"
        #@argumentoptional Float64 Lz_cutoff "--Lz-cutoff"
        #@argumentdefault Float64 3.141592653589793 angle_multiplier "--angle-multiplier" 
    end

    println("============================================================")
    println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
    println("               with four potential pins")
    println()
    println("In this version, there are four pins arranged into vertices of a regular tetrahedron.")
    println("They are identical in magnitude and size (shape)")
    println()
    println("The pins a negative by default, and positive if the 'positive' tag is included.")
    println("============================================================")

    if pin_size != 1
        print("Sorry, the code only works with k=1 for now.")
        return
    end

    # Reading basis input
    if n_el != nothing && n_orb != nothing
        println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
        basis = fullhilbertspace(n_el,n_orb)
        outname = @sprintf "%ie_%io" n_el n_orb
    elseif length(fname) > 0
        println("Reading basis vectors from [$(fname)]")
        
        basis = readwf(fname).basis

        outname = fname
    else
        println()
        println("WARNING: No input or incomplete input was specified. The program will now terminating.")
        println("Run the program with '-h' or '--help' tag to view possible arguments.")
        println()
        return
    end

    N_o = length(basis[1])
    d   = length(basis) # Dimension

    # Reading two-body interaction input
    v_list = Int32[]
    c_list = Float64[]


    if intname != "none"
        if length(intname) == 0
            println("Input m for Vₘ and the corresponding coefficient. ")
            println("Each pp term takes one line, with two numbers separated by a space.")
            println("Put a 0 to end")
            reading = true
            while reading
                data = readline()
                if data == "0"
                    reading = false
                else
                    try
                        pp = split(data)
                        push!(v_list,parse(Int32, pp[1]))
                        push!(c_list,parse(Float64,pp[2]))
                    catch
                        println("Invalid input. Try again or input 0 to end.")
                    end
                end
            end
        else
            println("Reading interaction from $(intname).")
            if isfile(intname)
                open(intname) do f
                    for line in map(s->split(s),readlines(f))
                        append!(v_list,parse(Int32,line[1]))
                        append!(c_list,parse(Float64,line[2]))
                    end
                end
            else
                print("Interaction file '$(intname)' not found. Terminating.")
                return false
            end
        end
    end




    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")

    if length(v_list) > 0
        @time V1_matrix = two_body(N_o, basis, v_list, c_list;verbose=true)
    else
        @time V1_matrix = spzeros(Complex{Float64},(d,d))
    end


    println("Constructing the pinning potentials")
    # three other vertices of the regular tetrahedron with one vertices at the north pole:
    # (2√2/3,0,-1/3), (-√2/3,√2/√3, -1/3), (-√2/3,-√2/√3,-1/3)
    # The θ angle is π/2 + arctan(1/2√2)
    θ₁ = π/2 + atan(1/√(8))
    ϕ₁ = 0.
    ϕ₂ = 2π/3
    ϕ₃ = 4π/3
    if positive_pin
        @time pin_matrix = sphere_wide_bump_matrix(basis, pin_size) + sphere_bump_matrix(basis,[θ₁,θ₁,θ₁],[ϕ₁,ϕ₂,ϕ₃];verbose=true)
        textappend = "_4pins_positive"
    else
        @time pin_matrix = -1*(sphere_wide_bump_matrix(basis, pin_size) + sphere_bump_matrix(basis,[θ₁,θ₁,θ₁],[ϕ₁,ϕ₂,ϕ₃];verbose=true))
        textappend = "_4pins"
    end

    if !isdir("wide_pin_k_$(pin_size)$(textappend)") mkdir("wide_pin_k_$(pin_size)$(textappend)") end

    #λ_range = vcat(collect(0:0.05:0.95),[0.99,0.999])
    #λ_range = vcat(collect(0:0.05:1),collect(1.5:0.5:5),[10,20,50,100,1000])
    #λ_range = [0.99,0.999]
    #λ_range = collect(0:1.:20.)
    #λ_range = collect(0:0.1:1.5)
    λ_range = λ_start:λ_step:λ_stop

    Threads.@threads for λ in λ_range
        try
            println("===== WORKING ON λ = $λ =====")
            H_matrix = V1_matrix + λ*pin_matrix + 1000*sparse(I, d, d)

            pinappendname = "$λ"

            #println("Hamiltonian matrix = ")
            #display(H_matrix)

            println("--------")

            if EDmethod == "Krylov"
                println("Diagonalizing with KrylovKit:")

                @time ϵ, ϕ, info = eigsolve(H_matrix, k,:SR)
            else
                if EDmethod != "Arpack"
                    println("Specified ED method not recognized. ARPACK will be used.")
                end

                println("Diagonalizing with ARPACK")

                @time ϵ, ϕ = eigs(H_matrix, nev=k,which=:SM)
            end

            #display(ϕ)

            # ====================== SAVE GROUND STATE =======================
            #println("Eigenvalues = ")
            #display(abs.(ϵ.-1000))

            println("--------")

            dirname = "wide_pin_k_$(pin_size)$(textappend)/$(outname)_$(intname)_$(pinappendname)_out"

            if !isdir(dirname) mkdir(dirname) end

            if EDmethod=="Krylov"
                extension="_krylov"
            else
                extension=""
            end

            open("$(dirname)/eigen$(extension).txt","w+") do f
                for i in 1:k
                    if EDmethod=="Krylov"
                        gs_coef = ϕ[i]  # Use this if using KrylovKit
                    else
                        gs_coef = ϕ[:,i] # Use this if using Arpack
                    end
                    write(f,"$(abs(ϵ[i])-1000 )\t")
                    ground_state = FQH_state(basis, gs_coef)
                    printwf(ground_state;fname="$(dirname)/g$(extension)_$(i-1)")
                    LZ = get_Lz_sphere(ground_state)
                    write(f,"$(LZ)\n")
                end
            end
        catch
        finally
            continue
        end
    end

end

@time main()