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
        @argumentdefault Int 2 pin_size "-k" "--pin_size"
        @argumentdefault String "" intname "-i" "--interaction-file"
        #@argumentdefault String "" pinname "-p" "--pin-potential-file"
        @argumentoptional Int n_el "-e" "--n_el"
        @argumentoptional Int n_orb "-o" "--n_orb"
        @argumentflag special "-s" "--special"
        @argumentflag positive_pin "-p" "--positive"
        @argumentdefault String "Arpack" EDmethod "--method"
        @argumentoptional Float64 L_z "-Lz" "--Lz"
        #@argumentdefault Float64 3.141592653589793 angle_multiplier "--angle-multiplier" 
    end

    println("============================================================")
    println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
    println("      with only a potential pin at the north pole")
    println()
    println("In this 'special' version, the north pin consists of two parts:")
    println("A negative part including the first k orbitals,")
    println("and a positive part including the (k+1)-th orbital.")
    println("============================================================")

    if L_z == nothing && n_el !=nothing
        Lzs = collect(-n_el:1:n_el)
    else
        Lzs = [L_z]
    end

    println(Lzs)
    for Lz in Lzs
        # Reading basis input
        if n_el != nothing && n_orb != nothing
            if Lz == nothing
                println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
                basis = fullhilbertspace(n_el,n_orb)
                outname = @sprintf "%ie_%io" n_el n_orb
            else
                println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals in the Lz=$(Lz) sector.")
                basis = fullhilbertspace(n_el,n_orb,Lz)
                outname = @sprintf "%ie_%io" n_el n_orb
            end
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
            @time V1_matrix = two_body(N_o, basis, v_list, c_list)#; quiet=true)
        else
            @time V1_matrix = spzeros(Complex{Float64},(d,d))
        end


        println("Constructing the pinning potentials")
        if special
            if positive_pin
                @time pin_north_matrix = 2*sphere_wide_bump_matrix(basis, pin_size) - sphere_wide_bump_matrix(basis,pin_size+1) 
                textappend = "_positive_special"
            else
                @time pin_north_matrix = -2*sphere_wide_bump_matrix(basis, pin_size) + sphere_wide_bump_matrix(basis,pin_size+1) 
                textappend = "_special"
            end
        else
            if positive_pin
                @time pin_north_matrix = sphere_wide_bump_matrix(basis, pin_size)
                textappend = "_positive"
            else
                @time pin_north_matrix = -1*sphere_wide_bump_matrix(basis, pin_size)
                textappend = ""
            end
        end

        #λ_range = vcat(collect(0:0.05:0.95),[0.99,0.999])
        λ_range = vcat(collect(0:0.05:1),collect(1.5:0.5:5),[10,20,50,100,1000])
        #λ_range = [0.99,0.999]
        for λ in λ_range
            try
                println("===== WORKING ON λ = $λ =====")
                H_matrix = V1_matrix + λ*pin_north_matrix + 1000*sparse(I, d, d)

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

                if !isdir("wide_pin_k_$(pin_size)$(textappend)_north") mkdir("wide_pin_k_$(pin_size)$(textappend)_north") end

                dirname = "wide_pin_k_$(pin_size)$(textappend)_north/$(outname)_$(intname)_$(pinappendname)_out"

                if !isdir(dirname) mkdir(dirname) end

                if EDmethod=="Krylov"
                    extension="_krylov"
                else
                    extension=""
                end

                if Lz == nothing
                    Lztext = ""
                else
                    Lztext = "_L_$(Lz)"
                end

                open("$(dirname)/eigen$(extension)$(Lztext).txt","w+") do f
                    for i in 1:k
                        if EDmethod=="Krylov"
                            gs_coef = ϕ[i]  # Use this if using KrylovKit
                        else
                            gs_coef = ϕ[:,i] # Use this if using Arpack
                        end
                        write(f,"$(abs(ϵ[i])-1000 )\t")
                        ground_state = FQH_state(basis, gs_coef)
                        printwf(ground_state;fname="$(dirname)/g$(extension)$(Lztext)_$(i-1)")
                        LZ = get_Lz_sphere(ground_state)
                        write(f,"$(LZ)\n")
                    end
                end
            catch
            finally
                continue
            end
        end
    end # End for Lz in Lzs
end

@time main()