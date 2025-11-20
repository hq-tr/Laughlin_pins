# Pinning potentials on the sphere

## Description
The routines in this directory perform exact diagonalization (ED) on the many-body Hilbert space in the lowest Landau level (LLL) on the sphere for the Hamiltonian of the form

$$H = V^{2bdy} + \lambda V_{\text{pins}}$$

where $$V^{2bdy}$$ is any combination of two-body Haldane pseudo potentials and $$V_{\text{pins}}$$ varies depending on the version of the script. $$\lambda$$ is a tuneable parameter, and ED is carried out concurrently (in embarrassingly parallel processes) for different values of $$\lambda$$ within the specified range. 

The current available scripts are:
* `sphere_FQHE_widebump_onlynorth_parallel.jl` -- Only one pin at the north pole. 
* `sphere_FQHE_widebump_samesign_parallel.jl` -- Two pins, one at the north pole and one at the south pole. The two pins are identical (same sign and same magnitude)
* `sphere_FQHE_widebump_4pins_parallel.jl` -- Four pins arranged at the vertices of a regular tetrahedron, with one pin at the north pole. All four pins are identical.

In every version and with the default setting, each pin attracts electrons within $k$ LLL orbital around it. $k$ is a tuneable parameter. They can also be made into positive pins, which repel electrons (attract holes or quasiholes) by adding the `--postive` flag.

Since the Hilbert space is the LLL, the Hamiltonian matrix is sparse. The default ED method is the Arnoldi method (a.k.a. implicitly restarted Lanczos). There is an option to use the method provided by KrylovKit.jl, which is recommend as it often seems slightly faster than ARPACK. This can be done by the tag `--method Krylov` (note that "Krylov" is case-sensitive). ARPACK is left as the default method because of legacy reasons.

## Command
These routines requires the [QHE_Julia](https://github.com/hq-tr/QHE_Julia) library. To run the code properly you have to modify all the lines the preamble of the script:

`include("/home/trung/_qhe-julia/<something>.jl")`

to change `/home/trung/_qhe-julia` into the path to wherever the QHE_Julia library is stored on your machine.

To view the available option tags:

`julia sphere_FQHE_widebump_onlynorth_parallel.jl -h`

Example usage:

`julia -t 20 sphere_FQHE_widebump_onlynorth_parallel.jl -e 6 -o 16 -n 15 -i v1.txt -p --method Krylov`

`julia -t 20 sphere_FQHE_widebump_samesign_parallel.jl -e 6 -o 16 -n 15 -i v1.txt -p --method Krylov`

`julia -t 20 sphere_FQHE_widebump_4pins_parallel.jl -e 6 -o 16 -n 15 -i v1.txt -p --method Krylov`

## Output and Plotting
All output files are included within a subdirectory named like `wide_pin_k_<k>_<options>` where `<k>` is an integer specifying the size of the pin (how many orbitals it affect) and `<option>` is a string specifying the specific script (e.g. `north`, `positive_north`, `4pins`, etc.). Note that if the `<option>` part of the directory name doesn't contain `positive`, it means that the results in that directory is for negative pins.

Each of the subdirectories contains further sub-subdirectories, each containing the results of a specific system size (number of electrons and number of orbitals) and value of $$\lambda$$. Each of those sub-subdirectories contains:
* One or several `.txt` file(s) containing the spectrum. The first column is the energy and the second column is the average angular momentum $$\langle L_z\rangle$$ of the eigenstate. Note that the second column is not integer of half-integer if the system breaks rotation symmetry (such as the case of the 4-pins script). 
* All the eigenstates up to the number of converged eigen pairs provided by the solver. The names of these files start with `g` and they have no extension, but they are just plain text files.

Some plotting scripts are included to visualize the results

* `plot_density_sphere.jl` -- Plot the density of the state on the sphere. Use by `julia plot_density_sphere.jl -f <filename>`
* `plot_gap_wide_pins.py` -- Plot the gap as a function of $$\lambda$$ 
* `plot_gap_4pins.py` -- Plot the gap as a function of $$\lambda$$, used specifically for the four-pin setup.


For the Python scripts, run `python3 <filename>.py -h` to check the available options.




