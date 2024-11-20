# this example code is designed for the spin Monte Carlo modelling of Langevin spin dynamics for a quantum magnet. The example here is a Heisenberg XXZ exchange parameter.

using Sunny, WGLMakie, Statistics, Plots, ColorSchemes, LinearAlgebra
using PyPlot, SciPy
using CSV, DataFrames

units = Units(:meV, :angstrom)
cryst = Crystal("Eu5Sn2As6.cif", symprec = 1e-5) # this crystal may be replaced with another .cif file for simulation

x_dim = 3
y_dim = 3
z_dim = 2 #system dimensions by unit cell


function save_data(vec1, vec2, vec3, vec4, filename, temp)
    output_directory = "output"

    if !isdir(output_directory)
        mkdir(output_directory)
    end
    filename = string("cryst", "_", filename, "_", "temp=$temp", "K")
    output_file_path = joinpath(output_directory, "$filename.csv")

    df = DataFrame(x_magn = vec1, y_magn = vec2, z_magn = vec3, field = vec4)
    if length(vec1)!=length(vec2)!=length(vec3)!=length(vec4)
        print("not same length!")
    end
    
    CSV.write(output_file_path, df)
end


end
spinvec = [1 => Moment(s=7/2, g=-2), 5 => Moment(s=7/2, g=-2), 9 => Moment(s=7/2, g=-2)] # magnitude of Spin vector (working example is Europium, with S=7/2
sys = System(cryst, dims = (x_dim, y_dim, z_dim), spinvec, :dipole; seed=0)

j1 = 1.0
j2 = 1.0
j3 = -1.0 #three exchange parameters for the XXZ Hamiltonian. This may be changed in when bonds are assigned for different spin models.

bond_list = reference_bonds(cryst, 6.0)

for bond in bond_list
    set_exchange!(sys, j1, bond)
end

# this function sets spins to an initial spin configuration. These can be editted to fit Neutron Scattering data, or other experimentally verified magnetic structures.
function set_spins()
    x_list = collect(1:1:x_dim)
    y_list = collect(1:1:y_dim)
    z_list = collect(1:1:z_dim)
    for i in x_list
        for j in y_list
            for k in z_list
                set_dipole!(sys, [0.260565, -0.437558, 0], (i, j, k, 1))
                set_dipole!(sys, [-0.594414, 0.210171, 0], (i, j, k, 2)) #2
                set_dipole!(sys, [-0.594414, 0.210171, 0], (i, j, k, 3)) #3
                set_dipole!(sys, [0.260565, -0.437558, 0], (i, j, k, 4)) #4
                set_dipole!(sys, [0.578943, 0.100423, 0], (i, j, k, 5)) #5
                set_dipole!(sys, [0.578943, 0.100423, 0], (i, j, k, 6)) #6
                set_dipole!(sys, [0.143311, 0.545155, 0], (i, j, k, 7)) #7
                set_dipole!(sys, [0.578943, 0.1000423, 0], (i, j, k, 8)) #8
                set_dipole!(sys, [-0.219852, -0.573847, 0], (i, j, k, 9)) #9
                set_dipole!(sys, [-0.578129, 0.1865, 0], (i, j, k, 10)) #10
            end
        end
    end
end


#this section of the code implements the Monte Carlo algorithm to sample thermodynamically favorable magnetic ground states at different applied magnetic field strengths. 
# the function mvsb sets the external field to a point defined in B_list, and then calculates the resulting energy change of spin flips to find a favorable state from the phase space. 
# as it stands this output is saved as a .csv file which can separately be plotted.

B_up = collect(0:.1:6)
B_back = collect(6:-.1:-6)
B_loop = collect(-6:.1:6) 
B_list = vcat(B_up, B_back, B_loop)

function mvsb(temp, list, dir, filename, damp)
    nsweeps = 10000
    sampler = Langevin(0.25, damping = damp, kT =temp)
    x_mag = Float64[]
    y_mag = Float64[]
    z_mag = Float64[]
    for i in 1:length(list)
        print(list[i])
        if dir==1
            set_field!(sys, [list[i], list[i], 0]*units.T)
        elseif dir==2
            set_field!(sys, [list[i], -list[i], 0]*units.T)
        elseif dir==3
            set_field!(sys, [0, 0, list[i]]*units.T)
        else
            print("3rd variable must be a 1, 2, or 3!")
        end
        for i in 1:nsweeps #montecarlo for each step
            step!(sys, sampler)
        end 

        total_mag = sum(sys.dipoles)/(length(sys.dipoles))
        append!(x_mag, total_mag[1])
        append!(y_mag, total_mag[2])
        append!(z_mag, total_mag[3])
            
            
    end
    temper = temp


    save_data(x_mag, y_mag, z_mag, list, filename, temper)
    if dir==1
        xlist = Float64[]
        for i in 1:length(x_mag)
            append!(xlist, dot([x_mag[i], y_mag[i], z_mag[i]], [1, 1, 0]))
        end    
        Plots.plot!(list, xlist, title = "Moment vs. H || x", label = "damp = $damp", palatte = "jet")

    elseif dir==2
        ylist = Float64[]
        for i in 1:length(y_mag)
            append!(ylist, dot([x_mag[i], y_mag[i], z_mag[i]], [1, -1, 0]))
        end
        Plots.plot!(list, ylist, title = "Moment vs. H || y", label = "damp = $damp", palatte = "chrome")
    elseif dir==3
#         for i in 1:length(z_mag)
#             for j in 1:length(z_mag)
#                if (i-j==1 || j-i==1) && z_mag[i]>0 && z_mag[j]<0 
#                    print(" Coercive field $i: $(list[i]) ")
#                end
#             end
#         end
        Plots.plot!(list, z_mag, title = "Moment vs. H || z", label = "T = $temper K", color = "blue")
    else
        print("Error - dir must be 1, 2, or 3")
    end
#     Plots.xlabel!("Field (T)")
#     Plots.ylabel!("Moment (Spin-Ang. Momentum)")


end
