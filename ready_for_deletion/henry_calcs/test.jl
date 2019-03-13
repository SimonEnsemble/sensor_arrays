#using Pkg
#Pkg.add("PorousMaterials")
#Pkg.add("NPZ")
#Pkg.update()
#Pkg.instantiate()
using PorousMaterials
using NPZ

frame = Framework("LICGOW.cif")
strip_numbers_from_atom_labels!(frame)
ljforcefield = LJForceField("UFF.csv", cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")
molecule = Molecule("He")
reps = replication_factors(frame.box, 14.0)
new_frame = replicate(frame, reps)
temperature = 298.0

box_length = 30.0
nb_pts = 16

box = Box(box_length, box_length, box_length, pi/2, pi/2, pi/2)
grid = Grid(box, (nb_pts, nb_pts, nb_pts), zeros(UInt64, nb_pts, nb_pts, nb_pts), :occupancy, [0.0, 0.0, 0.0])
grid_pts = collect(range(0, stop=box_length, length=nb_pts))

for i = 1:nb_pts, j = 1:nb_pts, k = 1:nb_pts
    x_grid_pt = [grid_pts[i], grid_pts[j], grid_pts[k]]
    translate_to!(molecule, x_grid_pt, new_frame.box)
    energy = vdw_energy(new_frame, molecule, ljforcefield)
    grid.data[i, j, k] = energy < 0.0
end

npzwrite("occu/" * split(frame.name, ".")[1] * "_occ.npy", grid.data)
