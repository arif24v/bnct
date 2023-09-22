import openmc
import os
import numpy as np

# Define materials

# Healthy Skin
reg_skin_material = openmc.Material()
reg_skin_material.add_nuclide('B10',3.2e-4)
reg_skin_material.add_nuclide('C12', 0.116751142398)
reg_skin_material.add_nuclide('C13', 0.001307857602)
reg_skin_material.add_nuclide('Ca40', 2.229643e-05)
reg_skin_material.add_nuclide('Ca42', 1.4881e-07)
reg_skin_material.add_nuclide('Ca43', 3.105e-08)
reg_skin_material.add_nuclide('Ca44', 4.7978e-07)
reg_skin_material.add_nuclide('Ca46', 9.2e-10)
reg_skin_material.add_nuclide('Ca48', 4.301e-08)
reg_skin_material.add_nuclide('Cl35', 0.000354578796)
reg_skin_material.add_nuclide('Cl37', 0.000113421204)
reg_skin_material.add_nuclide('H1', 0.6198694464951601)
reg_skin_material.add_nuclide('H2', 9.655350484e-05)
reg_skin_material.add_nuclide('K39', 0.000125898435)
reg_skin_material.add_nuclide('K40', 1.5795e-08)
reg_skin_material.add_nuclide('K41', 9.08577e-06)
reg_skin_material.add_nuclide('N14', 0.020513582493)
reg_skin_material.add_nuclide('N15', 7.5417507e-05)
reg_skin_material.add_nuclide('O16', 0.23977811121)
reg_skin_material.add_nuclide('O17', 9.109265e-05)
reg_skin_material.add_nuclide('O18', 0.00048079613999999995)
reg_skin_material.add_nuclide('P31', 6.6e-05)
reg_skin_material.add_nuclide('S32', 0.0002927254792)
reg_skin_material.add_nuclide('S33', 2.3059652e-06)
reg_skin_material.add_nuclide('S34', 1.29236492e-05)
reg_skin_material.add_nuclide('S36', 4.49064e-08)
reg_skin_material.add_nuclide('Zn64', 4.917e-07)
reg_skin_material.add_nuclide('Zn66', 2.773e-07)
reg_skin_material.add_nuclide('Zn67', 4.04e-08)
reg_skin_material.add_nuclide('Zn68', 1.845e-07)
reg_skin_material.add_nuclide('Zn70', 6.1e-09)
reg_skin_material.set_density('g/cm3', 1.5)

# Pure Boron-10
boron_material = openmc.Material()
boron_material.add_nuclide('B10', 1)
boron_material.set_density('g/cm3', 7)

# Pure H2O
water_material = openmc.Material()
water_material.add_element('H', 2)
water_material.add_element('O', 1)
water_material.set_density('g/cm3', 1)

# Tumor Material
skin_material = openmc.Material()
skin_material.add_nuclide('H1', 0.106)
skin_material.add_nuclide('B10', 19.3e-4)
skin_material.add_nuclide('B11', 85.3e-4)
skin_material.add_nuclide('C12', 0.14)
skin_material.add_nuclide('N14', 0.0184)
skin_material.add_nuclide('O16', 0.726)
skin_material.add_nuclide('Na22', 0.0014)
skin_material.add_nuclide('P31', 0.0039)
skin_material.add_element('Cl', 0.0014)
skin_material.add_element('K', 0.0039)
skin_material.set_density('g/cm3',3)




# Define materials collection and export to XML
my_materials = openmc.Materials([boron_material,water_material, skin_material, reg_skin_material])
my_materials.export_to_xml()

# Set cross sections path
os.environ['OPENMC_CROSS_SECTIONS'] = "/Users/arifv/Desktop/BNCT/lib80x_hdf5/cross_sections.xml"
openmc.config['cross_sections'] = "/Users/arifv/Desktop/BNCT/lib80x_hdf5/cross_sections.xml"
openmc.data.cross_sections = "/Users/arifv/Desktop/BNCT/lib80x_hdf5/cross_sections.xml"

# Define cylinders


#Change Tumor Radius, Change Neutron Source Energy, Change B-10 Concentration

# Define radii
inner_radius = 2
middle_radius = 5
outer_radius=12

# Transmission boundary (default) allows neutrons to cross
inner_cylinder = openmc.ZCylinder(r=inner_radius, boundary_type='transmission')
middle_cylinder = openmc.ZCylinder(r=middle_radius, boundary_type='transmission')
outer_cylinder = openmc.ZCylinder(r=outer_radius, boundary_type='vacuum')
#top_plane = openmc.ZPlane(z0=10.0, boundary_type='reflective')
#bottom_plane = openmc.ZPlane(z0=-10.0, boundary_type='reflective')

# Define cells
inner_slab = openmc.Cell(region=-inner_cylinder, fill=reg_skin_material)
middle_slab = openmc.Cell(region=+inner_cylinder & -middle_cylinder, fill=skin_material)
outer_slab = openmc.Cell(region=+middle_cylinder & -outer_cylinder, fill=water_material)

# Define geometry
geometry = openmc.Geometry([inner_slab, middle_slab, outer_slab])

# Export geometry to XML
geometry.export_to_xml()

# Define source
point = openmc.stats.Point((0, 0, 0))
#energy_distribution = openmc.stats.Discrete([1e6], [1])
lower = 1
upper = 1e3
energy_distribution = openmc.stats.PowerLaw(1, 1e3, -1)
source = openmc.Source(space=point, energy=energy_distribution)

# Define settings
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 0
settings.particles = 10000
settings.run_mode = 'fixed source'
settings.source = source

# Define filters
neutron_particle_filter = openmc.ParticleFilter(['neutron'])
inner_slab_filter = openmc.CellFilter(inner_slab)
middle_slab_filter = openmc.CellFilter(middle_slab)

#energy_bins = np.logspace(np.log10(lower*0.001), np.log10(upper*0.001), 50)

# as you increase the bins, the sum of all of them approaches the true reaction rate, 
# can express this as an integral to find total reaction rate over different neutron energies
# increase upper to get more accurate result
# https://www.w3resource.com/numpy/array-creation/logspace.php
# neutrons from 10^-5 to 10^10 eV with 50 equally space bins
energy_bins = np.logspace(-5, 4, 50)


# Define tallies

# Reaction Rates (not important)
boron_tally = openmc.Tally(1)
boron_tally.scores = ['(n,a)']
boron_tally.filters = [neutron_particle_filter, inner_slab_filter, openmc.EnergyFilter(energy_bins)]


nitrogen_tally = openmc.Tally(2)
nitrogen_tally.scores = ['(n,a)']
nitrogen_tally.filters = [neutron_particle_filter, middle_slab_filter]

# Calculate dose using dose coefficients
energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(particle="neutron", geometry="AP")
energy_function_filter_n = openmc.EnergyFunctionFilter(energy_bins_n, dose_coeffs_n)
energy_function_filter_n.interpolation = "cubic"  # cubic interpolation is recommended by ICRP


flux_inner_tally = openmc.Tally(name="flux_inner_tally")
flux_inner_tally.filters = [inner_slab_filter, neutron_particle_filter]
flux_inner_tally.scores=["flux"]

flux_middle_tally = openmc.Tally(name="flux_middle_tally")
flux_middle_tally.filters = [middle_slab_filter, neutron_particle_filter]
flux_middle_tally.scores=["flux"]

# (n,a) tally for tumor region
dose_cell_tally_na = openmc.Tally(name="neutron_na_dose_on_cell")
dose_cell_tally_na.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n]
dose_cell_tally_na.scores = ["(n,a)"]

# (n,gamma) tally for tumor region
dose_cell_tally_ga = openmc.Tally(name="neutron_ga_dose_on_cell")
dose_cell_tally_ga.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n]
dose_cell_tally_ga.scores = ["(n,gamma)"]

# (n,p) tally for tumor region
dose_cell_tally_p = openmc.Tally(name="neutron_p_dose_on_cell")
dose_cell_tally_p.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n]
dose_cell_tally_p.scores = ["(n,p)"]


# (n,a) tally for healthy skin region
reg_dose_cell_tally_na = openmc.Tally(name="reg_neutron_na_dose_on_cell")
reg_dose_cell_tally_na.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n]
reg_dose_cell_tally_na.scores = ["(n,a)"]

# (n,gamma) tally for healthy skin region
reg_dose_cell_tally_ga = openmc.Tally(name="reg_neutron_ga_dose_on_cell")
# Include the EnergyFunctionFilter as a filter
reg_dose_cell_tally_ga.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n]
reg_dose_cell_tally_ga.scores = ["(n,gamma)"]

# (n,p) for healthy skin region
reg_dose_cell_tally_p = openmc.Tally(name="reg_neutron_p_dose_on_cell")
# Include the EnergyFunctionFilter as a filter
reg_dose_cell_tally_p.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n]
reg_dose_cell_tally_p.scores = ["(n,p)"]


#cylindrical_mesh = openmc.CylindricalMesh.from_domain(
#    domain=geometry, # the corners of the mesh are being set automatically to surround the geometry
#    dimension=[1, 1, 1] # 100 voxels in each axis direction (r, z, phi)
#)

#mesh_filter = openmc.MeshFilter(cylindrical_mesh)
#cyl_tally = openmc.Tally()
#cyl_tally.scores = ['(n,a)']
#cyl_tally.filters = [mesh_filter, neutron_particle_filter, energy_function_filter_n]

#cyl_tally1 = openmc.Tally()
#cyl_tally1.scores = ['(n,p)']
#cyl_tally1.filters = [mesh_filter, neutron_particle_filter, energy_function_filter_n]

#cyl_tally2 = openmc.Tally()
#cyl_tally2.scores = ['(n,gamma)']
#cyl_tally2.filters = [mesh_filter, neutron_particle_filter, energy_function_filter_n]


my_tallies = openmc.Tallies([boron_tally, nitrogen_tally, dose_cell_tally_na, 
                             dose_cell_tally_ga, dose_cell_tally_p, reg_dose_cell_tally_na, 
                             reg_dose_cell_tally_ga, reg_dose_cell_tally_p, flux_middle_tally, flux_inner_tally])


settings.tallies = my_tallies




# Create the model and run simulation
model = openmc.model.Model(geometry, my_materials, settings, my_tallies)
openmc_exec = '/Users/arifv/opt/anaconda3/envs/new_env/bin/openmc'
!rm *.h5
model.run(openmc_exec=openmc_exec)


statepoint = openmc.StatePoint('statepoint.10.h5')

boron_tally_result = statepoint.tallies[boron_tally.id].mean
nitrogen_tally_result = statepoint.tallies[nitrogen_tally.id].mean
middle_flux = statepoint.tallies[flux_middle_tally.id].mean
inner_flux = statepoint.tallies[flux_inner_tally.id].mean

# Print out the (n,a) reaction rate
print("Inner total reaction rate:", np.sum(boron_tally_result))
print("(n,a) reaction rate:", nitrogen_tally_result)
print("middle flux:", middle_flux)
print("inner flux:", inner_flux)


# Define RBE values for various interactions

doseDict = {1: {'name': 'neutron_na_dose_on_cell', 'conversion': 20,'RBE':3.8}, #n,a
            2: {'name': 'neutron_ga_dose_on_cell', 'conversion': 1, 'RBE':1},   #n,gamma
            3: {'name': 'neutron_p_dose_on_cell', 'conversion': 1, 'RBE': 3.2}    #n,p
           }

reg_doseDict = {1: {'name': 'reg_neutron_na_dose_on_cell', 'conversion': 20,'RBE':1.3}, #n,a
                2: {'name': 'reg_neutron_ga_dose_on_cell', 'conversion': 1, 'RBE':1},   #n,gamma
                3: {'name': 'reg_neutron_p_dose_on_cell', 'conversion': 1, 'RBE': 3.2}    #n,p
               }

D_TOT = 0
reg_D_TOT = 0
neutrons_per_second = 1e12 # Subject to change

for i in range (1,4):

    # Use the StatePoint data to get the tally result, pSv-cm^3/source neutron
    #figure out what mean.flatten[0] does
    print(statepoint.get_tally(name=doseDict[i]['name']).mean)
    tally_result = statepoint.get_tally(name=doseDict[i]['name']).mean.flatten()[0]
    print(tally_result)

    #pSv-cm^3/second
    tally_result = tally_result * middle_flux * neutrons_per_second

    middle_slab_volume_cm3 = np.pi*(middle_radius**2 - inner_radius**2)

    #psv/second
    tally_result = tally_result/middle_slab_volume_cm3

    #60 seconds in a minute, 60 minutes in an hour, x hours
    #pSv/hour
    tally_dose_pSv = tally_result * 60 * 60

    tally_dose_Sv = tally_dose_pSv * 1e-12

    tally_dose_Gy = tally_dose_Sv/doseDict[i]['conversion']
    
    #Gy/hour alpha particles
    print("Total dose in Gy/hr:", tally_dose_Gy,doseDict[i]['name'])
    
    D_TOT= D_TOT + tally_dose_Gy * doseDict[i]['RBE']
for i in range (1,4):

    tally_result = statepoint.get_tally(name=reg_doseDict[i]['name']).mean.flatten()[0]

    #pSv-cm^3/second
    tally_result = tally_result * inner_flux * neutrons_per_second

    inner_slab_volume_cm3 = np.pi*(inner_radius**2)

    #psv/second
    tally_result = tally_result/inner_slab_volume_cm3

    #60 seconds in a minute, 60 minutes in an hour
    #pSv/hour
    tally_dose_pSv = tally_result * 60 * 60

    tally_dose_Sv = tally_dose_pSv * 1e-12

    tally_dose_Gy = tally_dose_Sv/reg_doseDict[i]['conversion']
    
    #Gy/hour alpha particles
    print("Total dose in Gy/hr:", tally_dose_Gy, reg_doseDict[i]['name'])
    
    reg_D_TOT= reg_D_TOT + tally_dose_Gy * reg_doseDict[i]['RBE']
print ("DTOT/hr:", D_TOT)
print ("reg_DTOT/hr:", reg_D_TOT)


boron_tally_results = statepoint.tallies[boron_tally.id].mean
energy_integrated_result = np.sum(boron_tally_results)  # Sum over all energy bins

# Print the energy-integrated result
print("Energy-Integrated Boron Tally Result:", energy_integrated_result)

# Access results within specific energy bins
for i, energy_bin_result in enumerate(boron_tally_results):
    print(f"Energy Bin {i}: {energy_bins[i]} to {energy_bins[i+1]} eV")
    print("Tally Result:", energy_bin_result)
