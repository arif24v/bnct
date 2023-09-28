import openmc
import os
import numpy as np
import neutronics_material_maker as nmm

# Define materials
# Tumor Material
tumor_material = openmc.Material()
tumor_material.add_nuclide('H1', 0.106)
tumor_material.add_nuclide('B10', 19.3e-4)
tumor_material.add_nuclide('B11', 85.3e-4)
tumor_material.add_nuclide('C12', 0.14)
tumor_material.add_nuclide('N14', 0.0184)
tumor_material.add_nuclide('O16', 0.726)
tumor_material.add_nuclide('Na22', 0.0014)
tumor_material.add_nuclide('P31', 0.0039)
tumor_material.add_element('Cl', 0.0014)
tumor_material.add_element('K', 0.0039)
tumor_material.set_density('g/cm3',3)

# Regular Skin Material
skin_material = nmm.Material.from_library(name='Skin (ICRP)', material_id=1)
skin_material = skin_material.openmc_material


# Define materials collection and export to XML
my_materials = openmc.Materials([tumor_material, skin_material])
my_materials.export_to_xml()

# Set cross sections path
os.environ['OPENMC_CROSS_SECTIONS'] = "/Volumes/Untitled/cross_test/endfb-viii.0-hdf5/cross_sections.xml"


# Change Tumor Radius, Change Neutron Source Energy, Change B-10 Concentration
# Define radii in cm
inner_radius = 0.4
middle_radius = 0.5

# Transmission boundary (default) allows neutrons to cross
inner_cylinder = openmc.ZCylinder(r=inner_radius, boundary_type='transmission')
middle_cylinder = openmc.ZCylinder(r=middle_radius, boundary_type='vacuum')
top_plane = openmc.ZPlane(z0=10.0, boundary_type='vacuum') #10 cm high
bottom_plane = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')

# Define cells
inner_slab = openmc.Cell(region=-inner_cylinder, fill=skin_material)
middle_slab = openmc.Cell(region=+inner_cylinder & -middle_cylinder, fill=tumor_material)

# Define geometry
geometry = openmc.Geometry([inner_slab, middle_slab])
geometry.export_to_xml()

# Define source
point = openmc.stats.Point((0, 0, 0))

energy_distribution = openmc.stats.PowerLaw(1, 1000, -1.001)
source = openmc.Source(space=point, energy=energy_distribution)

# Define settings
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 0
settings.particles = 10000
settings.run_mode = 'fixed source'
settings.source = source

# https://www.w3resource.com/numpy/array-creation/logspace.php
# neutrons from 10^-5 to 10^10 eV with 50 equally space bins
energy_bins = np.logspace(-5, 4, 50)

# Define filters
neutron_particle_filter = openmc.ParticleFilter(['neutron'])
inner_slab_filter = openmc.CellFilter(inner_slab)
middle_slab_filter = openmc.CellFilter(middle_slab)

# Define Dose Coefficients
energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(particle="neutron", geometry="AP")
energy_function_filter_n = openmc.EnergyFunctionFilter(energy_bins_n, dose_coeffs_n)
energy_function_filter_n.interpolation = "cubic"  # cubic interpolation is recommended by ICRP

# Define tallies
# Calculate dose using dose coefficients
inner_flux = openmc.Tally(name="inner_flux")
inner_flux.filters = [inner_slab_filter, neutron_particle_filter]
inner_flux.scores=["flux"]

middle_flux = openmc.Tally(name="middle_flux")
middle_flux.filters = [middle_slab_filter, neutron_particle_filter]
middle_flux.scores=["flux"]

tumor_na = openmc.Tally(name="tumor_na")
tumor_na.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n, 
                    openmc.EnergyFilter(energy_bins)]
tumor_na.scores = ["(n,a)"]

tumor_ga = openmc.Tally(name="tumor_ga")
tumor_ga.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n, 
                    openmc.EnergyFilter(energy_bins)]
tumor_ga.scores = ["(n,gamma)"]

check_tumor_ga = openmc.Tally(name="check_tumor_ga")
check_tumor_ga.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n]
check_tumor_ga.scores = ["(n,gamma)"]

tumor_np = openmc.Tally(name="tumor_np")
tumor_np.filters = [middle_slab_filter, neutron_particle_filter, energy_function_filter_n, 
                             openmc.EnergyFilter(energy_bins)]
tumor_np.scores = ["(n,p)"]

skin_na = openmc.Tally(name="skin_na")
skin_na.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n,
                   openmc.EnergyFilter(energy_bins)]
skin_na.scores = ["(n,a)"]

skin_ga = openmc.Tally(name="skin_ga")
skin_ga.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n,
                   openmc.EnergyFilter(energy_bins)]
skin_ga.scores = ["(n,gamma)"]


skin_np = openmc.Tally(name="skin_np")
skin_np.filters = [inner_slab_filter, neutron_particle_filter, energy_function_filter_n,
                   openmc.EnergyFilter(energy_bins)]
skin_np.scores = ["(n,p)"]


my_tallies = openmc.Tallies([inner_flux,middle_flux,tumor_na,tumor_np,tumor_ga,skin_na,skin_np,
                             skin_ga,check_tumor_ga])
settings.tallies = my_tallies

# Create the model and run simulation
model = openmc.model.Model(geometry, my_materials, settings, my_tallies)
openmc_exec = '/Users/arifv/opt/anaconda3/envs/new_env/bin/openmc'
!rm *.h5
model.run(openmc_exec=openmc_exec)


statepoint = openmc.StatePoint('statepoint.10.h5')


inner_flux = statepoint.tallies[inner_flux.id]
middle_flux = statepoint.tallies[middle_flux.id]
print("Inner Flux Before .mean method:",inner_flux)
print("Middle Flux Before .mean method:",middle_flux)
inner_flux = statepoint.tallies[inner_flux.id].mean
middle_flux = statepoint.tallies[middle_flux.id].mean
print("Inner Flux After .mean method:",inner_flux)
print("Middle Flux After .mean method:", middle_flux)

# Define RBE values for various interactions
tumor_doseDict = {1: {'name': 'tumor_na', 'conversion': 20,'RBE':3.8}, #n,a
            2: {'name': 'tumor_ga', 'conversion': 1, 'RBE':1},   #n,gamma
            3: {'name': 'tumor_np', 'conversion': 1, 'RBE': 3.2}    #n,p
           }

doseDict_check = {1: {'name': 'neutron_na_s_dose_on_cell', 'conversion': 20,'RBE':3.8}, #n,a
            2: {'name': 'neutron_ga_s_dose_on_cell', 'conversion': 1, 'RBE':1},   #n,gamma
            3: {'name': 'neutron_p_s_dose_on_cell', 'conversion': 1, 'RBE': 3.2}    #n,p
           }

skin_doseDict = {1: {'name': 'skin_na', 'conversion': 20,'RBE':1.3}, #n,a
                2: {'name': 'skin_ga', 'conversion': 1, 'RBE':1},   #n,gamma
                3: {'name': 'skin_np', 'conversion': 1, 'RBE': 3.2}    #n,p
               }

D_TOT = 0
reg_D_TOT = 0
neutrons_per_second = 1e12 # Subject to change


tumor_na_results = statepoint.get_tally(name='tumor_na').mean.flatten()
tumor_ga_results = statepoint.get_tally(name='tumor_ga').mean.flatten()
tumor_np_results = statepoint.get_tally(name='tumor_np').mean.flatten()

skin_na_results = statepoint.get_tally(name='skin_na').mean.flatten()
skin_ga_results = statepoint.get_tally(name='skin_ga').mean.flatten()
skin_np_results = statepoint.get_tally(name='skin_np').mean.flatten()


tumor_dose = 0
skin_dose = 0

for i in range (1,4):
    for j, tally_result in enumerate(statepoint.get_tally(name=tumor_doseDict[i]['name']).mean.flatten()):
        tally_result = tally_result * middle_flux.flatten() * neutrons_per_second
        middle_slab_volume_cm3 = np.pi*(middle_radius**2 - inner_radius**2)
        tally_result = tally_result/middle_slab_volume_cm3
        tally_dose_pSv = tally_result * 60 * 60
        tally_dose_Sv = tally_dose_pSv * (1e-12)
        tally_dose_Gy = tally_dose_Sv/tumor_doseDict[i]['conversion']
        tumor_dose= tumor_dose + tally_dose_Gy * tumor_doseDict[i]['RBE']
    print(i,"th",tumor_dose)
for i in range (1,4):
    for j, tally_result in enumerate(statepoint.get_tally(name=skin_doseDict[i]['name']).mean.flatten()):
        tally_result = tally_result * inner_flux.flatten() * neutrons_per_second
        inner_slab_volume_cm3 = np.pi*(inner_radius**2)
        tally_result = tally_result/inner_slab_volume_cm3
        tally_dose_pSv = tally_result * 60 * 60
        tally_dose_Sv = tally_dose_pSv * (1e-12)
        tally_dose_Gy = tally_dose_Sv/skin_doseDict[i]['conversion']
        skin_dose = skin_dose + tally_dose_Gy * skin_doseDict[i]['RBE']
    print(i,"th",skin_dose)  
print ("DTOT/hr:", tumor_dose)
print ("reg_DTOT/hr:", skin_dose)
