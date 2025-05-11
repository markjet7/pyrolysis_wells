
import biosteam as bst 

co2 = bst.Chemical("CO2")
# co2.Cn.g.add_method(29.0, Tmax=90)
co = bst.Chemical("CO", phase='g')
h2o = bst.Chemical("H2O")
h2 = bst.Chemical("H2", phase='g')
n2 = bst.Chemical("N2", phase='g')
o2 = bst.Chemical("O2", phase='g')
ch4 = bst.Chemical("CH4", phase='g')
so2 = bst.Chemical("SO2", phase='g')


cellulose = bst.Chemical('Cellulose',
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='s',
                  formula="C6H10O5", # Glucose monomer minus water, molecular weight is computed based on formula
                Hf=-975708.8)
lignin = bst.Chemical("Lignin",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='s',
                  formula="C9H10O3", # 
                  Hf=-975708.8)
hemicellulose = bst.Chemical("Hemicellulose",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='s',
                  formula="C5H8O4", # 
                  Hf=-975708.8)

phenol = bst.Chemical("Phenol",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='l',
                  formula="C6H6O", # 
                  Hf=-975708.8)

guaiacol = bst.Chemical("Guaiacol",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='l',
                  formula="C7H8O2", # 
                  Hf=-975708.8)

furfural = bst.Chemical("Furfural",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='l',
                  formula="C5H4O2", # 
                  Hf=-975708.8)

aceticAcid = bst.Chemical("AceticAcid",
                  Cp=1.364, # Heat capacity [kJ/kg]
                  rho=1540, # Density [kg/m3]
                  default=True, # Default other chemicals properties like viscosity to that of water at 25 C
                  search_db=False, # Not in database, so do not search the database
                  phase='l',
                  formula="CH3COOH", # 
                  Hf=-975708.8)

ash = bst.Chemical("Ash",
                   search_ID="SiO",
                   rho=1540,
                   default=True,
                   Cp=0,
                   Hf=0,
                   phase='s')

mws = [12, 60]
ut = [41.8+19.7,100 - 41.8 - 19.7]
n = [u/m for (u, m) in zip(ut, mws)]
# biochar = bst.Chemical("Biochar",
#                   Cp=1.364, # Heat capacity [kJ/kg]
#                   rho=1540, # Density [kg/m3]
#                   default=True, # Default other chemicals properties like viscosity to that of water at 25 C
#                   search_db=False, # Not in database, so do not search the database
#                   phase='s',
#                   formula="C5.12Si0.642", # 
#                   Hf=-975708.8)

biochar = bst.Chemical("Biochar", search_ID="C")

mws = [12, 1.01, 14, 32, 16]
ut = [61.9,	6.04,	1.95,	0.05,	30.06]
n = [u/m for (u, m) in zip(ut, mws)]
# heavyEnds = bst.Chemical('HeavyEnds',
#                   Cp=1.364, # Heat capacity [kJ/kg]
#                   rho=1540, # Density [kg/m3]
#                   default=True, # Default other chemicals properties like viscosity to that of water at 25 C
#                   search_db=False, # Not in database, so do not search the database
#                   phase='l',
#                   formula="C5.2H5.98N0.14S0.001O1.9", # 
#                   Hf=-975708.8)
heavyEnds = bst.Chemical("HeavyEnds", search_ID="Guiacol")


mws = [12, 1.01, 14, 32, 16]
ut = [21.12,	3.56,	1.29,	0.05,	73.98]
n = [u/m for (u, m) in zip(ut, mws)]
# lightEnds = bst.Chemical('LightEnds',
#                   Cp=1.364, # Heat capacity [kJ/kg]
#                   rho=1540, # Density [kg/m3]
#                   default=True, # Default other chemicals properties like viscosity to that of water at 25 C
#                   search_db=False, # Not in database, so do not search the database
#                   phase='l',
#                   formula="C1.76H3.52N0.092S0.002O4.62", # 
#                   Hf=-975708.8)
lightEnds = bst.Chemical("LightEnds", search_ID="AceticAcid")

# Fix viscosity errors
for chemical in (cellulose, lignin, hemicellulose, phenol, guaiacol, furfural, aceticAcid,
                #   heavyEnds, lightEnds, biochar,
                  ash):
    chemical.mu.add_method(1e-7)

chemicals = bst.Chemicals([
    co2, so2, co, h2o, h2, n2, o2, ch4, cellulose, lignin, hemicellulose, phenol, guaiacol, furfural, aceticAcid,heavyEnds, lightEnds, ash, biochar
])
#%%