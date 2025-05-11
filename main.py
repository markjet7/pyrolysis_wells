# %%
import biosteam as bst
import numpy as np

np.int = int
import warnings

import pickle 

# warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.append("/Users/markmw/Library/CloudStorage/OneDrive-IowaStateUniversity/Dropbox/Research/Charm/process/")

sys.path.append("/Users/markmw/Library/CloudStorage/OneDrive-IowaStateUniversity/General - SESA Home/Scripts/sesalca")
# change working directory to the directory of this script
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from _compounds import *
from _lca import lca
from _tea_charm import *
# from steamdrawio import draw

from kineticsML import model

import sesalca 

# from sesaCanvas import sesaCanvas 
# sesaCanvas.run()
# %%


bst.settings.set_thermo(chemicals)
# %%
from _utilities import *
from _units import *

# %%
bst.preferences.update(
    flow="tonnes/day", composition=False, N=100, background_color="white"
)

# %%
def create_system(composition={"c":0.365, "h":0.351, "l":0.188, "ash":0.19}):
    bst.main_flowsheet.clear()
    moisture_content = 0.25
    capacity = 10  # tonnes/day
    biomass = bst.Stream(
        "Biomass",
        Cellulose=composition["c"],
        Hemicellulose=composition["h"],
        Lignin=composition["l"],
        Ash=composition["ash"],
        H2O=0,
        units="kg/hr",
    )
    biomass.set_total_flow(capacity, units="tonnes/day")
    biomass.set_flow(moisture_content * biomass.get_total_flow("kg/hr"), "kg/hr", "H2O")
    dried_biomass = bst.Stream("Dried_biomass")
    grinder_power = 75  # kWh/tonne
    screen_unders = bst.Stream("Screen_unders")
    recycle = bst.Stream("Recycle", N2=50, units="tonnes/day")
    with bst.System("Handling_sys") as handling_sys:
        # dryer_air = bst.Stream("Dryer_air", O2=0.21, N2=0.79)

        dryer = bst.units.DrumDryer(
            "Dryer",
            T=343.15,
            P=101325,
            split={"Ash": 0.05},
            ins=(biomass, "Dryer_air", "Dryer_Natural_gas"),
            outs=(dried_biomass, "Hot_Air", "Dryer_Emissions"),
            moisture_content=0.12,
        )
        dryer.ins[2].price = 3*22000/1e6 * 2.205 # $4/mmbtu * 22000 btu/lb * 1/2.205 lb/kg

        screener_overs = bst.Stream("Screener_overs")

        feed_conveyor = bst.units.ConveyingBelt("Feed_conveyor", ins=dried_biomass)
        feed_bin = bst.units.StorageTank("Feed_bin", ins=feed_conveyor-0)

        grinder_mixer = bst.units.Mixer("Grinder_mixer", ins=(feed_bin-0, screener_overs))
        grinder = bst.units.HammerMill("Grinder", ins=grinder_mixer-0)
        grinder.add_power_utility(biomass.get_total_flow("tonnes/hr") * 75)
        screener = bst.units.Splitter(
            "Screener", split=0.9, outs=(screen_unders, screener_overs),
            ins=grinder-0
        )
    # handling_sys.diagram("thorough")

    pyrolysis_vapors = bst.Stream("Pyrolysis_vapors")
    gas_compressor = bst.units.IsentropicCompressor(
        "Recycle_cp", P=10 * 101325, eta=0.95, ins=(recycle)
    )
    with bst.System("Pyrolysis_sys") as pyrolysis_sys:
        pyrolysis_mixer = bst.units.Mixer("Pyrolysis_mixer", ins=(gas_compressor-0, screener-0))
        # pyrolysis = bst.units.MixTank('Pyrolysis')
        pyrolysis = Reactor("Pyrolysis", outs=(pyrolysis_vapors), ins=(pyrolysis_mixer-0))

        # (screener-0, gas_compressor - 0) - pyrolysis_mixer - 0 - pyrolysis
    # pyrolysis_sys.diagram("thorough")

    ncg = bst.Stream("NCG")
    clean_vapors = bst.Stream("Clean_vapors")
    quench_water = bst.Stream("Quench_water", T=300, P=101325, H2O=1)
    with bst.System("Recovery_sys") as recovery_sys:
        pyrolysis_cyclone = bst.units.Splitter(
            "Pyrolysis_cyclone",
            split={"Biochar": 0.98, "Ash": 0.98},
            ins=pyrolysis_vapors,
            outs=("", clean_vapors),
        )
        pyrolysis_quench = bst.units.Mixer("Quench", ins=(pyrolysis_cyclone-1, quench_water))
        pyrolysis_flash = bst.units.Flash(
            "Light_ends_col", P=101325, T=350, outs=(ncg, ""), ins=pyrolysis_quench-0
        )
        pyrolysis_flash.flash_inlet = False
        light_oil_storage = bst.units.StorageTank("Oil_storage", ins=(pyrolysis_flash-1))
        light_oil_storage_pump = bst.units.Pump("Oil_storage_pump", ins=light_oil_storage-0, outs=["Biooil"])

        biochar_storage = bst.units.StorageTank("Biochar_storage", ins=(pyrolysis_cyclone-0), outs=("Biochar"))

    # recovery_sys.diagram("thorough")

    with bst.System("Utilities_sys") as utilities_sys:
        boiler_mixer = bst.units.Mixer("Boiler_mixer", ins=(ncg))
        # boiler = bst.units.MixTank('Boiler', outs=('',))
        boiler = Boiler("HRSG", ins=(boiler_mixer-0,'', ''))
        boiler_splitter = bst.units.Splitter(
            "spBoiler", split=0.1, outs=('', "Flue_gas"), ins=(boiler-0)
        )
        boiler_splitter.outs[0] = gas_compressor.ins[0]

    # system = bst.System(
    #     "sys", path=(handling_sys, pyrolysis_sys, recovery_sys, utilities_sys)
    # )
    system = bst.main_flowsheet.create_system("Biomass_Fast_Pyrolysis")
    boiler.other_units = [u for u in system.units if u.ID != boiler.ID]

    # products: phenol, guaiacol, furfural, aceticAcid, ash, biochar
    biochar_yield = 0.288
    light_end_yield = 0.397
    heavy_end_yield = 0.192
    ncg_yield = 0.192
    pyrolysis_yields = {
        "Biochar": biochar_yield,
        "LightEnds": light_end_yield,
        "HeavyEnds": heavy_end_yield,
        "CO2": 0.118,
        "CO": 0.058,
        "CH4": 0.01,
    }

    def pyrolysis_specification():
        pyrolysis._run()
        biomass_in = sum(pyrolysis.ins[0].imass["Cellulose", "Hemicellulose", "Lignin"])

        predicted_yields = model.predict([[composition["c"], composition["h"], composition["l"]*3]])[0]

        pyrolysis_yields = {
            "Biochar": predicted_yields[1],
            "LightEnds": predicted_yields[0]*0.397/(0.397+0.288),
            "HeavyEnds": predicted_yields[0]*0.288/(0.397+0.288),
            "CO2": predicted_yields[2]* 0.118/(0.118+0.058+0.01),
            "CO": predicted_yields[2]* 0.058/(0.118+0.058+0.01),
            "CH4": predicted_yields[2]* 0.01/(0.118+0.058+0.01),
        }


        normalized = {k:v/sum(list(pyrolysis_yields.values())) for k,v in pyrolysis_yields.items()}
        print(normalized)

        for compound, y in normalized.items():
            pyrolysis.outs[0].imass[compound] = biomass_in * y + pyrolysis.ins[0].imass[compound]
        
        pyrolysis.outs[0].imass["H2O"] = pyrolysis.ins[0].imass["H2O"]
        pyrolysis.outs[0].imass["N2"] = pyrolysis.ins[0].imass["N2"]
        pyrolysis.outs[0].imass["CO2"] = pyrolysis.ins[0].imass["CO2"]
        pyrolysis.outs[0].imass["Ash"] = pyrolysis.ins[0].imass["Ash"]

        pyrolysis.outs[0].imass["Cellulose"] = pyrolysis.ins[0].imass["Cellulose"]*(1-min(1, (sum(list(pyrolysis_yields.values())))))
        pyrolysis.outs[0].imass["Hemicellulose"] = pyrolysis.ins[0].imass["Hemicellulose"]*(1-min(1, (sum(list(pyrolysis_yields.values())))))
        pyrolysis.outs[0].imass["Lignin"] = pyrolysis.ins[0].imass["Lignin"]*(1-min(1, (sum(list(pyrolysis_yields.values())))))

        pyrolysis.outs[0].T = 773.15


    pyrolysis.add_specification(pyrolysis_specification)
    gas_compressor.add_power_utility(0.1)
    system.operating_hours = 0.9 * 365 * 24

    system.flowsheet.stream["Dryer_air"].ID = "Airnoghg"
    system.flowsheet.stream["Dryer_Natural_gas"].ID = "Natural_gas"
    system.flowsheet.stream["Boiler_air"].ID = "Air2noghg"
    system.flowsheet.stream["Quench_water"].ID = "WaterProductionnoghg"

    return system

system = create_system()
#%%
# Clean Pine (CP)	Tulip Poplar (TP)	Hybrid Poplar (HP)	Corn Stover (CS)	Switchgrass (SG)	Switchgrass, Fast Pyrolysis at 450 °C (SG450)	Oriented Strand Board (OSB)

# {"Pine" -> pine, "Tulip Poplar" -> {46.32222222222222, 24.95555555555556, 24.438888888888886}, "Hybrid Poplar" -> {46.32222222222222, 24.95555555555556, 24.438888888888886}, "Corn Stover" -> {38.75, 25.799999999999997, 14.350000000000001}, "Switchgrass" -> {36.980000000000004, 27.380000000000003, 12.960000000000003}, "Oriented Strand Board" -> {}}

feedstocks = {
    "Clean Pine": {"c":0.365, "h":0.351, "l":0.188, "ash":0.19},
    "Tulip Poplar": {"c":0.463, "h":0.249, "l":0.244, "ash":0.044},
    "Hybrid Poplar": {"c":0.463, "h":0.249, "l":0.244, "ash":0.044},
    "Corn Stover": {"c":0.3875, "h":0.258, "l":0.1435, "ash":0.211},
    "Switchgrass": {"c":0.3698, "h":0.274, "l":0.1296, "ash":0.226},
    "Oriented Strand Board": {"c":0.5, "h":0.05, "l":0.4, "ash":0.05},
}
# Woody biomass (CP, TP, and HP)	$109.67 ($99.49)
# CS	$86.31 ($78.30)
# SG and SG450	$87.74 ($79.60)
# OSB	$64.64 ($58.64)a
prices = {
    "Clean Pine": 109.67,
    "Tulip Poplar": 109.67,
    "Hybrid Poplar": 109.67,
    "Corn Stover": 86.31,
    "Switchgrass": 87.74,
    "Oriented Strand Board": 64.64,
}
#%%
system.simulate()
#%%
def update_impact_factors():
    impactFactors = lca(system, feedstock="Clean Pine", product="Biooil")
    
    
    df = pd.DataFrame(impactFactors["flows"])
    df.to_csv("flowImpactFactors.csv")
    
    ffs = {}
    rps  = {}
    for f in feedstocks:
        rp = sesalca.recommend_process(f"{f} production")
        rps[f] = rp
        ffs[f] = sesalca.get_total_impacts(rp, "TRACI")
    ffs 
    df = pd.DataFrame(ffs)
    df.to_csv("feedstockImpactFactors.csv")

#%%
rps = {
    "Clean Pine": "bark chips, wet, measured as dry mass to generic market for residual hardwood, wet | residual softwood, wet | APOS, U",
    "Tulip Poplar": "bark chips, wet, measured as dry mass to generic market for residual hardwood, wet | residual hardwood, wet | APOS, U",
    "Hybrid Poplar": "bark chips, wet, measured as dry mass to generic market for residual hardwood, wet | residual hardwood, wet | APOS, U",
    "Corn Stover": "maize silage production | maize silage | APOS, U",
    "Switchgrass": "alfalfa/grass silage production | alfalfa-grass silage | APOS, U",
    "Oriented Strand Board": "oriented strand board production | residual wood, dry | APOS, U",
}

# ffs = {}
# for f in feedstocks:
#     ffs[f] = sesalca.get_total_impacts(rps[f], "TRACI")

# df = pd.DataFrame(ffs)
# df.to_csv("feedstockImpactFactors.csv")

#%%
flowImpactFactors = pd.read_csv("flowImpactFactors.csv", index_col=0).to_dict()
rawFeedstockImpactFactors = pd.read_csv("feedstockImpactFactors.csv", index_col=0).to_dict()
#%%
woodDensity = 500 # kg/m3
feedstockImpactFactors = {}
feedstockImpactFactors["Clean Pine"] = {k:v/woodDensity for (k,v) in rawFeedstockImpactFactors["Clean Pine"].items()}
feedstockImpactFactors["Tulip Poplar"] = {k:v/woodDensity for (k,v) in rawFeedstockImpactFactors["Tulip Poplar"].items()}
feedstockImpactFactors["Hybrid Poplar"] = {k:v/woodDensity for (k,v) in rawFeedstockImpactFactors["Hybrid Poplar"].items()}
feedstockImpactFactors["Corn Stover"] = {k:v for (k,v) in rawFeedstockImpactFactors["Corn Stover"].items()}
feedstockImpactFactors["Switchgrass"] = {k:v for (k,v) in rawFeedstockImpactFactors["Switchgrass"].items()}
feedstockImpactFactors["Oriented Strand Board"] = {k:v/woodDensity for (k,v) in rawFeedstockImpactFactors["Oriented Strand Board"].items()}
#%%
df = pd.DataFrame(flowImpactFactors)
df.to_clipboard()
#%%
# list(flowImpactFactors.keys())
# ['Unnamed: 0',
#  '2:natural gas, high pressure',
#  '3:electricity, high voltage',
#  '4:bundle, energy wood, measured as dry mass',
#  '5:heat, air-water heat pump 10kW',
#  '6:electricity, high voltage',
#  '7:transport, freight, lorry 28 metric ton, vegetable oil methyl ester 100%',
#  '8:heat, district or industrial, natural gas',
#  '9:transport, tractor and trailer, agricultural',
#  '10:natural gas, high pressure',
#  '12:dryer',
#  '13:charcoal',
#  '14:deep well closure']

# for p in system.products:
#     print(p.ID)

# Hot_Air
# Dryer_Emissions

# Flue_gas
# Biooil
# Biochar

# for f in system.feeds:
#     print(f.ID)

# Biomass
# Airnoghg
# Natural_gas
# WaterProductionnoghg
# Natural_gas
# Boiler_air
#%%
impactNames = ['human health - carcinogenics (kg benzene-Eq)',
 'environmental impact - ozone depletion (kg CFC-11-Eq)',
 'human health - respiratory effects, average (kg PM2.5-Eq)',
 'environmental impact - ecotoxicity (kg 2,4-D-Eq)',
 'environmental impact - eutrophication (kg N)',
 'environmental impact - global warming (kg CO2-Eq)',
 'environmental impact - photochemical oxidation (kg NOx-Eq)',
 'environmental impact - acidification (moles of H+-Eq)',
 'human health - non-carcinogenics (kg toluene-Eq)']
#%%
def get_results(feedstock, sensitivity={}):
    system = create_system(feedstocks[feedstock])
    system.simulate()
    impacts = {}
    distance = 32.19 # km
    impacts["Feedstock Transportation"] = {
        k:float(flowImpactFactors['9:transport, tractor and trailer, agricultural'][k])*distance*system.flowsheet.stream["Biomass"].F_mass/1000*sensitivity.get("Feedstock Transportation", 1) for k in impactNames
    }

    impacts["Oil Transportation"] = {
        k:float(flowImpactFactors['7:transport, freight, lorry 28 metric ton, vegetable oil methyl ester 100%'][k])*distance*system.flowsheet.stream["Biooil"].F_mass*32.19/1000*sensitivity.get("Oil Transportation", 1) for k in impactNames
    }

    impacts["Deep Well Closure"] = {
        k:float(flowImpactFactors['14:deep well closure'][k])*system.flowsheet.stream["Biooil"].F_mass/20*sensitivity.get("Deep Well Closure", 1) for k in impactNames
    }

    impacts["Electricity"] = {
        k:float(flowImpactFactors['3:electricity, high voltage'][k])*system.power_utility.rate/1000*3600*sensitivity.get("Electricity", 1) for k in impactNames
    }

    impacts["Utilities"] = {
        k:float(flowImpactFactors['8:heat, district or industrial, natural gas'][k])*system.get_heating_duty()/1000/system.operating_hours*sensitivity.get("Utilities", 1) for k in impactNames
    }

    impacts["Feedstock Production"] = {
        k:float(feedstockImpactFactors[feedstock][k])*system.flowsheet.stream["Biomass"].F_mass*sensitivity.get("Feedstock Production", 1) for k in impactNames
    }

    # df = pd.DataFrame(impacts)
    # # plot 'environmental impact - global warming (kg CO2-Eq)'
    # df.T["environmental impact - global warming (kg CO2-Eq)"].plot(kind="bar")

    oilCarbon = system.flowsheet.stream["Biooil"].get_atomic_flow("C")*12*44/12
    charCarbon = system.flowsheet.stream["Biochar"].get_atomic_flow("C")*12*44/12*0.7*sensitivity.get("Biochar Carbon", 1)

    kgCO2 = {
        "Feedstock": impacts["Feedstock Production"]["environmental impact - global warming (kg CO2-Eq)"],
        "Feedstock Transportation": impacts["Feedstock Transportation"]["environmental impact - global warming (kg CO2-Eq)"],
        "Oil Transportation": impacts["Oil Transportation"]["environmental impact - global warming (kg CO2-Eq)"],
        "Electricity": impacts["Electricity"]["environmental impact - global warming (kg CO2-Eq)"],
        "Utilities": impacts["Utilities"]["environmental impact - global warming (kg CO2-Eq)"],
        "Deep Well Closure": impacts["Deep Well Closure"]["environmental impact - global warming (kg CO2-Eq)"],
        "Biochar": -charCarbon,
        "Biooil": -oilCarbon,
    }

    tea = TEA(
        system=system
    )
    tea.system.flowsheet.stream["Biomass"].price = prices[feedstock] * 0.03/86.31*sensitivity.get("Feedstock Price", 1) # Adjusted to corn stover price
    msp = tea.solve_price(system.flowsheet.stream["Biooil"])
    table = {k:v/(system.flowsheet.stream["Biooil"].F_mass*tea.operating_days*24) for (k,v) in tea.mfsp_table(system.flowsheet.stream["Biooil"]).items()}

    carbon_abatement_cost = {k:(v * system.flowsheet.stream["Biooil"].F_mass)/abs(sum(list(kgCO2.values()))/1000) for k,v in table.items()}

    carbon_abatement_cost_no_char = {k:(v * system.flowsheet.stream["Biooil"].F_mass)/abs((sum(list(kgCO2.values()))-kgCO2["Biochar"])/1000) for k,v in table.items()}

    gwp = {k:v/system.flowsheet.stream["Biooil"].F_mass for k,v in kgCO2.items()}

    return {
        "msp": msp,
        "co2": co2,
        "gwp": gwp,
        "costs": table,
        "impacts": impacts,
        "output": system.flowsheet.stream["Biooil"].F_mass,
        "carbon_abatement_cost": carbon_abatement_cost,
        "carbon_abatement_cost_no_char": carbon_abatement_cost_no_char
    }

results = {f:get_results(f) for f in feedstocks.keys()}
pickle.dump(results, open("results.pkl", "wb"))

#%%
df = pd.DataFrame({r:results[r]["gwp"] for r in results })
df.T.plot(kind="bar", stacked=True)
df.T.to_csv("plots/gwp.csv")
#%%
df = pd.DataFrame({r:results[r]["costs"] for r in results })
df.T.plot(kind="bar", stacked=True)
df.T.to_csv("plots/costs.csv")
#%%
df = pd.DataFrame({r:results[r]["carbon_abatement_cost"] for r in results })
df.T.plot(kind="bar", stacked=True)
df.T.to_csv("plots/carbon_abatement_cost.csv")
#%%
df = pd.DataFrame({r:results[r]["carbon_abatement_cost_no_char"] for r in results })
df.T.plot(kind="bar", stacked=True)
df.T.to_csv("plots/carbon_abatement_cost_no_char.csv")

#%%
impacts = []
# flatten the impact structure
for f in results.keys():
    for s in results[f]["impacts"].keys():
        for i in results[f]["impacts"][s].keys():
            impacts.append([f, s, i, results[f]["impacts"][s][i]/results[f]["output"]])
df = pd.DataFrame(impacts)
df.to_csv("plots/impacts.csv")
df
#%%
results["Clean Pine"]["output"]
#%%
sensitivity = {
        "Feedstock Transportation": 1,
        "Oil Transportation": 1,
        "Deep Well Closure": 1,
        "Electricity": 1,
        "Utilities": 1,
        "Feedstock Production": 1,
        "Biochar Carbon": 1,
        "Feedstock Price": 1
    }
#%%
sensitivities = {}
for f in feedstocks.keys():
    sensitivities[f] = {}
    for s in sensitivity.keys():
        sensitivities[f][s] = {}
        for k,v in {"low":0.8, "high":1.2, "base":1}.items():
            sensitivity[s] = v
            results = get_results(f, sensitivity)
            sensitivities[f][s][k] = sum(list(results["carbon_abatement_cost"].values()))
        

#%%
# flatten the sensitivity structure
flatten = {}
for f in feedstocks.keys():
    for s in sensitivity.keys():
            flatten[f + ": " + s] = [f, s, sensitivities[f][s]["low"], sensitivities[f][s]["base"], sensitivities[f][s]["high"]]
#%%
df = pd.DataFrame(flatten, index=["Feedstock", "Sensitivity", "Low", "Base", "High"]).T
df 
df.to_csv("plots/sensitivity.csv")
#%%

# %%
# NREL Comparison (Table ES-4. Economic Summary (Modeled) for a 2020 Sensitivity Case with Hydroprocessing at the Biorefinery)
NREL = {
    "Feedstock Handling": 1 / capacity * 480000 * (capacity / 2000) ** 0.7,
    "Fast Pyrolysis*": 1 / capacity * 35000000 * (capacity / 2000) ** 0.7,
    "Quench and Recovery": 1 / capacity * 51590000 * (capacity / 2000) ** 0.7,
    "Utilities": 1 / capacity * (49260000 + 8400000) * (capacity / 2000) ** 0.7,
}

this = {
    "Feedstock Handling": 1
    / capacity
    * sum(
        [
            i.installed_cost
            for i in [feed_conveyor, feed_bin, grinder_mixer, grinder, screener]
        ]
    ),
    "Fast Pyrolysis*": 1 / capacity * sum([i.installed_cost for i in [pyrolysis]]),
    "Quench and Recovery": 1
    / capacity
    * sum(
        [
            i.installed_cost
            for i in [
                pyrolysis_cyclone,
                pyrolysis_quench,
                pyrolysis_light_ends,
                oil_storage_pump,
                oil_storage,
            ]
        ]
    ),
    "Utilities": 1 / capacity * sum([i.installed_cost for i in [boiler_mixer, boiler]]),
}

tea.cost_growth = 0.29
pioneer_capital = tea.FCI()

tea.learning_rate = 0.20
tea.units_manufactured = 1000
pioneer_learning_capital = tea.FCI()

tea.cost_growth = 1
learning_capital = tea.FCI()

tea.cost_growth = 1
tea.learning_rate = 0
tea.units_manufactured = 1
# %%
# change matplotlib theme to an attractive one
# plt.style.use("ggplot")
# # show the NREL and this cost comparison with a bar chart
# # show the NREL bars next to the this bars
# plt.figure(figsize=(14, 10))
# # increase font size
# plt.rcParams.update({"font.size": 16})
# plt.bar(NREL.keys(), NREL.values(), width=0.4, label="NREL-Scaled (2015)")
# # shift the this bars to the right by 0.5
# plt.bar(
#     [i + 0.4 for i in range(len(this))], this.values(), width=0.4, label="This study"
# )
# # add labels to the bars
# for i, v in enumerate(NREL.values()):
#     plt.text(i - 0.2, v, f"${round(v, -2)}")

# for i, v in enumerate(this.values()):
#     plt.text(i + 0.2, v, f"${round(v, -2)}")
# plt.legend()
# plt.ylabel(
#     "Pyrolysis biorefinery installed costs by section\n (USD/tonne/day capacity)"
# )
# plt.grid(axis="y")

# plt.savefig("capital_cost.png", bbox_inches="tight", dpi=300)

# %%

# %%
capital_cost_comparison = {
    "Kim et al. 2014\n(3 tpd)": [
        271 / 236.7 * 1 / capacity * (350000 * (capacity / 3) ** 0.7),
        271 / 236.7 * 1 / capacity * (500000 * (capacity / 3) ** 0.7),
        271 / 236.7 * 1 / capacity * (2050000 * (capacity / 3) ** 0.7),
        1 / capacity * (350000 * (capacity / 3) ** 0.7),
    ],
    "This Study\nPioneer\n20% Learning Rate\n1000 units\n(10 tpd)": [
        pioneer_learning_capital / capacity
    ],
    "This Study\nPioneer\n(10 tpd)": [pioneer_capital / capacity],
    "This Study\n20% Learning Rate\n1000 units\n(10 tpd)": [learning_capital / capacity],
    "NREL 2021\n(2000 tpd)": list(
        NREL.values()
    ),  # Assuming NREL is a dictionary where values are the data points
    "Ganguly et al.\n(250 tpd)": [1 / capacity * 19000000 * (capacity / 250) ** 0.7],
    "Badger 2011": [271 / 224.9 * 1 / capacity * 6030000 * (capacity / 100) ** 0.7],
    "This study\n(10 tpd)": list(
        this.values()
    ),  # Assuming 'this' is a dictionary where values are the data points
    "Sandia 2011\n(50 tpd)": [
        271 / 224.9 * 1 / capacity * 851000 * (capacity / 50) ** 0.7
    ],
}

capital_cost_totals = {}
for key, value in capital_cost_comparison.items():
    capital_cost_totals[key] = sum(value)

capital_cost_totals = pd.DataFrame(capital_cost_totals, index=[0])

# plot the capital costs totals sorted from highest to lowest using a horizontal bar chart label the bars with the total cost and names on the left
capital_cost_totals.sort_values(by=0, axis=1, inplace=True)
capital_cost_totals.transpose().plot.barh(
    figsize=(16, 12),
    legend=False,
    xlabel="Capital Cost (2021 USD/tonne/day capacity)",
    title=f"Pyrolysis Biorefinery\n(Feedstock Capacity: {capacity} tonne/day)",
    color="C1",
)

capital_cost_totals.to_csv("plots\capital_cost_comparison.csv")

plt.savefig("capital_cost_comparison.png", bbox_inches="tight", dpi=300)
# %%
capital_cost_totals.transpose().sort_values(by=0)

# %%
msp = tea.solve_price(oil_storage.outs[0])
print("Bio-oil minimum selling price is ${:.3f}/kg".format(msp))
tea.system.flowsheet.stream["Bio_oil"].price = msp
# %%
opex = pd.DataFrame(tea.mfsp_table(oil_storage.outs[0]), index=[0])

# change colorscheme to one with many colors 
plt.style.use("seaborn-v0_8-colorblind")
opex.to_csv("plots/opex.csv")
opex.plot(kind="bar", stacked=True, figsize=(8, 18))

# reverse the order of the legend laels
handles, labels = plt.gca().get_legend_handles_labels()
plt.gca().legend(
    handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1, 0.5)
)
plt.ylabel("Annual operating cost (USD)")
plt.xlabel("")
plt.xticks([])
# add labels to the bars
acc = 0
accn = 0
for i, v in enumerate(opex.values[0]):
    if v > 300000:
        plt.text(-0.05, acc + v / 2, "$" + str(v / 1e6)[0:3] + "M")
    if v > 0:
        acc += v
    if v < 0:
        plt.text(-0.05, accn + v / 2, "$" + str(v / 1e6)[0:4] + "M")
        accn += v

plt.savefig("opex.png", bbox_inches="tight", dpi=300)


# %%
def plot_msp(filename="msp.png"):
    opex2 = pd.DataFrame(tea.mfsp_table(oil_storage.outs[0]), index=[0])
    opex = opex2 / oil_storage.outs[0].get_total_flow("tonnes/year")
    opex.to_csv("./plots/" + filename.replace(".png", ".csv"))
    opex.plot(kind="bar", stacked=True, figsize=(8, 18))

    # reverse the order of the legend laels
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.gca().legend(
        handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1, 0.5)
    )
    plt.ylabel("Bio-oil Minimum Selling Price (USD/MT)")
    plt.xlabel("")
    plt.xticks([])
    # add labels to the bars
    acc = 0
    accn = 0
    for i, v in enumerate(opex.values[0]):
        if v > 10:
            plt.text(-0.075, acc + v / 2, "$" + str(v)[0:4] + "/MT")
        if v > 0:
            acc += v
        if v < 0:
            plt.text(-0.075, accn + v / 2, "$" + str(v)[0:5] + "/MT")
            accn += v

    plt.text(
        -0.475,
        msp * 1000 + 5,
        "$" + str(msp * 1000)[0:5] + "/MT",
        color="Black",
        backgroundcolor="White",
    )

    # add horizontal line
    plt.axhline(y=msp * 1000, color="black", linestyle="--")

    plt.savefig(filename, bbox_inches="tight", dpi=300)


plot_msp()


# %%
# Sensitivity Analysis
sensitivity = {}

tea.system.simulate()
tea.system.simulate()
baseline = tea.solve_price(oil_storage.outs[0])

base = tea.system.flowsheet.unit["Pyrolysis"].installed_costs["Reactors"]
tea.system.flowsheet.unit["Pyrolysis"].installed_costs["Reactors"] = base * 0.8
low = tea.solve_price(oil_storage.outs[0])
tea.system.flowsheet.unit["Pyrolysis"].installed_costs["Reactors"] = base * 1.2
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"Reactor cost (${base/1e6:2.2f} MM)"] = [high, low]
tea.system.flowsheet.unit["Pyrolysis"].installed_costs["Reactors"] = base
print(tea.solve_price(oil_storage.outs[0]))

# base = tea.system.flowsheet.stream["Biochar"].price
# tea.system.flowsheet.stream["Biochar"].price = base * 0.8
# low = tea.solve_price(oil_storage.outs[0])
# tea.system.flowsheet.stream["Biochar"].price = base * 1.2
# high = tea.solve_price(oil_storage.outs[0])
# sensitivity[f"Biochar price (${(base+0.02)*1000}/tonne)"] = [low, high]
# tea.system.flowsheet.stream["Biochar"].price = base

base = tea.biomass_collection_cost
tea.biomass_collection_cost = base * 0.8
low = tea.solve_price(oil_storage.outs[0])
tea.biomass_collection_cost = base * 1.2
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"Biomass collection cost (${base}/tonne)"] = [high, low]
tea.biomass_collection_cost = base
print(tea.solve_price(oil_storage.outs[0]))

base = tea.biomass_delivery_cost
tea.biomass_delivery_cost = base * 0.8
low = tea.solve_price(oil_storage.outs[0])
tea.biomass_delivery_cost = base * 1.2
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"Biomass delivery cost (${base:.2f}/tonne/mile)"] = [high, low]
tea.biomass_delivery_cost = base
print(tea.solve_price(oil_storage.outs[0]))

base = tea.bio_oil_delivery_cost
tea.bio_oil_delivery_cost = base * 0.8
low = tea.solve_price(oil_storage.outs[0])
tea.bio_oil_delivery_cost = base * 1.2
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"Bio-oil delivery cost (${base:0.3f}/tonne/mile)"] = [high, low]
tea.bio_oil_delivery_cost = base
print(tea.solve_price(oil_storage.outs[0]))

base = tea.system.flowsheet.stream["Biomass"].get_flow("tonnes/day", "H2O")
tea.system.flowsheet.stream["Biomass"].set_flow(base * 0.8, "tonnes/day", "H2O")
tea.system.simulate()
low = tea.solve_price(oil_storage.outs[0])
tea.system.flowsheet.stream["Biomass"].set_flow(base * 1.2, "tonnes/day", "H2O")
tea.system.simulate()
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"Biomass moisture ({25}%)"] = [high, low]
tea.system.flowsheet.stream["Biomass"].set_flow(base, "tonnes/day", "H2O")
tea.system.simulate()
print(tea.solve_price(oil_storage.outs[0]))

base = pyrolysis_yields["LightEnds"]
pyrolysis_yields["LightEnds"] = base * 0.8
tea.system.simulate()
low = tea.solve_price(oil_storage.outs[0])
pyrolysis_yields["LightEnds"] = base * 1.2
tea.system.simulate()
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"LightEnds yield ({base} kg/kg)"] = [low, high]
pyrolysis_yields["LightEnds"] = base
tea.system.simulate()
print(tea.solve_price(oil_storage.outs[0]))

base = pyrolysis_yields["HeavyEnds"]
pyrolysis_yields["HeavyEnds"] = base * 0.8
tea.system.simulate()
low = tea.solve_price(oil_storage.outs[0])
pyrolysis_yields["HeavyEnds"] = base * 1.2
tea.system.simulate()
high = tea.solve_price(oil_storage.outs[0])
sensitivity[f"HeavyEnds yield ({base} kg/kg)"] = [low, high]
pyrolysis_yields["HeavyEnds"] = base
tea.system.simulate()
print(tea.solve_price(oil_storage.outs[0]))

# base = pyrolysis_yields["Biochar"]
# pyrolysis_yields["Biochar"] = base * 0.8
# tea.system.simulate()
# low = tea.solve_price(oil_storage.outs[0])
# pyrolysis_yields["Biochar"] = base * 1.2
# tea.system.simulate()
# high = tea.solve_price(oil_storage.outs[0])
# sensitivity[f"Biochar yield ({base} kg/kg)"] = [low, high]
# pyrolysis_yields["Biochar"] = base
# tea.system.simulate()
# print(tea.solve_price(oil_storage.outs[0]))


# %%
sensitivityPD = pd.DataFrame(sensitivity, index=["Low", "High"]) - baseline
sensitivityPD.sort_values("Low", axis=1, inplace=True, key=abs)
sensitivityPD.to_csv("plots/sensitivity.csv")

plt.figure(figsize=(10, 8))
left_bars = plt.barh(sensitivityPD.columns, sensitivityPD.iloc[0], color="C0")
right_bars = plt.barh(sensitivityPD.columns, sensitivityPD.iloc[1], color="C1")
plt.legend((left_bars[0], right_bars[0]), ("Pessimistic", "Optimistic"))
plt.xlabel("Bio-oil Minimum Selling Price (USD/MT)")

x_min = min(sensitivityPD.min())
x_max = max(sensitivityPD.max())
x_range = [x_min * 1.4 + (x_max * 1.4 - x_min * 1.4) * i / 7 for i in range(7)]
plt.xticks(
    [i for i in x_range], ["$" + str((i + baseline) * 1000)[0:4] for i in x_range]
)
plt.xlim(x_min * 1.8, x_max * 1.4)

# add labels to the bars
for i, v in enumerate(sensitivityPD.iloc[0]):
    plt.text(v + 0.005, i - 0.1, "$" + str((v + baseline) * 1000)[0:3])
for i, v in enumerate(sensitivityPD.iloc[1]):
    plt.text(v - 0.015, i - 0.1, "$" + str((v + baseline) * 1000)[0:3])


plt.savefig("sensitivity.png", bbox_inches="tight", dpi=300)
# %%
from steamdrawio import *


# shapes["furnace"] = ["mxgraph.pid.vessels.furnace", "80", "100"]
# draw(system)
#%%


# %%
# LCA
lca_stover = 4.2231 / 1000  # kg CO2-eq/kg
lca_biooil_transport = 12.7 / 1000  # kg CO2-eq/kg
lca_power = 0.48  # kg/kWh
lca_natural_gas = 0.055/0.68 / 0.0283168  # kg CO2-eq/ft^3 / kg/Sm3 * ft3/m3 https://www.epa.gov/energy/greenhouse-gases-equivalencies-calculator-calculations-and-references.
biochar_carbon_permanence = 0.7
lca_biooil_transport = 12.7 / 1000  # kg CO2-eq/kg


def get_lca():

    oil_output = tea.system.flowsheet.stream["Bio_oil"].get_total_flow("kg/year")
    oil_carbon = (
        tea.system.flowsheet.stream["Bio_oil"].get_atomic_flow("C")
        * 12
        / tea.system.flowsheet.stream["Bio_oil"].get_total_flow("kg/hr")
    )
    char_output = tea.system.flowsheet.stream["Biochar"].get_total_flow("kg/year")
    char_carbon = (
        tea.system.flowsheet.stream["Biochar"].get_atomic_flow("C")
        * 12
        / tea.system.flowsheet.stream["Biochar"].get_total_flow("kg/hr")
    )

    ng_input = tea.system.flowsheet.stream["Dryer_Natural_gas"].get_total_flow("kg/year")
    ng_carbon = ng_input * lca_natural_gas

    power_carbon = (
        tea.system.flowsheet.stream["Biomass"].get_total_flow("tonnes/year")
        * grinder_power
        * lca_power
    )

    lca_results = {
        "Stover": lca_stover * tea.system.flowsheet.stream["Biomass"].get_total_flow("kg/year") / oil_output,
        "Biochar": -char_carbon
        * char_output
        * 44
        / 12
        * biochar_carbon_permanence
        / oil_output,
        "Power": power_carbon / oil_output,
        "Natural Gas": ng_carbon / oil_output,
        "Bio-oil": -oil_carbon * oil_output * 44 / 12 / oil_output,
        "Bio-oil transport": -oil_output * lca_biooil_transport / oil_output,
    }

    msp_lca = tea.solve_price(oil_storage.outs[0])
    carbon_abatement_cost = (msp_lca * oil_output) / (
        abs(sum(list(lca_results.values())) / 1000) * oil_output
    )

    lca_results_no_char = lca_results.copy()
    lca_results_no_char["Biochar"] = 0

    carbon_abatement_cost_no_char = (msp_lca * oil_output) / (
        abs(sum(list(lca_results_no_char.values())) / 1000) * oil_output
    )

    return (carbon_abatement_cost, carbon_abatement_cost_no_char)

print("Carbon abatement cost is ${:.3f}/kg CO2-eq".format(get_lca()[0]))
print("Carbon abatement cost without biochar is ${:.3f}/kg CO2-eq".format(get_lca()[1]))

# %%
# lca_sensitivity Analysis
lca_sensitivity = {}

tea.system.simulate()
baseline = get_lca()[1]

base = tea.system.flowsheet.stream["Biomass"].get_flow("tonnes/day", "H2O")
tea.system.flowsheet.stream["Biomass"].set_flow(base * 0.8, "tonnes/day", "H2O")
tea.system.simulate()
low = get_lca()[1]
tea.system.flowsheet.stream["Biomass"].set_flow(base * 1.2, "tonnes/day", "H2O")
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"Biomass moisture ({40}%)"] = [high, low]
tea.system.flowsheet.stream["Biomass"].set_flow(base, "tonnes/day", "H2O")
tea.system.simulate()

base = biochar_yield
pyrolysis_yields["Biochar"] = base * 0.8
tea.system.simulate()
low = get_lca()[1]
pyrolysis_yields["Biochar"] = base * 1.2
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"Biochar yield ({base} kg/kg)"] = [low, high]
pyrolysis_yields["Biochar"] = base
tea.system.simulate()

base = light_end_yield
pyrolysis_yields["LightEnds"] = base * 0.8
tea.system.simulate()
low = get_lca()[1]
pyrolysis_yields["LightEnds"] = base * 1.2
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"LightEnds yield ({base} kg/kg)"] = [low, high]
pyrolysis_yields["LightEnds"] = base
tea.system.simulate()

base = heavy_end_yield
pyrolysis_yields["HeavyEnds"] = base * 0.8
tea.system.simulate()
low = get_lca()[1]
pyrolysis_yields["HeavyEnds"] = base * 1.2
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"HeavyEnds yield ({base} kg/kg)"] = [low, high]
pyrolysis_yields["HeavyEnds"] = base
tea.system.simulate()

base = lca_stover
lca_stover = base * 0.8
tea.system.simulate()
low = get_lca()[1]
lca_stover = base * 1.2
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"Stover LCA ({base} kg CO2-eq/kg)"] = [low, high]
lca_stover = base

base = lca_biooil_transport
lca_biooil_transport = base * 0.8
tea.system.simulate()
low = get_lca()[1]
lca_biooil_transport = base * 1.2
tea.system.simulate()
high = get_lca()[1]
lca_sensitivity[f"Bio-oil transport LCA ({base} kg CO2-eq/kg)"] = [low, high]
lca_biooil_transport = base


# %%
lca_sensitivityPD = pd.DataFrame(lca_sensitivity, index=["Low", "High"]) - baseline
lca_sensitivityPD.sort_values("Low", axis=1, inplace=True, key=abs)

plt.figure(figsize=(10, 8))
left_bars = plt.barh(lca_sensitivityPD.columns, lca_sensitivityPD.iloc[0], color="C0")
right_bars = plt.barh(lca_sensitivityPD.columns, lca_sensitivityPD.iloc[1], color="C1")
plt.legend((left_bars[0], right_bars[0]), ("Optimistic", "Pessimistic"))
plt.xlabel("Marginal Abatement Cost (USD/MT CO2eq)")

x_min = min(lca_sensitivityPD.min())
x_max = max(lca_sensitivityPD.max())
x_range = [x_min * 1.4 + (x_max * 1.4 - x_min * 1.4) * i / 7 for i in range(7)]
plt.xticks([i for i in x_range], [str((i + baseline))[0:4] for i in x_range])
plt.xlim(x_min * 1.8, x_max * 1.4)

# add labels to the bars
for i, v in enumerate(lca_sensitivityPD.iloc[0]):
    plt.text(v + 2, i - 0.1, str((v + baseline))[0:3])
for i, v in enumerate(lca_sensitivityPD.iloc[1]):
    plt.text(v - 5, i - 0.1, str((v + baseline))[0:3])


plt.savefig("lca_sensitivity.png", bbox_inches="tight", dpi=300)

#%% OpenLCA lca analysis
steps = [
    "stover production",
    "stover collection",
    "stover transport",
    "stover storage",
    "Fast pyrolysis",
    "Bio-oil transportation",
    "Bio-oil well injection",
]

processes = []
for s in steps:
    print(s)
    r = recommend_process(s, 10)
    for k in r:
        print(k)
    processes.append(r)
    print("**************")

#%%
recommend_process("us electricity grid mixture", 20)
#%%
sesalca.get_process_quantitative_flow('electricity, high voltage, production mix | electricity, high voltage | APOS, U')


# %%
msps = {}
macs = {}
mspTables = {}
# Pioneer plant analysis

tea.cost_growth = 0.29
tea.pioneer = 0.2231
tea.learning_rate = 0
tea.units_manufactured = 1
msps["Pioneer Growth (29%)\nPerformance(22%)"] = (
    tea.solve_price(oil_storage.outs[0]) * 1000
)
macs["Pioneer Growth (29%)\nPerformance(22%)"] = get_lca()[1]
mspTables["Pioneer Growth (29%)\nPerformance(22%)"] = tea.mfsp_table(
    oil_storage.outs[0]
)

tea.cost_growth = 1
tea.pioneer = 0
tea.learning_rate = 0.20
tea.units_manufactured = 1000
msps["Learning Rate(20%)\nUnits (1000)"] = tea.solve_price(oil_storage.outs[0]) * 1000
macs["Learning Rate(20%)\nUnits (1000)"] = get_lca()[1]
mspTables["Learning Rate(20%)\nUnits (1000)"] = tea.mfsp_table(oil_storage.outs[0])

tea.cost_growth = 0.29
tea.pioneer = 0.2231
tea.learning_rate = 0.20
tea.units_manufactured = 1000
msps["Pioneer Growth (29%)\nPerformance(22%)\nLearning Rate(20%)\nUnits (1000)"] = (
    tea.solve_price(oil_storage.outs[0]) * 1000
)
macs[
    "Pioneer Growth (29%)\nPerformance(22%)\nLearning Rate(20%)\nUnits (1000)"
] = get_lca()[1]
mspTables[
    "Pioneer Growth (29%)\nPerformance(22%)\nLearning Rate(20%)\nUnits (1000)"
] = tea.mfsp_table(oil_storage.outs[0])

tea.cost_growth = 1
tea.pioneer = 1
tea.learning_rate = 0
tea.units_manufactured = 1
msps["Basecase"] = tea.solve_price(oil_storage.outs[0]) * 1000
macs["Basecase"] = get_lca()[1]
mspTables["Basecase"] = tea.mfsp_table(oil_storage.outs[0])

# %%
mspsPD = pd.DataFrame(msps, index=["MSP"]).transpose().sort_values("MSP")
mspsPD.plot.bar(
    legend=False,
    ylabel="Bio-oil Minimum Selling Price (MSP)\n($/tonne)",
    figsize=(8, 8),
)
mspsPD.to_csv("plots/msps.csv")
plt.savefig("msps_comparison.png", bbox_inches="tight", dpi=300)
# %%

macsPD = pd.DataFrame(macs, index=["MAC"]).transpose().sort_values("MAC")
macsPD.plot.bar(
    legend=False,
    ylabel="Marginal Abatement Cost (MAC)\n($/tonne CO2 removed)",
    figsize=(8, 8),
)


plt.savefig("mac_comparison.png", bbox_inches="tight", dpi=300)
# %%
mspTablesPD = pd.DataFrame(
    mspTables, index=mspTables["Basecase"].keys()
).transpose().sort_values(by="ROI") / tea.system.flowsheet.stream[
    "Bio_oil"
].get_total_flow(
    "tonnes/year"
)

mspTablesPD.plot.bar(
    stacked=True, figsize=(14, 10), ylabel="Bio-oil Minimum Selling Price ($/tonne)"
)
mspTablesPD.to_csv("plots/msps_table.csv") 
plt.savefig("msps_table_comparison.png", bbox_inches="tight", dpi=300)

# %%

plot_msp("pioneer_msp.png")

# %%
# Learning Rate Analysis

# %%

# system.save_report("charm_spreadsheet.xlsx")

#%%
# Summary Table 
biomass_production_stover_co2 = 44/12*12*biomass.get_atomic_flow("C")/1000*system.operating_hours/biomass.get_total_flow("tonnes/year")
biomass_production_emissions = biomass.get_total_flow("tonnes/year") * lca_stover/biomass.get_total_flow("tonnes/year")

stover_harvest_energy = 470 # MJ diesel/ha https://link.springer.com/article/10.1007/s12155-015-9664-4/tables/2
diesel_emissions = 74.14 * 9.488e-4 # kg CO2/MMBTU * MMBTU/MJ https://www.eia.gov/environment/emissions/co2_vol_mass.php
stover_production_yield = 5.4 # tons/acre
stover_harvest_area = capacity * system.operating_hours/24/stover_production_yield*0.405 #  (hectares) assumes 5.4 tons of corn stover per acre yield
biomass_collection_emissions = stover_harvest_energy*diesel_emissions/1000*stover_harvest_area/biomass.get_total_flow("tonnes/year")

stover_transport_energy = 70 # MJ/t 
stover_transport_emissions = stover_transport_energy*diesel_emissions/1000*biomass.get_total_flow("tonnes/year")/biomass.get_total_flow("tonnes/year")

bio_oil_co2 = oil_storage.outs[0].get_atomic_flow("C")*12*44/12/1000*system.operating_hours /biomass.get_total_flow("tonnes/year")
biochar_co2 = tea.system.flowsheet.stream["Biochar"].get_atomic_flow("C")*12*44/12/1000*system.operating_hours/biomass.get_total_flow("tonnes/year")

ng_emissions = tea.system.flowsheet.stream["Dryer_Natural_gas"].get_total_flow("tonnes/year")*lca_natural_gas/biomass.get_total_flow("tonnes/year")

bio_oil_transport_emissions = oil_storage.outs[0].get_total_flow("tonnes/year")*lca_biooil_transport/1000/biomass.get_total_flow("tonnes/year")

well_water_demand = 4 # million gallons
well_water_demand_kg = well_water_demand*3.78541e3 # tonnes

bio_oil_well_mud_fraction = 0.05 #13000 # kg/well https://content.cld.iop.org/journals/1748-9326/6/3/034014/revision1/erl381437suppdata.pdf?Expires=1708643500&Signature=VQo8jMC~yCF433XfDUmITlS8ZpHm8XYy6DuN1DGQgGpMTru1jd3r1n~mGSV8I84DweEJox6T7YRwRAE6wNo215FxYUKYxAU0GKMj0c7KtqTMYl19OJeK4EVwVZuHhJCFz~JHaBOw5t6xEsLZllO4z1D8QdizvgB-EQdU0txlVjtrBgXtqh1Ks9tjXbfPZVlzh4ptrwtHR8b1NbLFbe6oAL9Hn~qDMU6W01dBjWYeduZWK6pa7V8iV8ExAkZzsk-5eCgPn8k2FsNRigl2Cx1QfeW02x33c9TAsUvoJwnlpbokefWqvYmtHI50L-plnLcdbEvwhzVAgx34EWFZniCRrg__&Key-Pair-Id=KL1D8TIY3N7T8

bio_oil_well_demand = well_water_demand*bio_oil_well_mud_fraction # tonnes

wells_supplied_with_bio_oil = oil_storage.outs[0].get_total_flow("tonnes/year")/(well_water_demand_kg*bio_oil_well_mud_fraction)

well_hydraulic_fractioning_emissions = 370 # tonnes CO2-eq/well https://content.cld.iop.org/journals/1748-9326/6/3/034014/revision1/erl381437suppdata.pdf?Expires=1708643500&Signature=VQo8jMC~yCF433XfDUmITlS8ZpHm8XYy6DuN1DGQgGpMTru1jd3r1n~mGSV8I84DweEJox6T7YRwRAE6wNo215FxYUKYxAU0GKMj0c7KtqTMYl19OJeK4EVwVZuHhJCFz~JHaBOw5t6xEsLZllO4z1D8QdizvgB-EQdU0txlVjtrBgXtqh1Ks9tjXbfPZVlzh4ptrwtHR8b1NbLFbe6oAL9Hn~qDMU6W01dBjWYeduZWK6pa7V8iV8ExAkZzsk-5eCgPn8k2FsNRigl2Cx1QfeW02x33c9TAsUvoJwnlpbokefWqvYmtHI50L-plnLcdbEvwhzVAgx34EWFZniCRrg__&Key-Pair-Id=KL1D8TIY3N7T8

bio_oil_injection_cost_emissions = wells_supplied_with_bio_oil*well_hydraulic_fractioning_emissions/biomass.get_total_flow("tonnes/year")

all_emissions = [
    -biomass_production_emissions,
    -biomass_collection_emissions,
    -stover_transport_emissions,
    bio_oil_co2,
    biochar_co2*(biochar_carbon_permanence),
    -ng_emissions,
    -bio_oil_transport_emissions,
    -bio_oil_injection_cost_emissions
]

all_emissions_formatted = {
    "Biomass (Embodied)": biomass_production_stover_co2,
    "Biomass Production": -biomass_production_emissions,
    "Biomass Collection": -biomass_collection_emissions,
    "Stover Transport": -stover_transport_emissions,
    "Bio-oil CO2 (Embodied)": bio_oil_co2,
    "Biochar CO2 (Embodied)": biochar_co2,
    "Biochar CO2 (Embodied; 70% Permanent)": biochar_co2*(biochar_carbon_permanence),
    "Natural Gas": -ng_emissions,
    "Biogenic Emitted": biomass_production_stover_co2 - bio_oil_co2 - biochar_co2,
    "Bio-oil Transport": -bio_oil_transport_emissions,
    "Bio-oil Injection": -bio_oil_injection_cost_emissions
}

net_emissions = sum(all_emissions)
# %%
biomass_production_cash = 0 # 20 # $/ton
biomass_collection_cash = 10 # $/ton
biomass_transport_cash = 0 # $/ton

costs = tea.mfsp_table(oil_storage.outs[0])

bio_oil_cash_biomass = (sum(list(costs.values()))-costs["Bio-Oil Delivery"])/biomass.get_total_flow("tonnes/year")
bio_oil_cash_no_char = (sum(list(costs.values()))-costs["Bio-Oil Delivery"])/oil_storage.outs[0].get_total_flow("tonnes/year")
bio_oil_cash_with_char = (sum(list(costs.values()))-costs["Bio-Oil Delivery"]-tea.system.flowsheet.stream["Biochar"].get_total_flow("tonnes/year")*100)/oil_storage.outs[0].get_total_flow("tonnes/year")

bio_oil_transport = costs["Bio-Oil Delivery"]/oil_storage.outs[0].get_total_flow("tonnes/year")

water_injection_cost = 0.6 # $/gallon
# water_injection_cost = 7.4 # $/tonne Henderson, Claire, Acharya, Harish, Matis, Hope, Kommepalli, Hareesh, Moore, Brian, and Wang, Hua. Cost Effective Recovery of Low-TDS Frac Flowback Water for Re-use. United States: N. p., 2011. Web. doi:10.2172/1030557.

bio_oil_injection_cost = water_injection_cost

oil_biomass_ratio = oil_storage.outs[0].get_total_flow("tonnes/year")/biomass.get_total_flow("tonnes/year")

net_costs = (biomass_production_cash + biomass_collection_cash + biomass_transport_cash + bio_oil_cash_biomass + bio_oil_transport*oil_biomass_ratio + bio_oil_injection_cost*oil_biomass_ratio)

all_costs_formatted = {
    "Biomass Production": biomass_production_cash,
    "Biomass Collection": biomass_collection_cash,
    "Biomass Transport": biomass_transport_cash,
    "Bio-oil Cash ($/t biomass)": bio_oil_cash_biomass,
    "Bio-oil Cash (No Biochar $/t oil)": bio_oil_cash_no_char,
    "Bio-oil Cash (With Biochar $/t oil)": bio_oil_cash_with_char,
    "Bio-oil Transport": bio_oil_transport,
    "Bio-oil Injection": bio_oil_injection_cost,
    "Net Costs": net_costs
}
all_costs_formatted
# %%
net_costs_per_co2 = net_costs/net_emissions
net_costs_per_co2
# %%
net_costs_no_injection = (biomass_production_cash + biomass_collection_cash + biomass_transport_cash + bio_oil_cash_biomass + bio_oil_transport*oil_biomass_ratio)
net_emissions_no_injection = sum(all_emissions) + bio_oil_injection_cost_emissions
#%%
net_costs_per_co2_no_injection = net_costs_no_injection/net_emissions_no_injection
net_costs_per_co2_no_injection
# %%

summary= [
        ["Carbon Flows", "", ""],
        ["Biomass (Embodied)", "Emitted", ""],
        [biomass_production_stover_co2, -biomass_production_emissions, ""],
        ["Biomass Collection", "", ""],
        [biomass_collection_emissions, "", ""],
        ["Stover Transport", "", ""],
        [stover_transport_emissions, "", ""],
        ["Bio-oil CO2 (Embodied)", "Biochar CO2 (Embodied)", "Emitted"],
        [bio_oil_co2, biochar_co2, -ng_emissions],
        ["Emitted", "", ""],
        [-bio_oil_transport_emissions, "", ""],
        ["Emitted", "", ""],
        [-bio_oil_injection_cost_emissions, "", ""],
        [net_emissions, net_emissions_no_injection, ""],
        ["Cash Flows", "", ""],
        [biomass_production_cash, "", ""],
        [biomass_collection_cash, "", ""],
        [biomass_transport_cash, "", ""],
        [bio_oil_cash_biomass, bio_oil_cash_no_char, bio_oil_cash_with_char],
        [bio_oil_transport*oil_biomass_ratio, "", ""],
        [bio_oil_injection_cost*oil_biomass_ratio, "", ""],
        [net_costs, net_costs_no_injection, ""],
        [net_costs_per_co2, net_costs_per_co2_no_injection, ""],
]

tsv = ""
for row in summary:
    tsv += "\t".join([str(i) for i in row]) + "\n"
print(tsv)

# %%
# globals()
# %%
# import numpy 
# for key, value in globals().items():
#     if type(value) == numpy.float64 or type(value) == float:
#         print(key, value)
# %%
#%%
# Summary Table  2
tea.learning_rate = 0.2
tea.units_manufactured = 1000
msp = tea.solve_price(oil_storage.outs[0])
oil_storage.outs[0].price = msp        

biomass_production_stover_co2 = 44/12*12*biomass.get_atomic_flow("C")/1000*system.operating_hours/biomass.get_total_flow("tonnes/year")
biomass_production_emissions = biomass.get_total_flow("tonnes/year") * lca_stover/biomass.get_total_flow("tonnes/year")

stover_harvest_energy = 470 # MJ diesel/ha https://link.springer.com/article/10.1007/s12155-015-9664-4/tables/2
diesel_emissions = 74.14 * 9.488e-4 # kg CO2/MMBTU * MMBTU/MJ https://www.eia.gov/environment/emissions/co2_vol_mass.php
stover_production_yield = 5.4 # tons/acre
stover_harvest_area = capacity * system.operating_hours/24/stover_production_yield*0.405 #  (hectares) assumes 5.4 tons of corn stover per acre yield
biomass_collection_emissions = stover_harvest_energy*diesel_emissions/1000*stover_harvest_area/biomass.get_total_flow("tonnes/year")

stover_transport_energy = 70 # MJ/t 
stover_transport_emissions = stover_transport_energy*diesel_emissions/1000*biomass.get_total_flow("tonnes/year")/biomass.get_total_flow("tonnes/year")

bio_oil_co2 = oil_storage.outs[0].get_atomic_flow("C")*12*44/12/1000*system.operating_hours /biomass.get_total_flow("tonnes/year")
biochar_co2 = tea.system.flowsheet.stream["Biochar"].get_atomic_flow("C")*12*44/12/1000*system.operating_hours/biomass.get_total_flow("tonnes/year")

ng_emissions = tea.system.flowsheet.stream["Dryer_Natural_gas"].get_total_flow("tonnes/year")*lca_natural_gas/biomass.get_total_flow("tonnes/year")

bio_oil_transport_emissions = oil_storage.outs[0].get_total_flow("tonnes/year")*lca_biooil_transport/1000/biomass.get_total_flow("tonnes/year")

well_water_demand = 4 # million gallons
well_water_demand_kg = well_water_demand*3.78541e3 # tonnes

bio_oil_well_mud_fraction = 0.05 #13000 # kg/well https://content.cld.iop.org/journals/1748-9326/6/3/034014/revision1/erl381437suppdata.pdf?Expires=1708643500&Signature=VQo8jMC~yCF433XfDUmITlS8ZpHm8XYy6DuN1DGQgGpMTru1jd3r1n~mGSV8I84DweEJox6T7YRwRAE6wNo215FxYUKYxAU0GKMj0c7KtqTMYl19OJeK4EVwVZuHhJCFz~JHaBOw5t6xEsLZllO4z1D8QdizvgB-EQdU0txlVjtrBgXtqh1Ks9tjXbfPZVlzh4ptrwtHR8b1NbLFbe6oAL9Hn~qDMU6W01dBjWYeduZWK6pa7V8iV8ExAkZzsk-5eCgPn8k2FsNRigl2Cx1QfeW02x33c9TAsUvoJwnlpbokefWqvYmtHI50L-plnLcdbEvwhzVAgx34EWFZniCRrg__&Key-Pair-Id=KL1D8TIY3N7T8

bio_oil_well_demand = well_water_demand*bio_oil_well_mud_fraction # tonnes

wells_supplied_with_bio_oil = oil_storage.outs[0].get_total_flow("tonnes/year")/(well_water_demand_kg*bio_oil_well_mud_fraction)

well_hydraulic_fractioning_emissions = 370 # tonnes CO2-eq/well https://content.cld.iop.org/journals/1748-9326/6/3/034014/revision1/erl381437suppdata.pdf?Expires=1708643500&Signature=VQo8jMC~yCF433XfDUmITlS8ZpHm8XYy6DuN1DGQgGpMTru1jd3r1n~mGSV8I84DweEJox6T7YRwRAE6wNo215FxYUKYxAU0GKMj0c7KtqTMYl19OJeK4EVwVZuHhJCFz~JHaBOw5t6xEsLZllO4z1D8QdizvgB-EQdU0txlVjtrBgXtqh1Ks9tjXbfPZVlzh4ptrwtHR8b1NbLFbe6oAL9Hn~qDMU6W01dBjWYeduZWK6pa7V8iV8ExAkZzsk-5eCgPn8k2FsNRigl2Cx1QfeW02x33c9TAsUvoJwnlpbokefWqvYmtHI50L-plnLcdbEvwhzVAgx34EWFZniCRrg__&Key-Pair-Id=KL1D8TIY3N7T8


bio_oil_injection_cost_emissions = wells_supplied_with_bio_oil*well_hydraulic_fractioning_emissions/biomass.get_total_flow("tonnes/year")

all_emissions = [
    -biomass_production_emissions,
    -biomass_collection_emissions,
    -stover_transport_emissions,
    bio_oil_co2,
    biochar_co2*(biochar_carbon_permanence),
    -ng_emissions,
    -bio_oil_transport_emissions,
    -bio_oil_injection_cost_emissions
]

all_emissions_no_char = [
    -biomass_production_emissions,
    -biomass_collection_emissions,
    -stover_transport_emissions,
    bio_oil_co2,
    -ng_emissions,
    -bio_oil_transport_emissions,
    -bio_oil_injection_cost_emissions
]

all_emissions_formatted = {
    "Biomass (Embodied)": biomass_production_stover_co2,
    "Biomass Production": -biomass_production_emissions,
    "Biomass Collection": -biomass_collection_emissions,
    "Stover Transport": -stover_transport_emissions,
    "Bio-oil CO2 (Embodied)": bio_oil_co2,
    "Biochar CO2 (Embodied)": biochar_co2,
    "Biochar CO2 (Embodied; 70% Permanent)": biochar_co2*(biochar_carbon_permanence),
    "Natural Gas": -ng_emissions,
    "Biogenic Emitted": biomass_production_stover_co2 - bio_oil_co2 - biochar_co2,
    "Bio-oil Transport": -bio_oil_transport_emissions,
    "Bio-oil Injection": -bio_oil_injection_cost_emissions
}

net_emissions = sum(all_emissions)
# %%
biomass_production_cash = 0 # 20 # $/ton
biomass_collection_cash = 10 # $/ton
biomass_transport_cash = 0 # $/ton

costs = tea.mfsp_table(oil_storage.outs[0])

bio_oil_cash_biomass = (sum(list(costs.values()))-costs["Bio-Oil Delivery"])/biomass.get_total_flow("tonnes/year")
bio_oil_cash_no_char = (sum(list(costs.values()))-costs["Bio-Oil Delivery"])/oil_storage.outs[0].get_total_flow("tonnes/year")
bio_oil_cash_with_char = (sum(list(costs.values()))-costs["Bio-Oil Delivery"]-tea.system.flowsheet.stream["Biochar"].get_total_flow("tonnes/year")*100)/oil_storage.outs[0].get_total_flow("tonnes/year")

bio_oil_transport = costs["Bio-Oil Delivery"]/oil_storage.outs[0].get_total_flow("tonnes/year")

water_injection_cost = 0.6 # $/gallon
bio_oil_injection_cost = water_injection_cost

oil_biomass_ratio = oil_storage.outs[0].get_total_flow("tonnes/year")/biomass.get_total_flow("tonnes/year")

net_costs = (biomass_production_cash + biomass_collection_cash + biomass_transport_cash + bio_oil_cash_biomass + bio_oil_transport*oil_biomass_ratio + bio_oil_injection_cost*oil_biomass_ratio)

all_costs_formatted = {
    "Biomass Production": biomass_production_cash,
    "Biomass Collection": biomass_collection_cash,
    "Biomass Transport": biomass_transport_cash,
    "Bio-oil Cash ($/t biomass)": bio_oil_cash_biomass,
    "Bio-oil Cash (No Biochar $/t oil)": bio_oil_cash_no_char,
    "Bio-oil Cash (With Biochar $/t oil)": bio_oil_cash_with_char,
    "Bio-oil Transport": bio_oil_transport,
    "Bio-oil Injection": bio_oil_injection_cost,
    "Net Costs": net_costs
}
all_costs_formatted
# %%
net_costs_per_co2 = net_costs/net_emissions
net_costs_per_co2
# %%
net_costs_no_injection = (biomass_production_cash + biomass_collection_cash + biomass_transport_cash + bio_oil_cash_biomass + bio_oil_transport*oil_biomass_ratio)
net_emissions_no_injection = sum(all_emissions) + bio_oil_injection_cost_emissions
#%%
net_costs_per_co2_no_injection = net_costs_no_injection/net_emissions_no_injection
net_costs_per_co2_no_injection
# %%

summary= [
        ["Carbon Flows", "", ""],
        ["Biomass (Embodied)", "Emitted", ""],
        [biomass_production_stover_co2, -biomass_production_emissions, ""],
        ["Biomass Collection", "", ""],
        [biomass_collection_emissions, "", ""],
        ["Stover Transport", "", ""],
        [stover_transport_emissions, "", ""],
        ["Bio-oil CO2 (Embodied)", "Biochar CO2 (Embodied)", "Emitted"],
        [bio_oil_co2, biochar_co2, -ng_emissions],
        ["Emitted", "", ""],
        [-bio_oil_transport_emissions, "", ""],
        ["Emitted", "", ""],
        [-bio_oil_injection_cost_emissions, "", ""],
        [net_emissions, net_emissions_no_injection, ""],
        ["Cash Flows", "", ""],
        [biomass_production_cash, "", ""],
        [biomass_collection_cash, "", ""],
        [biomass_transport_cash, "", ""],
        [bio_oil_cash_biomass, bio_oil_cash_no_char, bio_oil_cash_with_char],
        [bio_oil_transport*oil_biomass_ratio, "", ""],
        [bio_oil_injection_cost*oil_biomass_ratio, "", ""],
        [net_costs, net_costs_no_injection, ""],
        [net_costs_per_co2, net_costs_per_co2_no_injection, ""],
]

tsv = ""
for row in summary:
    tsv += "\t".join([str(i) for i in row]) + "\n"
print(tsv)

tea.learning_rate = 0  
tea.units_manufactured = 1
# %%
# globals()
# %%
# import numpy 
# for key, value in globals().items():
#     if type(value) == numpy.float64 or type(value) == float:
#         print(key, value)
# %%
from steamdrawio import draw 
def labeling(stream):
    return f"{stream.ID}\nM: {stream.F_mass*24/1000:.2f} t/d\nC (%): {stream.get_atomic_flow('C')*12/stream.F_mass:.2f}"
draw(system, label_fn=labeling)
# %% 
?draw
# %%
for unit in system.units:
    if unit.utility_cost != None and unit.utility_cost*tea.operating_days*24/system.utility_cost > 0.01:
        print(unit.ID, unit.utility_cost*tea.operating_days*24/system.utility_cost)
    elif unit.utility_cost != None:
        print("Other")
# %%
