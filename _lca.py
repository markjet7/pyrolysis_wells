#%%
import sys
sys.path.append("/System/Volumes/Data/Users/markmw/Library/CloudStorage/OneDrive-IowaStateUniversity/General - SESA Home/Scripts/sesalca/")

import sesalca

def lca(system, feedstock="stover", product="Biooil", method="TRACI"):
    processes = [
        {
            "Name": f"{feedstock} production",
            "Amount": system.flowsheet.stream["Biomass"].F_mass,
            "Unit": "kg",
            "Type": "Input",
            "Reference": False,
            "Provider": None,
            "Avoided": False
        },
        {
            "Name": f"{feedstock} transport",
            "Amount": system.flowsheet.stream["Biomass"].F_mass*32.19/1000,
            "Unit": "t*km",
            "Type": "Input",
            "Reference": False,
            "Provider": 'transport, tractor and trailer, agricultural | transport, tractor and trailer, agricultural | APOS, U',
            "Avoided": False
        },
        {
            "Name": "oil transport",
            "Amount": system.flowsheet.stream["Biooil"].F_mass*32.19/1000,
            "Unit": "t*km",
            "Type": "Input",
            "Reference": False,
            "Provider": 'transport, freight, lorry 28 metric ton, vegetable oil methyl ester 100% | transport, freight, lorry 28 metric ton, vegetable oil methyl ester 100% | APOS, U',
            "Avoided": False
        },
        {
            "Name": "deep well closure",
            "Amount": system.flowsheet.stream["Biooil"].F_mass/20, # well emissions are in equivalent kg of input per meter of well depth
            "Unit": "m",
            "Type": "Output",
            "Reference": False,
            "Provider": 'market for deep well closure | deep well closure | APOS, U',
            "Avoided": False
        }
    ]

    process = sesalca.biosteam_to_lca(system, product, processes)
    impacts = sesalca.get_total_impacts(process, method)
    flows = sesalca.get_direct_flow_impacts(system.ID,method)
    result = {
        "impacts": impacts,
        "flows": flows
    }

    print("LCA Result: ", result)
    return result

#%%