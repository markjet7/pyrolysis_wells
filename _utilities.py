import biosteam as bst
# Utilities
HeatUtility = bst.HeatUtility
Gas_utility = bst.UtilityAgent(
    "natural_gas",
    T=1200,
    P=101325,
    T_limit=1100,
    H2O=1,
    heat_transfer_efficiency=0.85,
)
# bst.settings.heating_agents.append(Gas_utility)
# HeatUtility.heating_agents.append(Gas_utility)
HeatUtility.default_heating_agents()
HeatUtility.default_cooling_agents()

Cooling_utility = HeatUtility.get_agent("chilled_water")
Cooling_utility.regeneration_price = 0
Cooling_utility.heat_transfer_price = 0

Steam_utility = HeatUtility.get_agent("low_pressure_steam")
Steam_utility.regeneration_price = Steam_utility.regeneration_price * 0.1
HeatUtility.heating_agents.append(Gas_utility)