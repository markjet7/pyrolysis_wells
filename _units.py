from typing import Optional
import biosteam as bst
import thermosteam as tmo

from biosteam.units.decorators import cost


class Reactor(bst.Unit):
    def __init__(self, ID: str | None = "", ins=None, outs=None, thermo=None):
        super().__init__(ID, ins, outs, thermo)

    def _run(self):
        super()._run()

    def _design(self):
        feed = self.ins[0]
        output = self.outs[0]

        heat_duty = output.H - feed.H
        if heat_duty > 0:
            self.add_heat_utility(heat_duty, feed.T, output.T)
        else:
            self.add_heat_utility(0, feed.T, output.T)

    # https://www.nrel.gov/docs/fy15osti/62455.pdf
    def _cost(self):
        super()._cost()
        self.purchase_costs["Reactors"] = (
            2 * 3818000 * (self.ins[0].get_total_flow("lb/hr") / 2526000) ** 0.5
        )


class PyrolysisCompressor(bst.units.Compressor, new_graphics=False):
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.S = feed.S

        T2s = feed.T * (out.P / feed.P) ** ((1.4 - 1) / 1.4)

        # eta = (h2s - h1)/(h2 - h1) = (outs - feed)/(out - feed)
        # h2 = h1 + eta*(h2s - h1)
        # T2 = T1 + eta*(T2s - T1)
        out.T = feed.T + self.eta * (T2s - feed.T)
        if self.vle is True:
            out.vle(T=out.T, P=out.P)
        self.ideal_power, self.ideal_duty = self._calculate_ideal_power_and_duty()
        self.design_results["Ideal power"] = self.ideal_power / 3600.0  # kW
        self.design_results["Ideal duty"] = 0.0

    def _design(self):
        super()._design()
        self._set_power(self.design_results["Ideal power"] / self.eta)
        self.design_results["Compressors in parallel"] = 2


@cost(
    "Flow rate",
    "Hot process water softener system",
    CE=551,
    cost=78e3,
    S=235803,
    n=0.6,
    BM=1.8,
)
@cost("Flow rate", "Amine addition pkg", CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost("Flow rate", "Deaerator", CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost("Flow rate", "Boiler", CE=551, cost=28550e3, kW=1371, S=238686, n=0.6, BM=1.8)
@cost(
    "Ash disposal", "Baghouse bags", CE=551, cost=466183.0 / 4363.0, n=1.0, lifetime=5
)
class Boiler(bst.Unit):
    network_priority = 1
    _N_ins = 3
    _N_outs = 2
    _units = {"Flow rate": "kg/hr", "Work": "kW", "Ash disposal": "kg/hr"}

    def __init__(
        self,
        ID: str | None = "",
        ins=None,
        outs=None,
        thermo=None,
        boiler_efficiency=0.95,
        natural_gas_price=0.218,
        ash_disposal_price=0.010,
        other_units = None
    ):
        super().__init__(ID, ins, outs, thermo)
        self.boiler_efficiency = boiler_efficiency
        self.natural_gas_price = natural_gas_price
        self.ash_disposal_price = ash_disposal_price
        self.agent = bst.settings.get_heating_agent("low_pressure_steam")
        # self.natural_gas = self.ins[1]
        # self.ash_disposal = self.outs[1]
        # self.define_utility('Natural Gas', self.natural_gas)
        # self.define_utility("Ash Disposal", self.ash_disposal)
        self.other_units = other_units
        self.steam_utilities = set()
        self.steam_demand = self.agent.to_stream()
        self.ins[0].ID = "Boiler_fuel"
        self.ins[1].ID = "Natural_gas"
        self.ins[2].ID = "Boiler_air"

    def _run(self):
        self.ins[0].ID = "Boiler_fuel"
        self.ins[1].ID = "Natural_gas"
        self.ins[2].ID = "Boiler_air"

    def _cost(self):
        self._decorated_cost()

    def _design(self):
        units = self.other_units
        for u in units:
            for hu in u.heat_utilities:
                self.steam_utilities.add(hu)

        self.heat_demand = sum([i.duty for i in self.steam_utilities])

        fuel_needed = (
            self.heat_demand / self.boiler_efficiency
            - self.ins[0].HHV * self.boiler_efficiency
        )
        if fuel_needed > 0 and "CH4" in self.chemicals:
            natural_gas = bst.Stream("natural_gas", CH4=1, units="kmol/hr")

            natural_gas.set_total_flow(fuel_needed / natural_gas.HHV, units="kmol/hr")
            self.ins[1] = natural_gas

        combustion_rxns = self.chemicals.get_combustion_reactions()
        emissions = self.outs[0]
        carbon = 0
        hydrogen = 0
        oxygen = 0

        for f in self.ins:
            carbon += f.get_atomic_flow("C")
            hydrogen += f.get_atomic_flow("H")
            oxygen += f.get_atomic_flow("O")

        self.ins[2].imol["O2"] = 2 * carbon + 1 / 2 * hydrogen - 1 / 2 * oxygen
        self.ins[2].imol["N2"] = 3.76 * self.ins[2].imol["O2"]

        for f in self.ins:
            emissions.mol[:] += f.mol

        combustion_rxns.force_reaction(emissions)
        emissions.imol["O2"] = 0
        H_combustion = 0
        for f in self.ins:
            H_combustion += f.HHV

        self.H_content = self.boiler_efficiency * H_combustion

        duty_over_mol = 39000
        self.total_steam = self.H_content / duty_over_mol

        for c in self.chemicals:
            if c.get_phase() == "l":
                self.outs[1].imol[c.ID] = self.outs[0].imol[c.ID]
                self.outs[0].imol[c.ID] = 0
            elif c.get_phase() == "s":
                self.outs[1].imass[c.ID] = self.outs[0].imass[c.ID]
                self.outs[0].imass[c.ID] = 0

        self.design_results["Flow rate"] = self.total_steam
        self.design_results["Heat"] = self.H_content
        self.design_results["Ash disposal"] = self.outs[1].F_mass
