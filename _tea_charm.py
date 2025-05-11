import biosteam as bst
from numpy import *

def NPV_at_IRR(IRR, cashflow_array, duration_array):
    """Return NPV at given IRR and cashflow data."""
    return (cashflow_array/(1.+IRR)**duration_array).sum()

# Transportation 
# Flatbed truck dimensions
# Length: 48'
# Width: 8' 6"
# Height: 5'
# Load Height: 8' 6"
# Total Height Limit with Cargo: 13' 6"
# Weight Limit: 48,000 lbs

flat_bed_cost = 2.604 # $/mi
flat_bed_volume = 48 * 8.5 * 5 * 0.02832 # m3
flat_bed_weight = 48000/2205 # tonnes from 48000 lbs 
bio_oil_density = 1200 # kg/m3 https://task34.ieabioenergy.com/bio-oil/
biomass_density = 50 # kg/m3 https://www.mdpi.com/1996-1073/7/6/4019

class TEA:
    """
    Create a TEA object for techno-economic analysis of a biorefinery.
    """

    # __slots__ = (
    #     "labor_cost",
    #     "fringe_benefits",
    #     "maintenance",
    #     "property_tax",
    #     "property_insurance",
    #     "_FCI_cached",
    #     "labor_cost",
    #     "prices",
    #     "supplies",
    #     "maintanance",
    #     "administration",
    #     "biomass_collection_cost",
    #     "biomass_delivery_cost",
    #     "biomass_delivery_distance",
    #     "bio_oil_delivery_cost",
    #     "bio_oil_delivery_distance",
    #     "bio_oil_plug_value",
    #     "working_capital",
    #     "pioneer",
    #     "cost_growth",
    #     "learning_rate",
    #     "units_manufactured",
    # )

    def __init__(
        self,
        system,
        IRR=0,
        duration=(2018, 2038),
        depreciation="MACRS7",
        income_tax=0.21,
        operating_days=333,
        lang_factor=2.5,
        construction_schedule=(0.4, 0.6),
        WC_over_FCI=0.05,
        labor_cost=0,
        fringe_benefits=0.4,
        property_tax=0.001,
        property_insurance=0.005,
        supplies=0.2,
        maintenance=0.003,
        administration=0.005,
        finance_interest=0.4,
        finance_years=10,
        finance_fraction=0.07,
        startup_months=0,
        startup_FOCfrac=0,
        startup_VOCfrac=0,
        startup_salesfrac=0,
        working_capital=0.05,
        pioneer=1,
        cost_growth=1,
        learning_rate=0,
        manufactured_units=1,
    ):
        self.system = system
        self.IRR = IRR
        self.duration = duration
        self.depreciation = depreciation
        self.income_tax = income_tax
        self.operating_days = operating_days
        self.lang_factor = lang_factor
        self.construction_schedule = construction_schedule
        self.WC_over_FCI = WC_over_FCI
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies = supplies
        self.maintenance = maintenance
        self.administration = administration
        self.finance_fraction = finance_fraction
        self.finance_years = finance_years
        self.finance_interest = finance_interest
        self.startup_months = startup_months
        self.startup_FOCfrac = startup_FOCfrac
        self.startup_VOCfrac = startup_VOCfrac
        self.startup_salesfrac = startup_salesfrac
        self.working_capital = working_capital
        self.pioneer = pioneer
        self.cost_growth = cost_growth
        self.learning_rate = learning_rate
        self.units_manufactured = manufactured_units
        self.biomass_delivery_cost = 0.71  # $/dry-ton/mi
        self.bio_oil_plug_value = (
            0  # $/kg https://www.alibaba.com/showroom/oil-well-cement.html
        )
        capacity = self.system.flowsheet.stream["Dried_biomass"].F_mass*24/1000 # tpd
        if capacity < 100:
            self.biomass_collection_cost = 10 # $/dry-ton in-field collection and stacking
            self.biomass_delivery_distance = 0  # mi
        else:
            self.biomass_collection_cost = 40  # $/dry-ton includes windrowing, baling, in-field collection and stacking
            self.biomass_delivery_distance = 1.5*(capacity*self.operating_days/(5/(0.00156 *0.6*3.14)))  # mi
        self.biomass_delivery_cost = 2.604/20.9   # 0.71  # $/dry-ton/mi An Analysis of the Operational Costs of Trucking:2022 Update
        self.bio_oil_delivery_cost = 2.604/20.9
        # 2.5 / (
        #     2 * 20
        # )  # 0.71 * 0.0083 / 10  # $/dry-ton/mile times the ratio of biomass density (60 kg/m3) to bio-oil density (10 lb/gal Badger et al. 2011)
        # https://www.tcicapital.com/tci-insights/current-freight-trends/
        # https://ops.fhwa.dot.gov/freight/publications/size_regs_final_rpt/
        self.bio_oil_delivery_distance = 50  # mi
        self.bio_oil_plug_value = (
            0  # $/kg https://www.alibaba.com/showroom/oil-well-cement.html
        )
        self.bio_oil_injection_cost = 7.4 # $/tonne # 7.4 to 22.1 Henderson, Claire, Acharya, Harish, Matis, Hope, Kommepalli, Hareesh, Moore, Brian, and Wang, Hua. Cost Effective Recovery of Low-TDS Frac Flowback Water for Re-use. United States: N. p., 2011. Web. doi:10.2172/1030557.

        employee_costs = {
            "Plant Manager": [141569, 1],
            "Plant Engineer": [67414, 1],
            "Maintenance Supr": [54894, 1],
            "Maintenance Tech": [38522, 2],
            "Lab Manager": [53931, 1],
            "Lab Technician": [38522, 2],
            "Lab Tech-Enzyme": [38522, 0],
            "Shift Supervisor": [46227, 5],
            "Shift Operators": [38522, 5],
            "Shift Oper-Enzyme": [38522, 0],
            "Yard Employees": [26966, 12],
            "Clerks & Secretaries": [34670, 3],
        }  # NREL Reports

        self.labor_cost = (
            sum([salary * number for (salary, number) in employee_costs.values()])
            * system.flowsheet.stream["Biomass"].get_total_flow("tonnes/day")
            / 2000
        )

        self.prices = {
            "Feedstock": 0.03,  # $30 per tonne
            "Natural gas": 0.15,  # $3 per mmbtu
            "Electricity": 0.065,  # $0.1 per kWh
            "Water": 0.0005,  # $0.0005 per kg
            "Solid handling": 0.02,  # $0.02 per kg
            "Heavy oil": 0.05,  # $0.05 per kg
            "Light oil": 0.03,  # $0.03 per kg
            "Blend oil": 0,  # $0.03 per kg
            "Biochar": 0 / 1000,  # $0.1 per kg
        }

        self.system.flowsheet.stream["Biomass"].price = self.prices["Feedstock"]
        self.system.flowsheet.stream["Natural_gas"].price = self.prices["Natural gas"]
        self.system.flowsheet.unit["Pyrolysis_cyclone"].outs[0].price = self.prices["Biochar"] - self.prices["Solid handling"]

    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost  # installed_equipment_cost (generic number)

    def _TDC(self, DPI):
        return DPI

    def _FCI(self, TDC):
        self._FCI_cached = TDC
        return TDC

    def _FOC(self, FCI):
        plant_costs = FCI * (
            self.property_tax
            + self.property_insurance
            + self.maintenance
            + self.administration
        ) + self.labor_cost * (1 + self.fringe_benefits + self.supplies)
        logistic_costs = (
            self.biomass_collection_cost
            * self.system.flowsheet.stream["biomass"].get_total_flow("tonnes/year")
            + self.biomass_delivery_cost
            * self.biomass_delivery_distance
            * self.system.flowsheet.stream["biomass"].get_total_flow("tonnes/year")
            + self.bio_oil_delivery_cost
            * self.bio_oil_delivery_distance
            * self.system.flowsheet.stream["Bio_oil"].get_total_flow("tonnes/year")
        )
        return plant_costs + logistic_costs    

    def FCI(self):
        return (
            self.system.installed_equipment_cost
            / self.cost_growth
            * (self.units_manufactured) ** log2(1 - self.learning_rate)
        )


    def NPV(self):
        logistic_costs = (
            self.biomass_collection_cost
            * self.system.flowsheet.stream["Biomass"].get_total_flow("tonnes/year")
            + self.biomass_delivery_cost
            * self.biomass_delivery_distance
            * self.system.flowsheet.stream["Biomass"].get_total_flow("tonnes/year")
            + self.bio_oil_delivery_cost
            * self.bio_oil_delivery_distance
            * self.system.flowsheet.stream["Biooil"].get_total_flow("tonnes/year") 
            + self.bio_oil_injection_cost
            * self.system.flowsheet.stream["Biooil"].get_total_flow("tonnes/year")
        )

        material_cost = logistic_costs
        for feed in self.system.feeds:
            material_cost += feed.cost * self.operating_days * 24

        for unit in self.system.units:
            if unit.utility_cost != None:
                material_cost += unit.utility_cost * self.operating_days * 24
        # credit_cost = 0
        # for product in self.system.products:
        #     credit_cost -= product.cost * self.operating_days*24

        fixed_cost = (
            self.labor_cost * (1 + self.fringe_benefits + self.supplies)
            + (
                self.administration
                + self.maintenance
                + self.property_insurance
                + self.property_tax
            )
            * self.FCI()
        )

        # Discounted Cash Flow Rate of Return Analysis

        sales = (
            ones(self.duration[1] - self.duration[0])
            * sum([i.cost for i in self.system.products])
            * self.operating_days
            * 24
        )
        sales[0] = (1 - self.startup_months / 12) * sales[
            0
        ] + self.startup_months / 12 * self.startup_salesfrac * sales[0]

        total_variable_costs = ones(self.duration[1] - self.duration[0]) * material_cost
        total_variable_costs[0] = (1.0 - self.startup_months / 12) * (
            material_cost
        ) + self.startup_months / 12 * self.startup_VOCfrac * (material_cost)

        total_fixed_costs = ones(self.duration[1] - self.duration[0]) * (fixed_cost)
        total_fixed_costs[0] = (1.0 - self.startup_months / 12) * (
            fixed_cost
        ) + self.startup_months / 12 * self.startup_FOCfrac * (fixed_cost)

        total_product_costs = total_fixed_costs + total_variable_costs

        depreciation_periods = 7
        DDB = zeros(self.duration[1] - self.duration[0])
        SL = zeros(self.duration[1] - self.duration[0])
        Remaining = zeros(self.duration[1] - self.duration[0])
        Actual = zeros(self.duration[1] - self.duration[0])
        Salvage = 0

        for year in range(0, 7):
            if year < 1:
                DDB[year] = (self.FCI() - Salvage) * 200 / (100 * depreciation_periods)
                SL[year] = self.FCI() / (depreciation_periods - year)
                Remaining[year] = self.FCI() - Salvage - DDB[year]
                Actual[year] = SL[year] if SL[year] > DDB[year] else DDB[year]
            else:
                DDB[year] = Remaining[year - 1] * 200 / (100 * depreciation_periods)
                SL[year] = (
                    SL[year - 1]
                    if SL[year - 1] > DDB[year - 1]
                    else Remaining[year - 1] / (depreciation_periods - year)
                )
                Remaining[year] = Remaining[year - 1] - DDB[year]
                Actual[year] = SL[year] if SL[year] > DDB[year] else DDB[year]

        total_annual_sales = sales
        net_revenue = total_annual_sales - total_product_costs - Actual

        losses_forward = zeros(self.duration[1] - self.duration[0])
        taxable_income = zeros(self.duration[1] - self.duration[0])
        income_tax = zeros(self.duration[1] - self.duration[0])

        annual_cash_income = zeros(self.duration[1] - self.duration[0])
        annual_present_value = zeros(self.duration[1] - self.duration[0])
        for year in range(0, len(taxable_income)):
            if year >= 1:
                losses_forward[year] = (
                    taxable_income[year - 1] if taxable_income[year - 1] < 0 else 0
                )
            taxable_income[year] = losses_forward[year] + net_revenue[year]

            income_tax[year] = (
                0
                if taxable_income[year] * self.income_tax < 0
                else taxable_income[year] * self.income_tax
            )
            annual_cash_income[year] = (
                total_annual_sales[year] - total_product_costs[year] - income_tax[year]
            )
            annual_present_value[year] = annual_cash_income[year] / (1 + self.IRR) ** (
                year + 1
            )

        # pioneer plant performance
        for y, v in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
            if self.pioneer + v < 1:
                annual_present_value[y] = annual_present_value[y] * (self.pioneer + v)

        fixed_capital_investment = zeros(3)
        investment_stages = [0.08, 0.6, 0.32]
        total_capital_interest = zeros(3)

        for year in range(-2, 1):
            if year == -2:
                fixed_capital_investment[year + 2] = (
                    self.FCI() * investment_stages[year + 2]
                    + self.working_capital * self.FCI()
                )
            else:
                fixed_capital_investment[year + 2] = (
                    self.FCI() * investment_stages[year + 2]
                )

            total_capital_interest[year + 2] = fixed_capital_investment[year + 2] / (
                1 + self.IRR
            ) ** (year)

        total_capital_interest[2] = (
            fixed_capital_investment[2] + self.working_capital * self.FCI()
        ) / (1 + self.IRR) ** (year)

        NPV = (
            sum(annual_present_value)
            - sum(total_capital_interest)
            + self.working_capital
            * self.FCI()
            / (1 + self.IRR) ** (self.duration[1] - self.duration[0])
        )

        return NPV

    def solve_price(self, product):
        base_price = product.price

        def NPVPrice(price):
            product.price = price
            return self.NPV()

        from scipy.optimize import newton

        price = newton(NPVPrice, base_price)
        product.price = base_price
        return price

    def mfsp_table(self, product=None, solve=True):
        costs = {}
        costs[
            "Biomass Collection"
        ] = self.biomass_collection_cost * self.system.flowsheet.stream[
            "Biomass"
        ].get_total_flow(
            "tonnes/year"
        )
        # costs["Biomass Delivery"] = (
        #     self.biomass_delivery_cost
        #     * self.biomass_delivery_distance
        #     * self.system.flowsheet.stream["Biomass"].get_total_flow("tonnes/year")
        # )

        costs["Bio-Oil Delivery"] = (
            self.bio_oil_delivery_cost
            * self.bio_oil_delivery_distance
            * self.system.flowsheet.stream["Biooil"].get_total_flow("tonnes/year")
        )

        costs["Bio-Oil Injection"] = (
            self.bio_oil_injection_cost
            * self.system.flowsheet.stream["Biooil"].get_total_flow("tonnes/year")
        )

        for f in self.system.feeds:
            if abs(f.cost) > 0.1:
                costs[f.ID] = f.cost * self.operating_days * 24
            else:
                costs["Other"] = (
                    costs.get("Other", 0) + f.cost * self.operating_days * 24
                )
        # costs["Utilities"] = self.system.utility_cost
        costs["Utilities"] = 0 
        for unit in self.system.units:
            if unit.utility_cost != None and unit.utility_cost/self.system.utility_cost > 0.10:
                costs["Utility " + unit.ID] = unit.utility_cost * self.operating_days * 24
            elif unit.utility_cost != None:
                costs["Utilities"] += unit.utility_cost * self.operating_days * 24
            
        costs["O&M"] = self.FCI() * (
            self.property_tax
            + self.property_insurance
            + self.maintenance
            + self.administration
        ) + self.labor_cost * (1 + self.fringe_benefits + self.supplies)

        costs["Bio Oil Credit"] = (
            -self.bio_oil_plug_value
            * self.system.flowsheet.stream["Biooil"].get_total_flow("kg/hr")
            * self.operating_days
            * 24
        )

        revenues = 0
        if product == None:
            for f in self.system.products:
                if abs(f.cost) > 0:
                    costs[f.ID] = -f.cost * self.operating_days * 24
                    revenues += -f.cost * self.operating_days * 24
        else:
            for f in self.system.products:
                if product != None and f.ID != product.ID:
                    if abs(f.cost) > 0:
                        costs[f.ID] = -f.cost * self.operating_days * 24
                        revenues += -f.cost * self.operating_days * 24

        if product != None:
            if solve:
                price = self.solve_price(product)
            else:
                price = product.price
            sales = price * product.get_total_flow("kg/year")
            costs["Income Tax"] = self.income_tax * (
                sales - revenues - sum(list(costs.values()))
            )
            if costs["Income Tax"] < 0:
                costs["Income Tax"] = 0

            costs["ROI"] = sales - sum([v for v in costs.values()])
        else:
            costs["ROI"] = (
                self.FCI()
                * (1 + self.working_capital)
                * (self.IRR * (1 + self.IRR) ** (self.duration[1] - self.duration[0]))
                / ((1 + self.IRR) ** (self.duration[1] - self.duration[0]) - 1)
            )

        return costs
