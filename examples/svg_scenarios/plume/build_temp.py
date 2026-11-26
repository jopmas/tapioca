'''
Mantle plume at the centre of the domain, for custom_temp_field of plume.svg.

build_temperature receives the scenario with its initial temperature already built (steady-state conduction in the
lithosphere and the mantle adiabat below the LAB) and adds a plume with TemperatureFieldBuilder.add_plume: a head
below the lithosphere and a tail (conduit) that rises from the base of the domain, both DT_PLUME hotter than the
mantle around them.

Other anomalies (see TemperatureFieldBuilder):
    t.add_basal_anomaly(200.0, thickness=50e3)                         # hot thermal boundary layer at the base
    t.add_basal_anomaly(200.0, thickness=50e3, x0=800e3, width=200e3)  # only below the centre
    t.add_thermal_anomaly(-150.0, x0=400e3, z0=-300e3, rx=100e3, rz=50e3)  # cold blob (e.g. a detached slab)
z is positive upwards: -Lz at the base of the domain and 0 at its top (the top of the sticky air).
'''

DT_PLUME = 250.0      # temperature excess of the plume [oC]
HEAD_DEPTH = 300.0e3  # depth of the centre of the head below the surface [m]
HEAD_RADIUS = 80.0e3  # radius of the head [m]
TAIL_RADIUS = 25.0e3  # half width of the tail [m]


def _surface(scenario):
    '''z of the surface [m] (base of the sticky air).'''
    return -float(scenario.thickness_sa.max())


def build_temperature(scenario):
    '''Plume at x = Lx/2 added to the initial temperature of the scenario.'''
    t = scenario.temperature_builder
    t.add_plume(DT_PLUME, x0=scenario.Lx / 2, z_head=_surface(scenario) - HEAD_DEPTH,
                head_radius=HEAD_RADIUS, tail_radius=TAIL_RADIUS)
    return t


def infos(scenario):
    '''Lines of infos_<folder>.txt about the plume.'''
    return [f"plume at x = {scenario.Lx / 2e3:.0f} km: +{DT_PLUME:g} oC, head of radius {HEAD_RADIUS / 1e3:g} km "
            f"centred {HEAD_DEPTH / 1e3:g} km below the surface, tail of half width {TAIL_RADIUS / 1e3:g} km "
            f"from the base of the domain"]
