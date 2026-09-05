import festim as F
import numpy as np

n_B = 1.31e29

T0 = 300.0

D0 = 7.78e-7
E_D = 0.49023
p0 = 1.0e13

lambda_B = n_B ** (-1.0 / 3.0)
k0 = D0 / (lambda_B**2 * n_B)
E_k = E_D


def simulate_TDS(
    beta,
    thickness,
    inventory,
    fractions,
    energies,
    Kr0,
    E_r,
    Tf,
    return_trap_contribution=False,
    surface_model="arrhenius",
):
    kr0_festim = Kr0 / 2

    trap_density = inventory * n_B * np.asarray(fractions)
    n_traps = len(fractions)

    mobile = F.Species("D")
    trapped = [F.Species(f"D_trapped_{i + 1}", mobile=False) for i in range(n_traps)]

    model = F.HydrogenTransportProblem()
    model.show_progress_bar = False
    model.species = [mobile, *trapped]

    n_cells = max(2, int(np.ceil(thickness / 2e-9)))
    model.mesh = F.Mesh1D(np.linspace(0.0, thickness, n_cells + 1))

    left = F.SurfaceSubdomain1D(id=1, x=0.0)
    right = F.SurfaceSubdomain1D(id=2, x=thickness)
    film = F.VolumeSubdomain1D(
        id=3,
        borders=[0.0, thickness],
        material=F.Material(D_0=D0, E_D=E_D),
    )
    model.subdomains = [left, right, film]

    empty = [
        F.ImplicitSpecies(n=trap_density[i], others=[trapped[i]])
        for i in range(n_traps)
    ]

    model.reactions = [
        F.Reaction(
            reactant=[mobile, empty[i]],
            product=[trapped[i]],
            k_0=k0,
            E_k=E_k,
            p_0=p0,
            E_p=energies[i],
            volume=film,
        )
        for i in range(n_traps)
    ]

    model.initial_conditions = [
        F.InitialConcentration(
            species=trapped[i],
            value=trap_density[i],
            volume=film,
        )
        for i in range(n_traps)
    ]

    if surface_model == "sink":
        left_bc = F.FixedConcentrationBC(
            subdomain=left,
            value=0.0,
            species=mobile,
        )
    else:
        left_bc = F.SurfaceReactionBC(
            reactant=[mobile, mobile],
            gas_pressure=0.0,
            k_r0=kr0_festim,
            E_kr=E_r,
            k_d0=0.0,
            E_kd=0.0,
            subdomain=left,
        )

    model.boundary_conditions = [
        left_bc,
        F.ParticleFluxBC(
            subdomain=right,
            value=0.0,
            species=mobile,
        ),
    ]

    model.temperature = lambda t: T0 + beta * t

    model.settings = F.Settings(
        atol=1e10,
        rtol=1e-10,
        final_time=(Tf - T0) / beta,
        max_iterations=50,
    )

    model.settings.stepsize = F.Stepsize(
        initial_value=0.01,
        growth_factor=1.2,
        cutback_factor=0.9,
        target_nb_iterations=4,
        max_stepsize=1.0 / beta,
    )

    surface_flux = F.SurfaceFlux(field=mobile, surface=left)

    if return_trap_contribution:
        contents = [F.TotalVolume(field=species, volume=film) for species in trapped]
        model.exports = [surface_flux, *contents]
    else:
        contents = []
        model.exports = [surface_flux]

    model.initialise()
    model.run()

    time = np.asarray(surface_flux.t)
    temperature = T0 + beta * time
    flux = np.asarray(surface_flux.data)

    if not return_trap_contribution:
        return temperature, flux

    trap_contribution = [
        -np.diff(np.asarray(content.data)) / np.diff(time) for content in contents
    ]

    return temperature, flux, trap_contribution
