# Data sources
database(
    thermoLibraries = ['SiliconHydrideLibrary', 'primaryThermoLibrary'],
    reactionLibraries = [('Silicon_Giunta_1990', False), ('DolletSi2H4', False)],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['Silylene_Insertion', 'Silylene_to_Silene', 'H_Abstraction', 'H2_transfer'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='SiH4',
    reactive=True,
    structure=SMILES("[SiH4]")
)

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]")
)

species(
    label='SiH3',
    reactive=True,
    structure=SMILES("[SiH3]")
)

species(
    label='H',
    reactive=True,
    structure=SMILES("[H]")
)

species(
    label='site',
    reactive=True,
    structure=adjacencyList("1 X u0")
)

# Reaction systems
heterogeneousReactor(
    temperature=(1000,'K'),
    pressure=(0.5,'bar'),
    initialGasMoleFractions={
        "SiH4": 0.098,
        "SiH3": 0.001,
        "H": 0.001,
	"Ar": 0.9
    },
    initialSurfaceCoverages={
        "site": 1.0
    },
    terminationConversion={
        'SiH4': 0.9,
    },
    terminationTime=(1, 'hour'),
    areaToVolRatio=(1.0, 'm^-1'),
    surfaceSiteDensity=(2.72e-9, 'mol/cm^2')
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.001,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2000,'K',8),
    pressures=(0.01,20,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    saveRestartPeriod=None,
    drawMolecules=False,
    generatePlots=False,
    saveEdgeSpecies=True,
)


