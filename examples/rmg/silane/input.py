# Data sources
database(
    thermoLibraries = ['SiliconHydrideLibrary', 'primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['Silylene_Elimination'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='SiH4',
    reactive=True,
    structure=adjacencyList("""
	1 Si u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
	2 H u0 p0 c0 {1,S}
	3 H u0 p0 c0 {1,S}
	4 H u0 p0 c0 {1,S}
	5 H u0 p0 c0 {1,S}
	""")
)

species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]")
)

# Reaction systems
simpleReactor(
    temperature=(800,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "SiH4": 0.5,
	"H2": 0.5,
    },
    terminationConversion={
        'SiH4': 0.9,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

options(
    units='si',
    saveRestartPeriod=None,
    drawMolecules=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveConcentrationProfiles=True,
)
