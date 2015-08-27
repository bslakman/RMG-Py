#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy

import rmgpy.quantity

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius, SurfaceArrhenius
from rmgpy.thermo import ThermoData
from rmgpy.thermo.thermodata import SurfaceThermoData
from rmgpy.solver.heterogeneous import HeterogeneousReactor
from rmgpy.solver.base import TerminationTime, TerminationConversion
import rmgpy.constants as constants

################################################################################


class HeterogeneousTest(unittest.TestCase):
    def testSolve(self):
        """
        Test the heterogeneous reactor with two reactions: 
        Silylene elimination: SiH4 <=> SiH2 + H2
        Silane adsorption: SiH4 + X + X => SiH3X + HX
        """
        SiH4 = Species(
            molecule=[Molecule().fromAdjacencyList("""
                1 Si u0 p0 {2,S} {3,S} {4,S} {5,S}
                2 H u0 p0 {1,S}
                3 H u0 p0 {1,S}
                4 H u0 p0 {1,S}
                5 H u0 p0 {1,S}
                """)],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([10.27, 12.3, 14.14, 15.75, 18.33,
                                       20.2, 22.83], "cal/(mol*K)"),
                              H298=(8.2, "kcal/mol"),
                              S298=(48.913  , "cal/(mol*K)")))
        SiH2 = Species(
            molecule=[Molecule().fromAdjacencyList("""
                1 Si u0 p1 {2,S} {3,S}
                2 H u0 p0 {1,S}
                3 H u0 p0 {1,S}
                """)],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([8.35, 8.821, 9.381, 9.954, 10.921,
                                       11.651, 12.734], "cal/(mol*K)"),
                              H298=(62.875, "kcal/mol"),
                              S298=(50.932  , "cal/(mol*K)")))
        H2 = Species(
            molecule=[Molecule().fromAdjacencyList("""
                1 H u0 p0 {2,S}
                2 H u0 p0 {1,S}
                """)],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([6.895,6.975,6.994,7.009,7.081,
                                       7.219, 7.72], "cal/(mol*K)"),
                              H298=(0, "kcal/mol"),
                              S298=(31.233, "cal/(mol*K)")))

        X = Species(
            molecule=[Molecule().fromAdjacencyList("1 X u0 p0")],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([0., 0., 0., 0., 0., 0., 0.], "cal/(mol*K)"),
                              H298=(0.0, "kcal/mol"),
                              S298=(0.0, "cal/(mol*K)")))
        # for simplicity use same values for HX and SiH3X.
        HX = Species(
            molecule=[Molecule().fromAdjacencyList("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57], "cal/(mol*K)"),
                              H298=(-11.26, "kcal/mol"),
                              S298=(0.44, "cal/(mol*K)")))
        SiH3X = Species(
            molecule=[Molecule().fromAdjacencyList("""
                1 Si u0 p0 {2,S} {3,S} {4,S} {5,S}
                2 X u0 p0 {1,S}
                3 H u0 p0 {1,S}
                4 H u0 p0 {1,S}
                5 H u0 p0 {1,S}
                """)],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57], "cal/(mol*K)"),
                              H298=(-11.26, "kcal/mol"),
                              S298=(0.44, "cal/(mol*K)")))

        rxn1 = Reaction(reactants=[SiH4, X, X],
                        products=[SiH3X, HX],
                        reversible = False,
                        surface = True,
                        kinetics=SurfaceArrhenius(A=(1.148E19, 'cm^5/(mol^2*s)'),
                                           n=0.5,
                                           Ea=(3000, 'cal/mol'),
                                           T0=(1.0, 'K')))

        rxn2 = Reaction(reactants=[SiH4],
                        products=[SiH2, H2],
                        kinetics=Arrhenius(A=(3.3E15, '1/s'),
                                           n=-0.5,
                                           Ea=(55.9, 'kcal/mol'),
                                           T0=(1.0, 'K')))
        
        coreSpecies = [SiH4, SiH2, H2, X, SiH3X, HX]
        edgeSpecies = []
        coreReactions = [rxn1, rxn2]
        edgeReactions = []

        T = 1200
        P = 0.5E5
        rxnSystem = HeterogeneousReactor(
            T,
            P,
            initialGasMoleFractions={SiH4: 1.0},
            initialSurfaceCoverages={X: 1.0},
            areaToVolRatio=(364, 'm^-1'),
            surfaceSiteDensity=(2.72e-9, 'mol/cm^2'),
            termination=[]
            )

        rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies,
                                  edgeReactions)

        tlist = numpy.array([10 ** (i / 10.0)
                             for i in range(-100, 0)], numpy.float64)

        # Integrate to get the solution at each time point
        t = []
        y = []
        reactionRates = []
        speciesRates = []
        f = open('heterogeneousReactorOutput.txt', 'w')
        for t1 in tlist:
            rxnSystem.advance(t1)
            t.append(rxnSystem.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxnSystem.y.copy())
            reactionRates.append(rxnSystem.coreReactionRates.copy())
            speciesRates.append(rxnSystem.coreSpeciesRates.copy())
            f.write(str(t1) + "\t")
            f.write(str(rxnSystem.coreReactionRates) + "\t")
            f.write(str(rxnSystem.coreSpeciesRates) + "\t")
            f.write(str(rxnSystem.y) + "\t")
            f.write("\n") 
        f.close()

        # Convert the solution vectors to numpy arrays
        t = numpy.array(t, numpy.float64)
        y = numpy.array(y, numpy.float64)
        reactionRates = numpy.array(reactionRates, numpy.float64)
        speciesRates = numpy.array(speciesRates, numpy.float64)
        
        # Check that we've reached equilibrium
        self.assertAlmostEqual(reactionRates[-1, 0], 0.0, delta=1e-2)


        # Visualize the simulation results
        import pylab
        theta = y[:,(3,4)]/1.97e-3 
        f2 = open('theta.txt', 'w')
        f2.write(str(theta))
        f2.close()
        fig = pylab.figure(figsize=(6, 6))
        pylab.subplot(2, 1, 1)
        pylab.semilogx(t, theta)
        pylab.ylabel('Fractional coverage')
        pylab.legend(['empty site X', 'SiH3X'], loc=4)
        pylab.subplot(2, 1, 2)
        pylab.semilogx(t, y[:,0])
        pylab.legend(['SiH4'], loc=2)
        pylab.xlabel('Time(s)')
        pylab.ylabel('y_SiH4')
        fig.subplots_adjust(left=0.12, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        #pylab.show()
        pylab.savefig('heterogeneousTest.pdf')


        return


