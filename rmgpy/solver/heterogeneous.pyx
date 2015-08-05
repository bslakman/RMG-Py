################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains the :class:`HeterogeneousReactor` class, providing a reaction system
consisting of a isothermal, isobaric batch reactor which includes both gas and surface reactions.
"""

import numpy
cimport numpy

include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
else:
    from pydas.dassl cimport DASSL as DASx
    
from base cimport ReactionSystem
cimport cython

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

cdef class HeterogeneousReactor(ReactionSystem):
    """
    A reaction system consisting of an isothermal, isobaric batch
    reactor.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P
    cdef public double V
    cdef public dict initialGasMoleFractions
    cdef public list sensitiveSpecies
    cdef public double sensitivityThreshold
    # surface parameters
    cdef public dict initialSurfaceCoverages
    cdef public ScalarQuantity surfaceSiteDensity 
    cdef public ScalarQuantity areaToVolRatio
    cdef public double A

    cdef public numpy.ndarray reactantIndices
    cdef public numpy.ndarray productIndices
    cdef public numpy.ndarray networkIndices
    cdef public numpy.ndarray forwardRateCoefficients
    cdef public numpy.ndarray reverseRateCoefficients
    cdef public numpy.ndarray equilibriumConstants
    cdef public numpy.ndarray networkLeakCoefficients
    cdef public numpy.ndarray jacobianMatrix
    cdef public numpy.ndarray surfaceReactions

    def __init__(self, T, P, initialGasMoleFractions, initialSurfaceCoverages, areaToVolRatio, surfaceSiteDensity, termination, sensitiveSpecies=None, sensitivityThreshold=1e-3):
        ReactionSystem.__init__(self, termination)
        self.T = Quantity(T)
        self.P = Quantity(P)
        self.initialGasMoleFractions = initialGasMoleFractions
        self.initialSurfaceCoverages = initialSurfaceCoverages
        self.surfaceSiteDensity = surfaceSiteDensity

        self.V = 0 # will be set in initializeModel
        self.constantVolume = False
        self.sensitiveSpecies = sensitiveSpecies
        self.sensitivityThreshold = sensitivityThreshold
        self.areaToVolRatio = areaToVolRatio
        self.A = 0

        # These are helper variables used within the solver
        self.reactantIndices = None
        self.productIndices = None
        self.networkIndices = None
        self.forwardRateCoefficients = None
        self.reverseRateCoefficients = None
        self.equilibriumConstants = None
        self.networkLeakCoefficients = None
        self.jacobianMatrix = None

        self.surfaceReactions = None
        
    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialMoleFractions dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initialGasMoleFractions = {}
        for label, moleFrac in self.initialGasMoleFractions.iteritems():
            initialGasMoleFractions[speciesDict[label]] = moleFrac
        self.initialGasMoleFractions = initialGasMoleFractions

        initialSurfaceCoverages = {}
        for label, surfaceCoverage in self.initialSurfaceCoverages.iteritems():
            initialSurfaceCoverages[speciesDict[label]] = surfaceCoverage
        self.initialSurfaceCoverages = initialSurfaceCoverages
    
    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, atol, rtol, sensitivity, sens_atol, sens_rtol)

        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, l, index, neq
        cdef double V, areaToVolRatio
        cdef dict speciesIndex, reactionIndex
        cdef numpy.ndarray[numpy.int_t, ndim=2] reactantIndices, productIndices, networkIndices
        cdef numpy.ndarray[numpy.float64_t, ndim=1] forwardRateCoefficients, reverseRateCoefficients, equilibriumConstants, networkLeakCoefficients, atol_array, rtol_array, senpar
        
        pdepNetworks = pdepNetworks or []

        numCoreSpecies = len(coreSpecies)
        numCoreReactions = len(coreReactions)
        numEdgeSpecies = len(edgeSpecies)
        numEdgeReactions = len(edgeReactions)
        numPdepNetworks = len(pdepNetworks)

        surfaceSpecies = numpy.zeros_like(numCoreSpecies, dtype=bool)
        
        # Assign an index to each species (core first, then edge)
        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
            if spec.isSurfaceSpecies():
               surfaceSpecies[index] = True 
        for index, spec in enumerate(edgeSpecies):
            speciesIndex[spec] = index + numCoreSpecies
        # Assign an index to each reaction (core first, then edge)
        reactionIndex = {}
        for index, rxn in enumerate(coreReactions):
            reactionIndex[rxn] = index
        for index, rxn in enumerate(edgeReactions):
            reactionIndex[rxn] = index + numCoreReactions

        # Generate reactant and product indices
        # Generate forward and reverse rate coefficients k(T,P)
        reactantIndices = -numpy.ones((numCoreReactions + numEdgeReactions, 3), numpy.int )
        productIndices = -numpy.ones_like(reactantIndices)
        forwardRateCoefficients = numpy.zeros((numCoreReactions + numEdgeReactions), numpy.float64)
        reverseRateCoefficients = numpy.zeros_like(forwardRateCoefficients)
        equilibriumConstants = numpy.zeros_like(forwardRateCoefficients)
        surfaceReactions = numpy.zeros_like(forwardRateCoefficients, dtype=bool) 
        for rxnList in [coreReactions, edgeReactions]:
            for rxn in rxnList:
                j = reactionIndex[rxn]
                forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T.value_si, self.P.value_si)
                if rxn.isSurfaceReaction:
                    surfaceReactions[j] = True
                if rxn.reversible:
                    equilibriumConstants[j] = rxn.getEquilibriumConstant(self.T.value_si)
                    reverseRateCoefficients[j] = forwardRateCoefficients[j] / equilibriumConstants[j]
                for l, spec in enumerate(rxn.reactants):
                    i = speciesIndex[spec]
                    reactantIndices[j,l] = i
                for l, spec in enumerate(rxn.products):
                    i = speciesIndex[spec]
                    productIndices[j,l] = i

        networkIndices = -numpy.ones((numPdepNetworks, 3), numpy.int )
        networkLeakCoefficients = numpy.zeros((numPdepNetworks), numpy.float64)
        for j, network in enumerate(pdepNetworks):
            networkLeakCoefficients[j] = network.getLeakCoefficient(self.T.value_si, self.P.value_si)
            for l, spec in enumerate(network.source):
                i = speciesIndex[spec]
                networkIndices[j,l] = i

        self.reactantIndices = reactantIndices
        self.productIndices = productIndices
        self.forwardRateCoefficients = forwardRateCoefficients
        self.reverseRateCoefficients = reverseRateCoefficients
        self.equilibriumConstants = equilibriumConstants
        self.networkIndices = networkIndices
        self.networkLeakCoefficients = networkLeakCoefficients
        self.surfaceReactions = surfaceReactions
        self.surfaceReactions = surfaceSpecies

        # Set initial conditions
        t0 = 0.0
        # Compute number of equations    
        if sensitivity:    
            # Set DASPK sensitivity analysis to ON
            self.sensitivity = True
            # Compute number of variables
            neq = numCoreSpecies*(len(forwardRateCoefficients)+numCoreSpecies+1)
            
            atol_array = numpy.ones(neq, numpy.float64)*sens_atol
            atol_array[:numCoreSpecies] = atol
            
            rtol_array = numpy.ones(neq, numpy.float64)*sens_rtol
            rtol_array[:numCoreSpecies] = rtol
            
            senpar = numpy.zeros(len(forwardRateCoefficients)+numCoreSpecies, numpy.float64)
            
        else:
            neq = numCoreSpecies
            
            atol_array = numpy.ones(neq,numpy.float64)*atol
            rtol_array = numpy.ones(neq,numpy.float64)*rtol
            
            senpar = numpy.zeros(len(forwardRateCoefficients), numpy.float64)
            
        y0 = numpy.zeros(neq, numpy.float64)
        for spec, moleFrac in self.initialGasMoleFractions.iteritems():
            y0[speciesIndex[spec]] = moleFrac
            
        # Use ideal gas law to compute volume
        V = constants.R * self.T.value_si * numpy.sum(y0[:numCoreSpecies]) / self.P.value_si
        self.V = V # volume in m^3
        for j in range(numCoreSpecies):
            self.coreSpeciesConcentrations[j] = y0[j] / V
        
        # Compute area with user-specified ratio
        self.A = self.V * self.areaToVolRatio
        totalSurfaceSites = V * areaToVolRatio * self.surfaceSiteDensity # total surface sites in reactor
        
        for spec, coverage in self.initialSurfaceCoverages.iteritems():
            y0[speciesIndex[spec]] = totalSurfaceSites * coverage

        for j, isSurfaceSpecies in enumerate(self.surfaceSpecies): # should only go up to core species
            if isSurfaceSpecies:
               self.coreSpeciesConcentrations[j] = y0[j] / V / areaToVolRatio # moles per m2 of surface
            else:
                self.coreSpeciesConcentrations[j] = y0[j] / V # moles per m3 of gas       

        # Initialize the model
        dydt0 = - self.residual(t0, y0, numpy.zeros(neq, numpy.float64), senpar)[0]
        DASx.initialize(self, t0, y0, dydt0, senpar, atol_array, rtol_array)

    @cython.boundscheck(False)
    def residual(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt, numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)):

        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip, inet
        cdef numpy.ndarray surfaceReactions, surfaceSpecies
        cdef numpy.ndarray[numpy.float64_t, ndim=1] res, kf, kr, knet, delta, equilibriumConstants
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, z, first, second, third
        cdef double k, V, reactionRate, areaToVolRatio
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibriumConstants = self.equilibriumConstants

        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        
        inet = self.networkIndices
        knet = self.networkLeakCoefficients

        surfaceReactions = self.surfaceReactions
        surfaceSpecies = self.surfaceSpecies

        numCoreSpecies = len(self.coreSpeciesRates)
        numCoreReactions = len(self.coreReactionRates)
        numEdgeSpecies = len(self.edgeSpeciesRates)
        numEdgeReactions = len(self.edgeReactionRates)
        numPdepNetworks = len(self.networkLeakRates)
        
        
        res = numpy.zeros(numCoreSpecies, numpy.float64)

        coreSpeciesConcentrations = numpy.zeros_like(self.coreSpeciesConcentrations)
        coreSpeciesRates = numpy.zeros_like(self.coreSpeciesRates)
        coreReactionRates = numpy.zeros_like(self.coreReactionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        
        # Use ideal gas law to compute volume
        V = constants.R * self.T.value_si * numpy.sum(y[:numCoreSpecies]) / self.P.value_si
        self.V = V

        for j in range(numCoreSpecies):
            if surfaceSpecies[j]:
                C[j] = (y[j] / V) / areaToVolRatio 
            else:
                C[j] = y[j] / V
            coreSpeciesConcentrations[j] = C[j]
        


        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j,0] >= numCoreSpecies or ir[j,1] >= numCoreSpecies or ir[j,2] >= numCoreSpecies:
                reactionRate = 0.0
            elif ir[j,1] == -1: # only one reactant
                reactionRate = k * C[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                reactionRate = k * C[ir[j,0]] * C[ir[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]]
            k = kr[j]
            if ip[j,0] >= numCoreSpecies or ip[j,1] >= numCoreSpecies or ip[j,2] >= numCoreSpecies:
                pass
            elif ip[j,1] == -1: # only one reactant
                reactionRate -= k * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                reactionRate -= k * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                reactionRate -= k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]

            if self.surfaceReactions[j]: 
               reactionRate = reactionRate * self.vacantSiteFraction # k should include sticking coeff?
            
            # Set the reaction and species rates
            if j < numCoreReactions:
                # The reaction is a core reaction
                coreReactionRates[j] = reactionRate

                # Need to convert the gas and surface reaction rates to molar
                # flow rate by multiplying by area or volume accordingly
                if self.surfaceReactions[j]: 
                   molarFlowRate = reactionRate * self.A 
                else:
                   molarFlowRate = reactionRate * self.V

                # Add/substract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j,0]
                coreSpeciesRates[first] -= molarFlowRate
                second = ir[j,1]
                if second != -1:
                    coreSpeciesRates[second] -= molarFlowRate
                    third = ir[j,2]
                    if third != -1:
                        coreSpeciesRates[third] -= molarFlowRate
                first = ip[j,0]
                coreSpeciesRates[first] += molarFlowRate 
                second = ip[j,1]
                if second != -1:
                    coreSpeciesRates[second] += molarFlowRate
                    third = ip[j,2]
                    if third != -1:
                        coreSpeciesRates[third] += molarFlowRate 

            else:
                # The reaction is an edge reaction
                edgeReactionRates[j-numCoreReactions] = reactionRate

                # Need to convert the gas and surface reaction rates to molar
                # flow rate by multiplying by area or volume accordingly
                if self.surfaceReactions[j]: 
                   molarFlowRate = reactionRate * self.A 
                else:
                   molarFlowRate = reactionRate * self.V
                
                # Add/substract the total reaction rate from each species rate
                # Since it's an edge reaction its reactants and products could
                # be either core or edge species
                # We're only interested in the edge species
                first = ir[j,0]
                if first >= numCoreSpecies: edgeSpeciesRates[first-numCoreSpecies] -= molarFlowRate 
                second = ir[j,1]
                if second != -1:
                    if second >= numCoreSpecies: edgeSpeciesRates[second-numCoreSpecies] -=  molarFlowRate 
                    third = ir[j,2]
                    if third != -1:
                        if third >= numCoreSpecies: edgeSpeciesRates[third-numCoreSpecies] -= molarFlowRate 
                first = ip[j,0]
                if first >= numCoreSpecies: edgeSpeciesRates[first-numCoreSpecies] += molarFlowRate
                second = ip[j,1]
                if second != -1:
                    if second >= numCoreSpecies: edgeSpeciesRates[second-numCoreSpecies] += molarFlowRate 
                    third = ip[j,2]
                    if third != -1:
                        if third >= numCoreSpecies: edgeSpeciesRates[third-numCoreSpecies] += molarFlowRate 

        for j in range(inet.shape[0]):
            k = knet[j]
            if inet[j,1] == -1: # only one reactant
                reactionRate = k * C[inet[j,0]]
            elif inet[j,2] == -1: # only two reactants
                reactionRate = k * C[inet[j,0]] * C[inet[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * C[inet[j,0]] * C[inet[j,1]] * C[inet[j,2]]
            networkLeakRates[j] = reactionRate

        self.coreSpeciesConcentrations = coreSpeciesConcentrations
        self.coreSpeciesRates = coreSpeciesRates
        self.coreReactionRates = coreReactionRates
        self.edgeSpeciesRates = edgeSpeciesRates
        self.edgeReactionRates = edgeReactionRates
        self.networkLeakRates = networkLeakRates

        res = coreSpeciesRates # Already in units of moles/time
        
        if self.sensitivity:
            delta = numpy.zeros(len(y), numpy.float64)
            delta[:numCoreSpecies] = res
            if self.jacobianMatrix is None:
                jacobian = self.jacobian(t,y,dydt,0,senpar)
            else:
                jacobian = self.jacobianMatrix
            dgdk = self.computeRateDerivative()
            for j in range(numCoreReactions+numCoreSpecies):
                for i in range(numCoreSpecies):
                    for z in range(numCoreSpecies):
                        delta[(j+1)*numCoreSpecies + i] += jacobian[i,z]*y[(j+1)*numCoreSpecies + z] 
                    delta[(j+1)*numCoreSpecies + i] += dgdk[i,j]

        else:
            delta = res
        delta = delta - dydt
        
        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1
    
    @cython.boundscheck(False)
    def jacobian(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt, double cj, numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)):
        """
        Return the analytical Jacobian for the reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip
        cdef numpy.ndarray[numpy.float64_t, ndim=1] kf, kr, C, equilibriumConstants
        cdef numpy.ndarray[numpy.float64_t, ndim=2] pd
        cdef int numCoreReactions, numCoreSpecies, i, j
        cdef double k, V, Ctot, deriv, corr
        
        ir = self.reactantIndices
        ip = self.productIndices
        equilibriumConstants = self.equilibriumConstants

        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)
        
        pd = -cj * numpy.identity(numCoreSpecies, numpy.float64)
        
        V = constants.R * self.T.value_si * numpy.sum(y[:numCoreSpecies]) / self.P.value_si
        
        Ctot = self.P.value_si /(constants.R * self.T.value_si)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        for j in range(numCoreSpecies):
            C[j] = y[j] / V

        for j in range(numCoreReactions):
           
            k = kf[j]
            if ir[j,1] == -1: # only one reactant
                deriv = k
                pd[ir[j,0], ir[j,0]] -= deriv
                
                pd[ip[j,0], ir[j,0]] += deriv                
                if ip[j,1] != -1:
                    pd[ip[j,1], ir[j,0]] += deriv
                    if ip[j,2] != -1:
                        pd[ip[j,2], ir[j,0]] += deriv
                
                                
            elif ir[j,2] == -1: # only two reactants
                corr = - k * C[ir[j,0]] * C[ir[j,1]] / Ctot
                if ir[j,0] == ir[j,1]:  # reactants are the same
                    deriv = 2 * k * C[ir[j,0]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr    
                    
                else:
                    # Derivative with respect to reactant 1
                    deriv = k * C[ir[j, 1]]
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= deriv                        
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    
                    # Derivative with respect to reactant 2
                    deriv = k * C[ir[j, 0]] 
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr     
                            
                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr               
                    
                    
            else: # three reactants!! (really?)
                corr = - 2* k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]] / Ctot
                if (ir[j,0] == ir[j,1] & ir[j,0] == ir[j,2]):
                    deriv = 3 * k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,0]] -= 3 * deriv                                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 3 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr        
                    
                elif ir[j,0] == ir[j,1]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv                  
                    pd[ir[j,2], ir[j,0]] -= deriv    
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    
                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,2]] -= 2 * deriv                  
                    pd[ir[j,2], ir[j,2]] -= deriv                                                                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr    
                    
                    
                elif ir[j,1] == ir[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j,1]] * C[ir[j,1]] 
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= 2 * deriv
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv  
                    # derivative with respect to reactant 2
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,1]]
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= 2 * deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= 2 * corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                
                elif ir[j,0] == ir[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,1]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv                  
                    pd[ir[j,1], ir[j,0]] -= deriv    
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,1]] -= 2 * deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,1], i] -= corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                                
                else:
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j,1]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= deriv
                    pd[ir[j,2], ir[j,0]] -= deriv
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv     
                                    
                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j,0]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv   
                    pd[ir[j,2], ir[j,1]] -= deriv
                    
                    pd[ip[j,0], ir[j,1]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv 
                                 
                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j,0]] * C[ir[j,1]]             
                    pd[ir[j,0], ir[j,2]] -= deriv                    
                    pd[ir[j,1], ir[j,2]] -= deriv   
                    pd[ir[j,2], ir[j,2]] -= deriv
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                    
            
            
            k = kr[j]         
            if ip[j,1] == -1: # only one reactant
                deriv = k
                pd[ip[j,0], ip[j,0]] -= deriv
                
                pd[ir[j,0], ip[j,0]] += deriv                
                if ir[j,1] != -1:
                    pd[ir[j,1], ip[j,0]] += deriv
                    if ir[j,2] != -1:
                        pd[ir[j,2], ip[j,0]] += deriv
                
                                
            elif ip[j,2] == -1: # only two reactants
                corr = -k * C[ip[j,0]] * C[ip[j,1]] / Ctot
                if ip[j,0] == ip[j,1]:
                    deriv = 2 * k * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                 
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr   
                    
                else:
                    # Derivative with respect to reactant 1
                    deriv = k * C[ip[j, 1]]
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    
                    # Derivative with respect to reactant 2
                    deriv = k * C[ip[j, 0]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv              
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                      
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in range(numCoreSpecies):
                         pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr              
                    
                    
            else: # three reactants!! (really?)
                corr = - 2 * k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]] / Ctot
                if (ip[j,0] == ip[j,1] & ip[j,0] == ip[j,2]):
                    deriv = 3 * k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,0]] -= 3 * deriv          
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 3 * corr
                    
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr       
                    
                elif ip[j,0] == ip[j,1]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                    
                    pd[ip[j,2], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to reactant 3
                    deriv = k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,2]] -= 2 * deriv                    
                    pd[ip[j,2], ip[j,2]] -= deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,2], i] -= corr

                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr     
                    
                    
                elif ip[j,1] == ip[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = k * C[ip[j,1]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= 2 * deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv                 
                    
                    # derivative with respect to reactant 2
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= 2 * deriv   
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr                    
                                
                elif ip[j,0] == ip[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                  
                    pd[ip[j,1], ip[j,0]] -= deriv    
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,1]] -= 2 * deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,1], i] -= corr

                    pd[ir[j,0], ip[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr    
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr    
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr     
                                
                else:
                    # derivative with respect to reactant 1
                    deriv = k * C[ip[j,1]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= deriv
                    pd[ip[j,2], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv     
                                    
                    # derivative with respect to reactant 2
                    deriv = k * C[ip[j,0]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv   
                    pd[ip[j,2], ip[j,1]] -= deriv
                    
                    pd[ir[j,0], ip[j,1]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv 
                                 
                    # derivative with respect to reactant 3
                    deriv = k * C[ip[j,0]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,2]] -= deriv                    
                    pd[ip[j,1], ip[j,2]] -= deriv   
                    pd[ip[j,2], ip[j,2]] -= deriv 
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                        pd[ip[j,2], i] -= corr
                    
                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr  

        self.jacobianMatrix = pd + cj * numpy.identity(numCoreSpecies, numpy.float64)
        return pd
    
    @cython.boundscheck(False)
    def computeRateDerivative(self):
        """
        Returns derivative vector df/dk_j where dy/dt = f(y, t, k) and
        k_j is the rate parameter for the jth core reaction.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip
        cdef numpy.ndarray[numpy.float64_t, ndim=1] kf, kr, C, deriv
        cdef numpy.ndarray[numpy.float64_t, ndim=2] rateDeriv
        cdef double fderiv, rderiv, flux, V
        cdef int j, numCoreReactions, numCoreSpecies
        
        cdef double RT_inverse, gderiv
        
        ir = self.reactantIndices
        ip = self.productIndices
        
        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients    
        
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)      
        
        # Use stored volume, since this function is only called from residual function. 
        RT_inverse = 1/(constants.R * self.T.value_si)
        V = self.V

        C = self.coreSpeciesConcentrations

        rateDeriv = numpy.zeros((numCoreSpecies,numCoreReactions+numCoreSpecies), numpy.float64)
        
        for j in range(numCoreReactions):
            if ir[j,1] == -1: # only one reactant
                fderiv = C[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                fderiv = C[ir[j,0]] * C[ir[j,1]]                             
            else: # three reactants!! (really?)
                fderiv = C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]]          
                
            if ip[j,1] == -1: # only one reactant
                rderiv = kr[j] / kf[j] * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                rderiv = kr[j] / kf[j] * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                rderiv = kr[j] / kf[j] * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]
            
            flux = fderiv - rderiv
            gderiv = rderiv * kf[j] * RT_inverse
            
            deriv = numpy.zeros(numCoreSpecies, numpy.float64) # derivative for reaction j with respect to dG_species i

            deriv[ir[j,0]] += gderiv
            if ir[j,1] != -1: # only two reactants
                deriv[ir[j,1]] += gderiv
                if ir[j,2] != -1: # three reactants!! (really?)
                    deriv[ir[j,2]] += gderiv
            
            deriv[ip[j,0]] -= gderiv
            if ip[j,1] != -1: # only two reactants
                deriv[ip[j,1]] -= gderiv
                if ip[j,2] != -1: # three reactants!! (really?)
                    deriv[ip[j,2]] -= gderiv
            
            rateDeriv[ir[j,0], j] -= flux
            rateDeriv[ir[j,0], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
            if ir[j,1] != -1:
                rateDeriv[ir[j,1], j] -= flux
                rateDeriv[ir[j,1], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
                if ir[j,2] != -1:
                    rateDeriv[ir[j,2], j] -= flux
                    rateDeriv[ir[j,2], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
                
            rateDeriv[ip[j,0], j] += flux
            rateDeriv[ip[j,0], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
            if ip[j,1] != -1:
                rateDeriv[ip[j,1], j] += flux
                rateDeriv[ip[j,1], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
                if ip[j,2] != -1:
                    rateDeriv[ip[j,2], j] += flux  
                    rateDeriv[ip[j,2], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
                        
        rateDeriv = V * rateDeriv

        return rateDeriv
