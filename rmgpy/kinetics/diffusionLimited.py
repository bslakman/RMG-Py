
import logging
import math
#from rmgpy.species import Species
#from rmgpy.data.solvation import SolventData, SoluteData, SoluteGroups, SolvationDatabase


class LiquidKinetics():
    """
    Object for correcting a reaction's intrinsic rate, based on its solvent, reaction 
    family and intrinsic rate in the gas phase. Created and called by the 
    diffusionLimiter, since if solvent is present, both affect the rate.
    """
    
    def __init__(self, reaction, solventData):
        self.reaction = reaction
        self.solventData = solventData
    
    def correctIntrinsicRate(self, T):
        
        gasKinetics = self.reaction.kinetics
        k = gasKinetics.getRateCoefficient(T, P=1E5)
       
        solvationKineticsDatabase = self.reaction.family.solvationCorrections 
        x = solvationKineticsDatabase.estimateBarrierCorrection()
        return k*math.exp(-x/8.314/T)
      
        # """
        # If it is an H-abstraction reaction, correct the intrinsic rate here (best place?)
        # """
        # if self.reaction.family.label == 'H_Abstraction':
        #     # get log10 of ratio of intrinsic constants k_solv/k_gas
        #     correction = self.solventData.getHAbsCorrection()
        #     return (10**correction)*k
        # else:
        #     return k
        # if self.reaction.family.label == 'R_Addition_MultipleBond':
        #     # get x, the correction to the barrier
        #     x = 0
        #     # Assume Arrhenius, so correct k accordingly
        #     return k*math.exp(-x/8.314/T)
        # if self.reaction.family.label == 'Beta_Scission':
        #     # rate increases with Dimroth-Reichardt
        #     return k
        # return k

class DiffusionLimited():

    def __init__(self):
    # default is false, enabled if there is a solvent
        self.enabled = False

    def enable(self, solventData, solvationDatabase, comment=''):
    # diffusionLimiter is enabled if a solvent has been added to the RMG object.
        logging.info("Enabling diffusion-limited kinetics...")
        diffusionLimiter.enabled = True
        diffusionLimiter.database = solvationDatabase
        diffusionLimiter.solventData = solventData

    def getSolventViscosity(self, T):
        return self.solventData.getSolventViscosity(T)
            
                
    def getEffectiveRate(self, reaction, T):
        """
        Return the ratio of k_eff to k_intrinsic, which is between 0 and 1.
        
        It is 1.0 if diffusion has no effect.
        
        For 1<=>2 reactions, the reverse rate is limited.
        For 2<=>2 reactions, the faster direction is limited.
        For 2<=>1 or 2<=>3 reactions, the forward rate is limited.
        """
        liquidKinetics = LiquidKinetics(reaction, self.solventData)
        
        reactants = len(reaction.reactants)
        products = len(reaction.products)
        
        k_forward = liquidKinetics.correctIntrinsicRate(T)      
        Keq = reaction.getEquilibriumConstant(T) # Kc
        k_reverse = k_forward / Keq
                          
        if reactants == 1:
            if products == 1:
                k_eff = k_forward
            else: # two products; reverse rate is limited
                k_diff = self.getDiffusionLimit(T, reaction, forward=False)
                k_eff_reverse = k_reverse*k_diff/(k_reverse+k_diff)
                k_eff = k_eff_reverse * Keq
        else: # 2 reactants
            if products == 1 or products == 3:
                k_diff = self.getDiffusionLimit(T, reaction, forward=True)
                k_eff = k_forward*k_diff/(k_forward+k_diff)
            else: # 2 products
                if Keq > 1.0: # forward rate is faster and thus limited
                    k_diff = self.getDiffusionLimit(T, reaction, forward=True)
                    k_eff = k_forward*k_diff/(k_forward+k_diff)
                else: # reverse rate is faster and thus limited
                    k_diff = self.getDiffusionLimit(T, reaction, forward=False)
                    k_eff_reverse = k_reverse*k_diff/(k_reverse+k_diff)
                    k_eff = k_eff_reverse * Keq
        return k_eff        
    
    def getDiffusionLimit(self, T, reaction, forward=True):
        """
        Return the diffusive limit on the rate coefficient, k_diff.
        
        This is the upper limit on the rate, in the specified direction.
        (ie. forward direction if forward=True [default] or reverse if forward=False)
        """
        if forward:
            reacting = reaction.reactants
        else:
            reacting = reaction.products
        assert len(reacting)==2, "Can only calculate diffusion limit in a bimolecular direction"
        radii = 0.0
        diffusivities = 0.0
        for spec in reacting:
            try:
                soluteData = spec.soluteData
            except AttributeError:
                logging.warning("Species didn't have soluteData.")
                soluteData = self.database.getSoluteData(spec)
            # calculate radius with the McGowan volume and assuming sphere
            radius = ((75*soluteData.V/3.14159)**(1/3))/100
            diff = soluteData.getStokesDiffusivity(T, self.getSolventViscosity(T))
            radii += radius
            diffusivities += diff
        N_a = 6.022e23 # Avogadro's Number
        k_diff = 4*3.14159*radii*diffusivities*N_a
        return k_diff


# module level variable. There should only ever be one. It starts off disabled
diffusionLimiter = DiffusionLimited()
