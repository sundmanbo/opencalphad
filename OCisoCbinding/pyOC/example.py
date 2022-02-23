#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os


# ### Imports for pyOC

# In[2]:


import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat


# ### Setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way

# In[3]:


oc.setVerbosity(True)


# ### reading database (.tdb file)

# In[4]:


tdbFile=os.environ.get('OCPUBLICDATA')+'/steel7.TDB'
oc.readtdb(tdbFile)


# ### Play with phase status

# In[5]:


oc.setPhasesStatus(('* ',),phStat.Suspended)
phaseNames=('FCC_A1','M23C6','M6C')
oc.setPhasesStatus(phaseNames,phStat.Entered)


# ### Set pressure and temperature

# In[6]:


oc.setPressure(1E5)
oc.setTemperature(1173)


# ### Set element molar amounts

# In[7]:


elementMolarAmounts = {
	'C' : 0.04,
	'CR' : 0.06,
	'MO': 0.05,
	'SI': 0.003,
	'V': 0.01,
	'FE': 1.0-0.04-0.06-0.05-0.003-0.01
}
oc.setElementMolarAmounts(elementMolarAmounts)


# ### Calculate equilibrium without the grid-minimizer (equilibrium record is 'eq2')

# In[8]:


oc.changeEquilibriumRecord('eq2')
oc.calculateEquilibrium(gmStat.Off)


# ### Calculate equilibrium with the grid-minimizer (equilibrium record is 'default equilibrium')

# In[9]:


oc.changeEquilibriumRecord()
oc.calculateEquilibrium()


# ### Retrieving Gibbs energies and comparing them

# In[10]:


oc.changeEquilibriumRecord()
G=oc.getGibbsEnergy() # a scalar
oc.changeEquilibriumRecord('eq2')
G2=oc.getGibbsEnergy() # a scalar
print('G={0:e} (vs. without grid-minimizer: {1:e})'.format(G,G2))
oc.changeEquilibriumRecord()


# ### Retrieving chemical potentials

# In[11]:


mu=oc.getChemicalPotentials() # a dictionary (keys are element names, values are chemical potentials)
print('mu_FE= ',mu['FE'])


# ### Retrieving equilibrium phases composition

# In[12]:


phasesAtEquilibrium=oc.getPhasesAtEquilibrium() # a container class defined in pyOC.py
phaseElementComposition=phasesAtEquilibrium.getPhaseElementComposition() # a dictionary (keys are phase names) of dictionaries (keys are the element names, values are molar fractions)
print('n_C in FCC_A1#1 = ',phaseElementComposition['FCC_A1#1']['C'])
phaseSites=phasesAtEquilibrium.getPhaseSites() # a dictionary (keys are phase names, values are arrays of number of sites whose sizes depend on the number of sublattices)
print('a_i in FCC_A1#1 = ',phaseSites['FCC_A1#1'])
phaseConstituentComposition=phasesAtEquilibrium.getPhaseConstituentComposition() # a dictionary (keys are phase names) of dictionaries per sublattice (keys are the constituent names, values are molar fractions)
print('y_V^0 in liquid = ',phaseConstituentComposition['FCC_A1#1']['sublattice 0']['V'])


# ### Retrieving constituent composition (this info is 'updated' each time the getPhasesAtEquilibrium method is called by adding constituents that were not present before)

# In[13]:


constituentsDescription = oc.getConstituentsDescription() # a dictionary (keys are the constituent names) of dictionaries (keys are 'mass', 'charge' and 'elements', values are respectively the constituent molar mass, the constituent eletrical charge and a dictionary - keys are element names, values are stoichiometric coefficients)
print('m_V = ',constituentsDescription['V']['mass'])
print('q_V = ',constituentsDescription['V']['charge'])
print('stoi^V_V = ',constituentsDescription['V']['elements']['V'])

