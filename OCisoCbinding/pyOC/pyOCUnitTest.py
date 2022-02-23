import unittest
import numpy as np
import os
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat

@unittest.skipUnless(os.path.exists(os.environ.get('OCPRIVATEDATA','')+'/feouzr.tdb'), 'requires feouzr database')
class test_feouzr(unittest.TestCase):
	def setUp(self):
		oc.setVerbosity(False)
		# tdb filepath
		
		tdbFile=os.environ['OCPRIVATEDATA']+'/feouzr.tdb'
		# reading tdb
		elems=('O', 'U', 'ZR')
		oc.readtdb(tdbFile,elems)
		# set pressure
		oc.setPressure(1E5)

#Some global data, reference state SER ......................:
#T=   3000.00 K (  2726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
#N=   1.0000E+00 moles, B=   1.1041E+02 g, RT=   2.4944E+04 J/mol
#GS= -4.67120E+05 J, GS/N=-4.6712E+05 J/mol, HS=-1.0408E+05 J, SS= 1.210E+02 J/K

#Some data for components ...................................:
#Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
#O                 4.1492E-01  0.41492 -2.6690E+01  2.5636E-12  SER (default)   
#U                 3.4330E-01  0.34330 -1.4184E+01  6.9174E-07  SER (default)   
#ZR                2.4178E-01  0.24178 -1.1513E+01  9.9991E-06  SER (default)   

#Some data for phases .......................................:
#Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
#LIQUID.................. E  2.832E-01  0.00E+00  2.34E-01    1.21  0.00E+00  X:
# U      4.44834E-01  ZR     3.82781E-01  O      1.72385E-01
# Constitution: There are     5 constituents:
# U1           4.90503E-01  O2ZR1        5.71591E-02  O1           3.11648E-07
# ZR1          4.05352E-01  O2U1         4.69863E-02

#LIQUID_AUTO#2........... E  7.168E-01  0.00E+00  3.51E-01    2.04  0.00E+00  X:
# O      5.10761E-01  U      3.03177E-01  ZR     1.86062E-01
# Constitution: There are     5 constituents:
# O2U1         3.87533E-01  U1           2.32157E-01  O1           1.80264E-06
# ZR1          2.45847E-01  O2ZR1        1.34461E-01
	def test_LiquidWithMiscibilityGap(self):
		# set temperature
		oc.setTemperature(3000)
		# set element molar amounts
		elementMolarAmounts = {
			'U' : 0.343298,
			'O' : 0.414924,
			'ZR': 0.241778
		}
		oc.setElementMolarAmounts(elementMolarAmounts)
		# calculate equilibrium
		oc.calculateEquilibrium(gmStat.On)
		# retrieving Gibbs energy
		G=oc.getGibbsEnergy()
		np.testing.assert_allclose(G, -4.67120E+05, rtol=1e-5, atol=0)
		# retrieving mu data
		mu=oc.getChemicalPotentials()
		self.assertListEqual(list(mu.keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(mu.values()), [-2.6690E+01*2.4944E+04, -1.4184E+01*2.4944E+04, -1.1513E+01*2.4944E+04], decimal=-2)
		# retrieving equilibrium phases composition
		phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
		phaseElementComposition=phasesAtEquilibrium.getPhaseElementComposition()
		self.assertListEqual(list(phaseElementComposition.keys()), ['LIQUID#1', 'LIQUID_AUTO#2'])
		self.assertListEqual(list(phaseElementComposition['LIQUID#1'].keys()), ['O', 'U', 'ZR'])
		self.assertListEqual(list(phaseElementComposition['LIQUID_AUTO#2'].keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(phaseElementComposition['LIQUID#1'].values()), [5.10761E-01,3.03177E-01,1.86062E-01], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseElementComposition['LIQUID_AUTO#2'].values()), [1.72385E-01,4.44834E-01,3.82781E-01], decimal=6)
		phaseSites=phasesAtEquilibrium.getPhaseSites()
		self.assertListEqual(list(phaseSites.keys()), ['LIQUID#1', 'LIQUID_AUTO#2'])
		np.testing.assert_array_almost_equal(np.array(list(phaseSites.values())).ravel(), [1.0, 1.0], decimal=6)
		phaseConstituentComposition=phasesAtEquilibrium.getPhaseConstituentComposition()
		self.assertListEqual(list(phaseConstituentComposition.keys()), ['LIQUID#1', 'LIQUID_AUTO#2'])
		self.assertListEqual(list(phaseConstituentComposition['LIQUID#1'].keys()), ['O1', 'O2U1', 'O2ZR1', 'U1', 'ZR1'])
		self.assertListEqual(list(phaseConstituentComposition['LIQUID_AUTO#2'].keys()), ['O1', 'O2U1', 'O2ZR1', 'U1', 'ZR1'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['LIQUID_AUTO#2'].values()), [3.11648E-07,4.69863E-02,5.71591E-02,4.90503E-01,4.05352E-01], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['LIQUID#1'].values()), [1.80264E-06,3.87533E-01,1.34461E-01,2.32157E-01,2.45847E-01], decimal=6)       	
		# retrieving constituent composition
		constituentsDescription = oc.getConstituentsDescription()
		## to be tested
	
#Some global data, reference state SER ......................:
#T=   3000.00 K (  2726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
#N=   1.0000E+00 moles, B=   1.1041E+02 g, RT=   2.4944E+04 J/mol
#GS= -4.66921E+05 J, GS/N=-4.6692E+05 J/mol, HS=-1.0182E+05 J, SS= 1.217E+02 J/K

#Some data for components ...................................:
#Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
#O                 4.1492E-01  0.41492 -2.6734E+01  2.4518E-12  SER (default)   
#U                 3.4330E-01  0.34330 -1.4119E+01  7.3804E-07  SER (default)   
#ZR                2.4178E-01  0.24178 -1.1495E+01  1.0177E-05  SER (default)   

#Some data for phases .......................................:
#Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
#LIQUID.................. E  1.000E+00  0.00E+00  5.85E-01    1.71  0.00E+00  X:
# O      4.14924E-01  U      3.43298E-01  ZR     2.41778E-01
# Constitution: There are     5 constituents:
# U1           3.35427E-01  O2U1         2.51330E-01  O1           1.03405E-06
# ZR1          3.09983E-01  O2ZR1        1.03259E-01
	def test_LiquidWithoutMiscibilityGap(self):
		# set temperature
		oc.setTemperature(3000)
		# set element molar amounts
		elementMolarAmounts = {
			'U' : 0.343298,
			'O' : 0.414924,
			'ZR': 0.241778
		}
		oc.setElementMolarAmounts(elementMolarAmounts)
		# keep only liquid phase
		oc.setPhasesStatus(('* ',),phStat.Suspended)
		oc.setPhasesStatus(('LIQUID',),phStat.Entered, 1.0)
		# calculate equilibrium
		oc.calculateEquilibrium(gmStat.Off)
		# retrieving Gibbs energy
		G=oc.getGibbsEnergy()
		np.testing.assert_allclose(G, -4.66921E+05, rtol=1e-5, atol=0)
		# retrieving mu data
		mu=oc.getChemicalPotentials()
		self.assertListEqual(list(mu.keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(mu.values()), [-2.6734E+01*2.4944E+04, -1.4119E+01*2.4944E+04, -1.1495E+01*2.4944E+04], decimal=-2)
		# retrieving equilibrium phases composition
		phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
		phaseElementComposition=phasesAtEquilibrium.getPhaseElementComposition()
		self.assertListEqual(list(phaseElementComposition.keys()), ['LIQUID'])
		self.assertListEqual(list(phaseElementComposition['LIQUID'].keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(phaseElementComposition['LIQUID'].values()), [4.14924E-01,3.43298E-01,2.41778E-01], decimal=6)
		phaseSites=phasesAtEquilibrium.getPhaseSites()
		self.assertListEqual(list(phaseSites.keys()), ['LIQUID'])
		np.testing.assert_array_almost_equal(np.array(list(phaseSites.values())).ravel(), [1.0], decimal=6)
		phaseConstituentComposition=phasesAtEquilibrium.getPhaseConstituentComposition()
		self.assertListEqual(list(phaseConstituentComposition.keys()), ['LIQUID'])
		self.assertListEqual(list(phaseConstituentComposition['LIQUID'].keys()), ['O1', 'O2U1', 'O2ZR1', 'U1', 'ZR1'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['LIQUID'].values()), [1.03405E-06,2.51330E-01,1.03259E-01,3.35427E-01,3.09983E-01], decimal=6)	

#Some global data, reference state SER ......................:
#T=   3000.00 K (  2726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
#N=   1.0000E+00 moles, B=   1.1041E+02 g, RT=   2.4944E+04 J/mol
#GS= -4.62227E+05 J, GS/N=-4.6223E+05 J/mol, HS=-8.1712E+04 J, SS= 1.268E+02 J/K

#Some data for components ...................................:
#Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
#O                 4.1492E-01  0.41492 -2.6954E+01  1.9681E-12  SER (default)   
#U                 3.4330E-01  0.34330 -1.3372E+01  1.5588E-06  SER (default)   
#ZR                2.4178E-01  0.24178 -1.1402E+01  1.1178E-05  SER (default)   

#Some data for phases .......................................:
#Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
#C1_FCC.................. E  3.461E-01  0.00E+00  1.42E-01    2.43  0.00E+00  X:
# O      5.88380E-01  U      3.35030E-01  ZR     7.65900E-02
# Constitution: There are     5 constituents:
# O2U1         6.34573E-01  ZR1          1.05930E-01  O1           5.93474E-07
# U1           1.79357E-01  O2ZR1        8.01390E-02

#HCP_A3.................. E  6.539E-01  0.00E+00  4.43E-01    1.48  0.00E+00  X:
# U      3.47675E-01  ZR     3.29218E-01  O      3.23107E-01
# Constitution: Sublattice  1 with     2 constituents and     0.500000 sites
# O1     9.54677E-01  VA     4.53225E-02
#               Sublattice  2 with     2 constituents and     1.000000 sites
# U1     5.13633E-01  ZR1    4.86367E-01

	def test_SuspendingLiquid(self):
		# set temperature
		oc.setTemperature(3000)
		# set element molar amounts
		elementMolarAmounts = {
			'U' : 0.343298,
			'O' : 0.414924,
			'ZR': 0.241778
		}
		oc.setElementMolarAmounts(elementMolarAmounts)
		# keep only liquid phase
		oc.setPhasesStatus(('LIQUID',),phStat.Suspended)
		# calculate equilibrium
		oc.calculateEquilibrium(gmStat.On)
		# retrieving Gibbs energy
		G=oc.getGibbsEnergy()
		np.testing.assert_allclose(G, -4.62227E+05, rtol=1e-5, atol=0)
		# retrieving mu data
		mu=oc.getChemicalPotentials()
		self.assertListEqual(list(mu.keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(mu.values()), [-2.6954E+01*2.4944E+04, -1.3372E+01*2.4944E+04, -1.1402E+01*2.4944E+04], decimal=-2)
		# retrieving equilibrium phases composition
		phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
		phaseElementComposition=phasesAtEquilibrium.getPhaseElementComposition()
		self.assertListEqual(list(phaseElementComposition.keys()), ['C1_FCC', 'HCP_A3'])
		self.assertListEqual(list(phaseElementComposition['C1_FCC'].keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(phaseElementComposition['C1_FCC'].values()), [5.88380E-01,3.35030E-01,7.65900E-02], decimal=6)
		self.assertListEqual(list(phaseElementComposition['HCP_A3'].keys()), ['O', 'U', 'ZR'])
		np.testing.assert_array_almost_equal(list(phaseElementComposition['HCP_A3'].values()), [3.23107E-01,3.47675E-01,3.29218E-01], decimal=6)


@unittest.skipUnless(os.path.exists(os.environ.get('OCPUBLICDATA','')+'/steel7.TDB'), 'requires steel7 database')
class test_steel7(unittest.TestCase):
	def setUp(self):
		oc.setVerbosity(False)
		# tdb filepath
		tdbFile=os.environ['OCPUBLICDATA']+'/steel7.TDB'
		# reading tdb
		oc.readtdb(tdbFile)
		# set pressure
		oc.setPressure(1E5)
		
#Some global data, reference state SER ......................:
#T=   1173.00 K (   899.85 C), P=  1.0000E+05 Pa, V=  6.2399E-06 m3
#N=   1.0000E+00 moles, B=   5.5735E+01 g, RT=   9.7529E+03 J/mol
#G= -5.76738E+04 J, G/N=-5.7674E+04 J/mol, H= 3.1856E+04 J, S= 7.633E+01 J/K
#Some data for components ...................................:
#Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
#C                 4.0000E-02  0.04000 -4.1298E+00  1.6087E-02  SER (default)   
#CR                6.0000E-02  0.06000 -7.0639E+00  8.5546E-04  SER (default)   
#FE                8.3700E-01  0.83700 -5.6716E+00  3.4423E-03  SER (default)   
#MO                5.0000E-02  0.05000 -7.6456E+00  4.7813E-04  SER (default)   
#SI                3.0000E-03  0.00300 -2.0099E+01  1.8668E-09  SER (default)   
#V                 1.0000E-02  0.01000 -1.3475E+01  1.4053E-06  SER (default)   
#Some data for phases .......................................:
#Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
#FCC_A1#1................ E  1.921E-02  5.45E-09  1.04E-02    1.84  0.00E+00  X:
# C      4.56653E-01  MO     1.44330E-01  FE     2.08416E-03  SI     6.67740E-10
# V      3.62669E-01  CR     3.42642E-02
#Constitution: Sublattice  1 with     5 constituents and     1.000000 sites
# V      6.67472E-01  CR     6.30613E-02  FE     3.83579E-03  SI     1.22894E-09
# MO     2.65631E-01
#               Sublattice  2 with     2 constituents and     1.000000 sites
# C      8.40444E-01  VA     1.59556E-01
#FCC_A1_AUTO#2........... E  8.566E-01  6.20E-06  8.45E-01    1.01  0.00E+00  X:
# FE     9.21696E-01  C      1.35199E-02  SI     3.50221E-03  V      1.46572E-03
# CR     5.02328E-02  MO     9.58308E-03
#Constitution: Sublattice  1 with     5 constituents and     1.000000 sites
# FE     9.34328E-01  MO     9.71442E-03  SI     3.55021E-03  V      1.48581E-03
# CR     5.09212E-02
#               Sublattice  2 with     2 constituents and     1.000000 sites
# VA     9.86295E-01  C      1.37051E-02
#M23C6................... E  2.978E-02  3.53E-08  1.03E-03   29.00  0.00E+00  X:
# FE     3.83821E-01  C      2.06897E-01  V      1.28728E-04  SI     0.00000E+00
# CR     3.41462E-01  MO     6.76913E-02
#Constitution: Sublattice  1 with     3 constituents and    20.000000 sites
# FE     5.26893E-01  CR     4.73095E-01  V      1.17857E-05
#               Sublattice  2 with     4 constituents and     3.000000 sites
# MO     6.54349E-01  FE     1.97647E-01  CR     1.46838E-01  V      1.16580E-03
#               Sublattice  3 with     1 constituents and     6.000000 sites
# C      1.00000E+00
#M6C..................... E  9.442E-02  0.00E+00  1.35E-02    7.00  0.00E+00  X:
# MO     3.91919E-01  C      1.42857E-01  V      1.88041E-02  SI     0.00000E+00
# FE     3.81335E-01  CR     6.50840E-02
#Constitution: Sublattice  1 with     1 constituents and     2.000000 sites
# FE     1.00000E+00
#               Sublattice  2 with     1 constituents and     2.000000 sites
# MO     1.00000E+00
#               Sublattice  3 with     4 constituents and     2.000000 sites
# MO     3.71718E-01  FE     3.34674E-01  CR     2.27794E-01  V      6.58144E-02
#               Sublattice  4 with     1 constituents and     1.000000 sites
#C      1.00000E+00
	def test_melting(self):
		# set temperature
		oc.setTemperature(1173)
		# set element molar amounts
		elementMolarAmounts = {
			'C' : 0.04,
			'CR' : 0.06,
			'MO': 0.05,
			'SI': 0.003,
			'V': 0.01,
			'FE': 1.0-0.04-0.06-0.05-0.003-0.01
		}
		oc.setElementMolarAmounts(elementMolarAmounts)
		# calculate equilibrium
		oc.calculateEquilibrium(gmStat.On)
		# retrieving Gibbs energy
		G=oc.getGibbsEnergy()
		np.testing.assert_allclose(G, -5.76738E+04, rtol=1e-5, atol=0)
		# retrieving mu data
		mu=oc.getChemicalPotentials()
		self.assertListEqual(list(mu.keys()), ['C', 'CR', 'FE', 'MO', 'SI', 'V'])
		np.testing.assert_array_almost_equal(list(mu.values()), [-4.1298E+00*9.7529E+03, -7.0639E+00*9.7529E+03, -5.6716E+00*9.7529E+03, -7.6456E+00*9.7529E+03, -2.0099E+01*9.7529E+03, -1.3475E+01*9.7529E+03], decimal=-2)
		# retrieving equilibrium phases composition
		phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
		phaseElementComposition=phasesAtEquilibrium.getPhaseElementComposition()
		self.assertListEqual(list(phaseElementComposition.keys()), ['FCC_A1#1', 'M23C6', 'M6C', 'FCC_A1_AUTO#2'])
		self.assertListEqual(list(phaseElementComposition['FCC_A1#1'].keys()), ['C', 'CR', 'FE', 'MO', 'SI', 'V'])
		self.assertListEqual(list(phaseElementComposition['FCC_A1_AUTO#2'].keys()), ['C', 'CR', 'FE', 'MO', 'SI', 'V'])
		self.assertListEqual(list(phaseElementComposition['M23C6'].keys()), ['C', 'CR', 'FE', 'MO', 'SI', 'V'])
		self.assertListEqual(list(phaseElementComposition['M6C'].keys()), ['C', 'CR', 'FE', 'MO', 'SI', 'V'])
		np.testing.assert_array_almost_equal(list(phaseElementComposition['FCC_A1#1'].values()), [4.56653E-01,3.42642E-02,2.08416E-03,1.44330E-01,6.67740E-10,3.62669E-01], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseElementComposition['FCC_A1_AUTO#2'].values()), [1.35199E-02,5.02328E-02,9.21696E-01,9.58308E-03,3.50221E-03,1.46572E-03], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseElementComposition['M23C6'].values()), [2.06897E-01,3.41462E-01,3.83821E-01,6.76913E-02,0.0,1.28728E-04], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseElementComposition['M6C'].values()), [1.42857E-01,6.50840E-02,3.81335E-01,3.91919E-01,0.0,1.88041E-02], decimal=6)
		phaseSites=phasesAtEquilibrium.getPhaseSites()
		self.assertListEqual(list(phaseSites.keys()), ['FCC_A1#1', 'M23C6', 'M6C', 'FCC_A1_AUTO#2'])
		np.testing.assert_array_almost_equal(phaseSites['FCC_A1#1'], [1.0, 1.0])
		np.testing.assert_array_almost_equal(phaseSites['FCC_A1_AUTO#2'], [1.0, 1.0])
		np.testing.assert_array_almost_equal(phaseSites['M23C6'], [20.0, 3.0, 6.0])
		np.testing.assert_array_almost_equal(phaseSites['M6C'], [2.0, 2.0, 2.0, 1.0])
		phaseConstituentComposition=phasesAtEquilibrium.getPhaseConstituentComposition()
		self.assertListEqual(list(phaseConstituentComposition.keys()), ['FCC_A1#1', 'M23C6', 'M6C', 'FCC_A1_AUTO#2'])
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1#1'].keys()), ['sublattice 0', 'sublattice 1'])
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1#1']['sublattice 0'].keys()), ['CR', 'FE', 'MO', 'SI', 'V'])
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1#1']['sublattice 1'].keys()), ['C', 'VA'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['FCC_A1#1']['sublattice 0'].values()), [6.30613E-02,3.83579E-03,2.65631E-01,1.22894E-09,6.67472E-01], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['FCC_A1#1']['sublattice 1'].values()), [8.40444E-01,1.59556E-01], decimal=6)
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1_AUTO#2'].keys()), ['sublattice 0', 'sublattice 1'])
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1_AUTO#2']['sublattice 0'].keys()), ['CR', 'FE', 'MO', 'SI', 'V'])
		self.assertListEqual(list(phaseConstituentComposition['FCC_A1_AUTO#2']['sublattice 1'].keys()), ['C', 'VA'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['FCC_A1_AUTO#2']['sublattice 0'].values()), [5.09212E-02,9.34328E-01,9.71442E-03,3.55021E-03,1.48581E-03], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['FCC_A1_AUTO#2']['sublattice 1'].values()), [1.37051E-02,9.86295E-01], decimal=6)
		self.assertListEqual(list(phaseConstituentComposition['M23C6'].keys()), ['sublattice 0', 'sublattice 1', 'sublattice 2'])
		self.assertListEqual(list(phaseConstituentComposition['M23C6']['sublattice 0'].keys()), ['CR', 'FE', 'V'])
		self.assertListEqual(list(phaseConstituentComposition['M23C6']['sublattice 1'].keys()), ['CR', 'FE', 'MO', 'V'])
		self.assertListEqual(list(phaseConstituentComposition['M23C6']['sublattice 2'].keys()), ['C'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M23C6']['sublattice 0'].values()), [4.73095E-01,5.26893E-01,1.17857E-05], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M23C6']['sublattice 1'].values()), [1.46838E-01,1.97647E-01,6.54349E-01,1.16580E-03], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M23C6']['sublattice 2'].values()), [1.0], decimal=6)
		self.assertListEqual(list(phaseConstituentComposition['M6C'].keys()), ['sublattice 0', 'sublattice 1', 'sublattice 2', 'sublattice 3'])
		self.assertListEqual(list(phaseConstituentComposition['M6C']['sublattice 0'].keys()), ['FE'])
		self.assertListEqual(list(phaseConstituentComposition['M6C']['sublattice 1'].keys()), ['MO'])
		self.assertListEqual(list(phaseConstituentComposition['M6C']['sublattice 2'].keys()), ['CR', 'FE', 'MO', 'V'])
		self.assertListEqual(list(phaseConstituentComposition['M6C']['sublattice 3'].keys()), ['C'])
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M6C']['sublattice 0'].values()), [1.0], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M6C']['sublattice 1'].values()), [1.0], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M6C']['sublattice 2'].values()), [2.27794E-01,3.34674E-01,3.71718E-01,6.58144E-02], decimal=6)
		np.testing.assert_array_almost_equal(list(phaseConstituentComposition['M6C']['sublattice 3'].values()), [1.0], decimal=6)
		# retrieving constituent composition
		constituentsDescription = oc.getConstituentsDescription()
		## to be tested
		
		# unset temperature
		oc.raw().pytqsetc('T=NONE',0,0,0.0,oc.eq())
		# set liquid phase as fixed at 0.0
		oc.setPhasesStatus(('LIQUID',),phStat.Fixed,0.0)
		# calculate equilibrium associated with melting temperature
		oc.calculateEquilibrium(gmStat.Off)
		meltingTemperature=oc.getScalarResult('T')
		np.testing.assert_allclose(meltingTemperature, 1501.45395, rtol=1e-5, atol=0)
		
if __name__ == '__main__':

    unittest.main()
