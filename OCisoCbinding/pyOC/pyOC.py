import logging, sys
import json
import numpy as np
import rawpyOC
from rawpyOC import rawopencalphad as oc
from enum import IntEnum

class PhaseStatus(IntEnum):
	Suspended = -3
	Dormant   = -2
	Entered   =  0
	Fixed     =  2
	
class GridMinimizerStatus(IntEnum):
	On  = 0
	Off = -1

class PhasesAtEquilibrium(object):

	def __init__(self, phaseMolarAmounts, phaseElementComposition, phaseSites, phaseConstituentComposition):
		self.__phaseMolarAmounts=phaseMolarAmounts
		self.__phaseElementComposition=phaseElementComposition
		self.__phaseSites=phaseSites
		self.__phaseConstituentComposition=phaseConstituentComposition
	
	def getPhaseMolarAmounts(self):
		return self.__phaseMolarAmounts
		
	def getPhaseElementComposition(self):
		return self.__phaseElementComposition
		
	def getPhaseSites(self):
		return self.__phaseSites
		
	def getPhaseConstituentComposition(self):
		return self.__phaseConstituentComposition
		
	def __str__(self):
		return 'phase molar amounts:\n'+json.dumps(self.__phaseMolarAmounts, indent=4)+'\nphase element composition:\n'+json.dumps(self.__phaseElementComposition, indent=4)+'\nphase sites:\n'+json.dumps(self.__phaseSites, indent=4)+'\nphase constituent composition:\n'+json.dumps(self.__phaseConstituentComposition, indent=4)

class OpenCalphad(object):

	_defaultEquilibriumName='default equilibrium'
	_maxNbPhases=400
	_maxNbElements=100
	_maxNbSublattices=10
	_maxNbConstituents=100

	def __init__(self):
		self.__equilibriumNamesInOC = {}
		self.__constituentsDescription = {}
		self.__logger = logging.getLogger('OpenCalphad')
		self.__logger .setLevel(logging.INFO)
		ch = logging.StreamHandler()
		ch.setStream(sys.stderr)
		ch.setLevel(logging.INFO)
		formatter = logging.Formatter('(%(name)s being verbose): %(message)s')
		ch.setFormatter(formatter)
		self.__logger.addHandler(ch)
				
	def setVerbosity(self, isVerbose):
		if (isVerbose):
			level= logging.DEBUG
		else:
			level = logging.INFO
		self.__logger.setLevel(level)
		for handler in self.__logger.handlers:
		    handler.setLevel(level)
		oc.pytqquiet(~isVerbose)
	
	def raw(self):
		return oc
		
	def eq(self):
		return self.__eq

	def readtdb(self,tdbFilePath, elements=None):
		self.__eq = oc.pytqini(1)
		if self.__logger.getEffectiveLevel() is not logging.DEBUG:
			oc.pytqquiet(True)
		self.__eqName = OpenCalphad._defaultEquilibriumName
		eqNameInOC='%s' % self.__eqName.upper().replace(' ','_')
		self.__equilibriumNamesInOC[self.__eqName]=eqNameInOC
		self.__logger.debug('reading %s', tdbFilePath)
		if elements is None:
			oc.pytqrfil(tdbFilePath,self.__eq)
		else:
			oc.pytqrpfil(tdbFilePath,len(elements),(''.join(['%-2s']*len(elements)))%elements,self.__eq)
		comp = oc.pytqgcom(self.__eq)
		self.__componentNames=np.array([comp.compnames[:,i].tostring().decode().strip() for i in range(comp.n)])
		self.__nbComponents=comp.n
		self.__logger.debug('component (%d) names: %s', self.__nbComponents, self.__componentNames)

	def getComponentNames(self):
		return self.__componentNames
		
	def getConstituentsDescription(self):
		self.__logger.debug('constituents description:\n'+json.dumps(self.__constituentsDescription, indent=4))
		return self.__constituentsDescription
		
	def setPhasesStatus(self, phaseNames, phaseStatus, phaseAmount=0.0):
		phaseList=(' '.join(['%s']*len(phaseNames))) % phaseNames
		self.__logger.debug('modifying phases %s to status %s', phaseList, phaseStatus)
		oc.pytqphsts2(phaseList,phaseStatus,phaseAmount,self.__eq)
		
	def setTemperature(self,temperature):
		self.__logger.debug('setting temperature to %5.2f K', temperature)
		oc.pytqsetc('T',0,0,temperature,self.__eq)
		
	def setPressure(self,pressure):
		self.__logger.debug('setting pressure to %3.2e Pa', pressure)
		oc.pytqsetc('P',0,0,pressure,self.__eq)
		
	def setElementMolarAmounts(self,elementMolarAmounts):
		for el, n in elementMolarAmounts.items():
			i=np.where(self.__componentNames==el)[0]
			self.__logger.debug('setting molar amount %5.4f for element %s (%d)',n,el,i)
			oc.pytqsetc('N',i+1,0,n,self.__eq)
	
	def changeEquilibriumRecord(self,eqName=None,copiedEqName=None):
		if eqName is None:
			eqName=OpenCalphad._defaultEquilibriumName
		if copiedEqName is None:
			copiedEqName=OpenCalphad._defaultEquilibriumName
		eqNameInOC=self.__equilibriumNamesInOC.get(eqName,'')
		if (eqNameInOC==''):
			eqNameInOC='%s' % eqName.upper().replace(' ','_')
			self.__equilibriumNamesInOC[eqName]=eqNameInOC
			eq=oc.pytqselceq(self.__equilibriumNamesInOC[copiedEqName])
			self.__logger.debug('creating and selecting new equilibrium record \'%s\' (\'%s\')', eqName, eqNameInOC)
			iCopiedEq,self.__eq=oc.pytqcceq(eqNameInOC,eq)
		else:
			self.__logger.debug('selecting equilibrium record \'%s\' (\'%s\')', eqName, eqNameInOC)
			self.__eq=oc.pytqselceq(eqNameInOC)
		self.__eqName = eqName
	
	def calculateEquilibrium(self,gridMinimizerStatus=GridMinimizerStatus.On):
		self.__logger.debug('calculating equilibrium with grid minimizer %s', gridMinimizerStatus)
		if (self.__logger.isEnabledFor(level=logging.DEBUG)):
			oc.pytqlc(6,self.__eq)
		oc.pytqce('',gridMinimizerStatus,0,0.0,self.__eq)
		if (self.__logger.isEnabledFor(level=logging.DEBUG)):
			oc.pytqlr(6,self.__eq)

	def getErrorCode(self):
		return oc.pygeterr()

	def resetErrorCode(self):
		return oc.pyseterr(0)
			
	def getScalarResult(self,symbol):
		value=np.empty(1)
		oc.pytqgetv(symbol,0,0,1,value,self.__eq)
		self.__logger.debug('retrieving %s: %e',symbol,value[0])
		return value[0]
			
	def getGibbsEnergy(self):
		return self.getScalarResult('G')
		
	def getComponentAssociatedResult(self, symbol):
		values={}
		value=np.empty(1)
		for i in range(self.__nbComponents):
			oc.pytqgetv(symbol,i+1,0,1,value,self.__eq)
			values[self.__componentNames[i]] = value[0]
		self.__logger.debug('retrieving %s:\n%s',symbol,json.dumps(values, indent=4))
		return values

	def getChemicalPotentials(self):
		return self.getComponentAssociatedResult('MU')	
		
	def getPhasesAtEquilibrium(self):
		tmpNbPhases=np.empty(OpenCalphad._maxNbPhases)
		tmpNbElements=np.empty(OpenCalphad._maxNbElements)
		tmpiNbElements=np.empty(OpenCalphad._maxNbElements, dtype=np.int32)
		tmpiNbSublattices=np.empty(OpenCalphad._maxNbSublattices, dtype=np.int32)
		tmpNbSublattices=np.empty(OpenCalphad._maxNbSublattices)
		tmpiNbConstituents=np.empty(OpenCalphad._maxNbSublattices*OpenCalphad._maxNbConstituents, dtype=np.int32)
		tmpNbConstituents=np.empty(OpenCalphad._maxNbSublattices*OpenCalphad._maxNbConstituents)
		tmp5=np.empty(5)
		#
		nbPhases=oc.pytqgetv('NP',-1,0,OpenCalphad._maxNbPhases,tmpNbPhases,self.__eq)
		phaseMolarAmounts={}
		phaseElementComposition={}
		phaseSites={}
		phaseConstituentComposition={}
		for i in range(nbPhases):
			if (tmpNbPhases[i]>0.0):
				# 
				phaseName=oc.pytqgpn(i+1,self.__eq).decode().strip()
				phaseMolarAmounts[phaseName] = tmpNbPhases[i]
				#
				iph, ics = oc.pytqgpi2(phaseName,self.__eq)
				#
				phaseElementComposition[phaseName] = {}
				nbElements=oc.pytqgetv('X',i+1,-1,OpenCalphad._maxNbElements,tmpNbElements,self.__eq)
				for j in range(nbElements):
					phaseElementComposition[phaseName][self.__componentNames[j]]=tmpNbElements[j]
				#
				phaseConstituentComposition[phaseName]={}
				phaseName=oc.pytqgpn(i+1,self.__eq).decode().strip()
				nbSublattices = oc.pytqgphc1(i+1, tmpiNbSublattices, tmpiNbConstituents, tmpNbConstituents, tmpNbSublattices, tmp5, self.__eq)
				phaseSites[phaseName] = tmpNbSublattices[0:nbSublattices].tolist()
				count = 0
				for j in range(nbSublattices):
					sublatticeConstituentComposition = {}
					offset = count
					for k in np.nditer(tmpiNbConstituents[offset:offset+tmpiNbSublattices[j]]):
						constituentName = oc.pytqgpcn2(iph,count+1).decode().strip()
						sublatticeConstituentComposition[constituentName] = tmpNbConstituents[count]
						if not constituentName in self.__constituentsDescription:
							nspel, smass, qsp = oc.pytqgpcs(k, tmpiNbElements, tmpNbElements)
							self.__constituentsDescription[constituentName] = {}
							self.__constituentsDescription[constituentName]['mass'] = smass
							self.__constituentsDescription[constituentName]['charge'] = qsp
							self.__constituentsDescription[constituentName]['elements'] = { self.__componentNames[tmpiNbElements[l]-1] : tmpNbElements[l] for l in range(nspel) if (tmpiNbElements[l]>0)}
						count += 1
					if (nbSublattices==1):
						phaseConstituentComposition[phaseName] = sublatticeConstituentComposition
					else:
						phaseConstituentComposition[phaseName]["sublattice {0:d}".format(j)] = sublatticeConstituentComposition
		phasesAtEquilibrium=PhasesAtEquilibrium(phaseMolarAmounts,phaseElementComposition,phaseSites,phaseConstituentComposition)
		self.__logger.debug('phases at equilibrium:\n%s',phasesAtEquilibrium)
		return phasesAtEquilibrium
		
opencalphad = OpenCalphad()
