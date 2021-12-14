#Young adult-born granule cell , multi-compertmental, conductance-based neuron model 

import neuron
from neuron import h
from neuron import h, rxd
from neuron import h, gui
import math
import numpy as np
from neuron.units import ms, mV

#Defining compartments as segments.

class YoungAdultBornGranuleCell:
    def __init__(self):
        self._setup_morphology()
        self._setup_biophysics()
    
    def _setup_morphology(self):
        
        ###Topology###
        
        #Soma#
        
        self.soma = h.Section(name='soma', cell=self)
        
        self.soma.nseg = 1
        self.soma.diam = 10
        self.soma.cm = 1 # Membrane capacitance in micro Farads / cm^2  Mongiat et.al 2009
        self.soma.L = 10
        self.soma.Ra = 4000 #Beining et al (2017)
        
        self.soma.push()
        h.pt3dclear()
        h.pt3dadd(0.0, 0, 0.0, 10)
        h.pt3dadd(10.0, 10, 0.0, 10)
        h.pop_section()
        
        #Dendrite Branches
        
        self.dend0 = h.Section(name='dend0', cell=self)      #proximal Dendrite 0
        
        self.dend0.push()
        h.pt3dclear()
        h.pt3dadd(0, 0, 0, 0.09)
        h.pt3dadd(-10, 0, 0, 0.09)
       
        h.pop_section()
        
        self.dend1 = h.Section(name='dend1', cell=self)     #proximal dendrite 1
        
        self.dend1.push()
        h.pt3dclear()
        h.pt3dadd(-10, 0, 0, 0.09)
        h.pt3dadd(-20, -15, 0, 0.09)
        
        self.dend2 = h.Section(name='dend2', cell=self)   #Medial Dendrite 1
        
        self.dend2.push()
        h.pt3dclear()
        h.pt3dadd(-10, 0, 0, 0.09)
        h.pt3dadd(-20, 15, 0, 0.09)
        
        self.dend3 = h.Section(name='dend3', cell=self) #Medial Dendrite 2
        
        self.dend3.push()
        h.pt3dclear()
        h.pt3dadd(-20, -30, 0, 0.09)
        h.pt3dadd(-60, -15, 0, 0.09)
        
        self.dend4 = h.Section(name='dend4', cell=self) # Medial Dendrite 3
        
        self.dend4.push()
        h.pt3dclear()
        h.pt3dadd(-60, 30, 0, 0.09)
        h.pt3dadd(-100, 60, 0, 0.09)
        
        self.dend5 = h.Section(name='dend5', cell=self) #Medial Dendrite 4
        
        self.dend5.push()
        h.pt3dclear()
        h.pt3dadd(-60, 60, 0, 0.09)
        h.pt3dadd(-100, 105, 0, 0.09)
        
        self.dend6 = h.Section(name='dend6', cell=self) #Distal Dendrite 1
        
        self.dend6.push()
        h.pt3dclear()
        h.pt3dadd(-60, -60, 0, 0.09)
        h.pt3dadd(-120, -120, 0, 0.09)
        
        self.dend7 = h.Section(name='dend7', cell=self) #Distal Dendrite 2
        
        self.dend7.push()
        h.pt3dclear()
        h.pt3dadd(-60, 95, 0, 0.09)
        h.pt3dadd(-120, 100, 0, 0.09)
        
        self.dend8 = h.Section(name='dend8', cell=self) #Medial Dendrite 5
        
        self.dend8.push()
        h.pt3dclear()
        h.pt3dadd(-60, -30, 0, 0.09)
        h.pt3dadd(-100, -45, 0, 0.09)
                  
        
        #Soma-Dendrit Connections
        
        self.dend0.connect(self.soma(0))
        self.dend1.connect(self.dend0(1))
        self.dend2.connect(self.dend0(1))
        self.dend3.connect(self.dend1(1))
        self.dend4.connect(self.dend3(1))
        self.dend5.connect(self.dend2(1))
        self.dend6.connect(self.dend3(1))
        self.dend7.connect(self.dend2(1))
        self.dend8.connect(self.dend1(1))
        
        # ###Axon Morp.###
        
        self.axonhl = h.Section(name='axonhl', cell=self)      #Axonhillock-First section after soma
        self.axonhl.nseg = 1
        self.axonhl.diam = 1.5
        self.axonhl.cm = 1
        self.axonhl.L = 10
        self.axonhl.Ra = 200

        #Axon Initial Segment

        self.axonais = h.Section(name= 'axonais', cell=self)   #Axon initiation point with high density of Nav 
        self.axonais.nseg = 1
        self.axonais.diam = 1
        self.axonais.cm = 1
        self.axonais.L = 30
        self.axonais.Ra = 122           #Masoli et al (2015)

        #Axon Prox#

        self.axon1 = h.Section(name='axon1', cell=self)
        self.axon1.nseg = 1
        self.axon1.diam = 0.7
        self.axon1.cm = 1
        self.axon1.L = 100
        self.axon1.Ra = 122  

      
    
        #Bleb#
        self.bleb = h.Section(name='bleb', cell=self)      #bleb
        self.bleb.nseg = 1
        self.bleb.diam = 1.5
        self.bleb.cm = 1
        self.bleb.L = 10
        self.bleb.Ra = 122
        
       ##connection of axon##
        self.soma.connect(self.axonhl(0))
        self.axonhl.connect(self.axonais(0))
        self.axonais.connect(self.axon1(0))
        
        
        
        #connection to bleb#
        self.axon1.connect(self.bleb(0))
        
        
        #Shortcuts of all sec.
        self.somaall = [self.soma]
        self.bleball = [self.bleb]
        self.axonall = [self.axon1]
        self.dendall = [self.dend0, self.dend1, self.dend2, self.dend3, self.dend4, self.dend5, self.dend6, self.dend7, self.dend8]
      
        
        self.proxall = [self.dend1, self.dend2]
        self.medialall =  [self.dend3, self.dend4, self.dend5, self.dend8]
        self.distalall = [self.dend7, self.dend6]
        
    def _setup_biophysics(self):
        
        #Channel mechanisms:

  #pas = passive leak conductance
   #KM = M conductance                  (DGC_M.mod)
   #KA = A-type conductance             (Aradi_KA.mod)
 #fKDR = fast delayed rectifier         (Aradi_KDRf.mod)
 #sKDR = fast delayed rectifier         (Aradi_KDRs.mod)
#CaDep = SK and BK channels             (Aradi_CadepK.mod)
# sAHP = slow AHP conductance           (DGC_sAHP.mod)
   #Na = transient Na conductance       (Aradi_Na.mod)
  # Ca = T, N and L -type Ca channels   (Aradi_Ca.mod)
############################################################################
        
        #Biophysical Parameters#
        
        celsius = 32
        v_init = -62
        ihold = 57 #pA
        amplitude = 450 #pA
        dur1 = 50         //ms
        dur2 = 100
        dur3 = 1000
        #Passive Properties#
        Ra0 = 4000         #Ohm.cm
        cm0 = 1           #uF/cm2
        gpas = 2.5e-05    #S/cm2
        Rm0 = 1e-3/gpas   #kOhm.cm2
        Epas = -75        #mV
        EK = -100         #mV
        #The M conductance
        gMaxon = 20       #pS/um2; 
        kKM = 9          
        v0erevKM = 65
        kVM = 40
        gammaKM = 0.5
        taudivKM = 1
        Dtaumult1 = 30
        Dtaumult2 = 30
        tau0mult = 1
        VshiftKM = 0
       
        #//The SK and BK conductances
        gsksoma = 0.0001
        gskprox = 0.00001
        gskGCLs = 0.00001
        erevSK = EK
        tauskdiv = 2
        erevBK = EK
        BKmult = 0.3

        ##Calcium Channels#
        CaTmult = 1.75
        CaNmult = 1
        CaLmult = 1
        ca0 = .00007
        tauctdiv = 1
        taucadiv = 1
        Vshift = 0
        
        #//The sAHP conductance
        gbarsAHP = 10
        tau1RefsAHP = 400
        tau2sAHP = 200
        c1infsAHP = 0.25
        oinfsAHP = 0.5
        CaRefsAHP = 0.002
        cahsAHP = 0.01
        kcasAHP = 0.001
        nsAHP = 4
        msAHP = 1
        
        #//The delayed rectifier and A conductances
        KDRmult = 0.8
        V0KDR = 23
        taumultKDR = 1
        gKAs = 0.012
        gKAa = 0.004
        
        #//The transient Na conductance
        gNaT_mult = 1.5
        taumultNa = 1
        htaumultNa = 1
        ENa = 42
        
                ###Inserting Channels###  
        
        #Passive for all sec.#
        
    
            
        for sec in self.axonall:
            sec.insert('pas')
            sec.g_pas = gpas
            sec.e_pas = Epas
        
        self.axonais.insert('pas')
        self.axonais.g_pas = gpas
        self.axonais.e_pas = Epas
        
            
        for sec in self.somaall:
            sec.insert('pas')
            sec.g_pas = gpas
            sec.e_pas = Epas
            
        for sec in self.bleball:
            sec.insert('pas')
            sec.g_pas = gpas
            sec.e_pas = Epas
            
        ###Dendrites###
        ##passive##
        
        for sec in self.dendall:
            sec.insert('pas')
            sec.g_pas = 4e-05
            sec.e_pas = Epas
            
        #Proximal branches#
        
        for sec in self.proxall:
            sec.insert('Ca')
            sec.gtcabar_Ca = CaTmult*0.001
            sec.gncabar_Ca = CaNmult*0.001
            sec.glcabar_Ca = CaLmult*0.015*50/100
        
        for sec in self.proxall:
            sec.insert('CadepK')
            sec.gbkbar_CadepK = 0
            sec.gskbar_CadepK = gskprox
        
        for sec in self.proxall:
            sec.insert('Na')
            sec.gbar_Na = gNaT_mult*0.0000013
        
        for sec in self.proxall:
            sec.insert('fKDR')
            sec.gbar_fKDR = 0.004
        
        for sec in self.proxall:
            sec.insert('sKDR')
            sec.gbar_sKDR = 0.003
            
        #medial dendrites#
        
        for sec in self.medialall:
            sec.insert('Ca')
            sec.gtcabar_Ca = CaTmult*0.002
            sec.gncabar_Ca = CaNmult*0.001
            sec.glcabar_Ca = CaLmult*0.001
            
        for sec in self.medialall:
            sec.insert('CadepK')
            sec.gbkbar_CadepK = BKmult*0.0012
            sec.gskbar_CadepK = 0.0
        
        for sec in self.medialall:
            sec.insert('Na')
            sec.gbar_Na = gNaT_mult*0.008
        
        for sec in self.medialall:
            sec.insert('fKDR')
            sec.gbar_fKDR = 0.001
        
        for sec in self.medialall:
            sec.insert('sKDR')
            sec.gbar_sKDR = 0.003
            
        
        #Distal Dendrites#
        
        for sec in self.distalall:
            sec.insert('Ca')
            sec.gtcabar_Ca = CaTmult*0.002
            sec.gncabar_Ca = CaNmult*0.001
            sec.glcabar_Ca = 0.0
            
        for sec in self.distalall:
            sec.insert('CadepK')
            sec.gbkbar_CadepK = BKmult*0.0012
            sec.gskbar_CadepK = 0.0
            
        for sec in self.distalall:
            sec.insert('fKDR')
            sec.gbar_fKDR = 0.001
        
        for sec in self.distalall:
            sec.insert('sKDR')
            sec.gbar_sKDR = 0.004
            
            
                                        
                
        ##Axon##
        
        for sec in self.axonall:
            sec.insert('Na')
            sec.gbar_Na = gNaT_mult*0.15625 ### 0.04  (S-H)
            
        self.axonais.insert('Na')
        self.axonais.gbar_Na = gNaT_mult*0.0625
        
        self.axonais.insert('fKDR')
        self.axonais.gbar_fKDR = 0.0032
        
        for sec in self.axonall:
            sec.insert('fKDR')
            sec.gbar_fKDR = 0.00128
            
        for sec in self.axonall:
            sec.insert('KA')
            sec.gbar_KA = gKAa ## 0.004
        
        for sec in self.axonall:
            sec.insert('KM')
            sec.k_KM = kKM
            sec.v0erev_KM = v0erevKM
            sec.kV_KM = kVM
            sec.gamma_KM = gammaKM
            sec.taudiv_KM = taudivKM
            sec.Dtaumult1_KM = Dtaumult1
            sec.Dtaumult2_KM = Dtaumult2
            sec.Vshift_KM = VshiftKM
            sec.tau0mult_KM = tau0mult
            sec.gbar_KM = gMaxon
        
        #Soma#
        self.soma.insert('Ca')    #Ca channels (T,N,L-type)  
        self.soma.gtcabar_Ca = CaTmult*0.00003
        self.soma.gncabar_Ca = CaNmult*0.003
        self.soma.glcabar_Ca = CaLmult*0.0075 #downregulated
    
       
        self.soma.insert('CadepK') #Calcium dependent Potassium Channels
        self.soma.gbkbar_CadepK = BKmult*0.00003
        self.soma.gskbar_CadepK = 0
        
        self.soma.insert('Na') #Sodium channels
        self.soma.gbar_Na = gNaT_mult*0.078125
        
        self.soma.insert('fKDR')
        self.soma.gbar_fKDR = 0.0016
        
        self.soma.insert('sKDR')
        self.soma.gbar_sKDR = 0.003
        
        self.soma.insert('KA')
        self.soma.gbar_KA = gKAs*0.9
        
        
        self.soma.insert('sAHP')
        self.soma.gbar_sAHP = gbarsAHP
        
    ###Adjustment for density mechanism##
                       
          ##AXONALL#
    
        for sec in self.axonall:
            if h.ismembrane('KM'):
                sec.erev_KM = EK
        
        for sec in self.axonall:
            if h.ismembrane('Na'):
                sec.ena = ENa
                sec.taumult_Na = taumultNa
        
        for sec in self.somaall:
            if h.ismembrane('Na'):
                sec.ena = ENa
                sec.taumult_Na = taumultNa
                
        
                
        for sec in self.somaall:
            if h.ismembrane('sAHP'):
                sec.ek = EK
             
                
    def __repr__(self):
        return 'YoungAdultBornGranuleCell'