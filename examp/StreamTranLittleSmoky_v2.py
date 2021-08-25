
#######################################################################
#######################################################################
# Example input deck - runs a two stage simulation.  First basic model 
# domain and parameters are defined.  After setting up the simulation,
# groundwater inflow is calculated using Rn and discharge measurments.
# The groundater age is then estimated using an exponential (or dispersion) 
# mixing model from the SF6 measurments.
# This input requires the tracer_tools package to be in the python path.
#######################################################################
#######################################################################
from StreamTran import *
from pylab import *
import numpy as np
import cfc_tools as cfc
from get_atm_conc import convert2aqueous
import pdb
import datetime as dt
import noble_gas_tools as ng
import matplotlib as mpl
import pickle as p

#mpl.rcParams['figure.figsize']=[10.0,8.0];
#mpl.rcParams['text.usetex'] = True;
#mpl.rcParams['font.size'] = 20.0;

#initialize the stream transport simulation
sTran = StreamTranSim('LittleSmoky')
FieldData = {}
#cacluate groundwater inflow?
run_tran = 1
#calculate groundwater age?
calc_age = 1


#######################################################################
#######################################################################
#inputs
#######################################################################
#######################################################################

#stream temperature
T = 10.475 #mean water temperature for samples
T_gw = 4.8 #mean recharge temperature from paskapoo GW
#stream elevation
E = 760. #mean sampling elevation
E_gw = 878. #mean sampling elevation for now

#initialize tracer input time series.  
#make tracer input time series (cfc,sf6,He) using atmospheric equilibrium
df_cfc_atm = cfc.get_gas_conc()
P = cfc.lapse_rate(E)
P_gw = cfc.lapse_rate(E_gw)
df_cfc = convert2aqueous(df_cfc_atm,T,P=P,addHe=True,addAr39=False,addKr81=False) #cfc should pmol/kg after this un
df_cfc_gw = convert2aqueous(df_cfc_atm,T_gw,P=P_gw,addHe=True,addAr39=False,addKr81=False) #cfc should pmol/kg after this un

#sampling time - don't worry about this right now.
t_i = df_cfc.index[-1]

##################################
#discharge field data - gauged data
xd = np.array([0.01,149.2])*1000.
yd = np.array([1.0,4.5])
DischData=(xd,yd)

##################################
#groundwater age data - one age for all for now
xd = np.array([40.])*1000.
yd = np.array([4.*365])  #mean age of groundwater wells sampled
C_tau = (xd,yd)
sTran.C_tau = C_tau


#######################################################################
#tracer info
#######################################################################
#BELOW TRACERS ARE ADDED TO THE SIMULATION.  ADDITIONAL TRACERS CAN BE ADDED OR REMOVED AT WILL.  
#TRACERS ARE KEYED BY THEIR NAME, SO THE NAME IS IMPORTANT, THE ORDER IS NOT...
##################################
#tracer - 1
tracer = 'Rn'
k_exch = 2. #m/d  #THIS IS ESTIMATED BETTER USING RAYMOND ET AL. 2012 JUST BEFORE SIMULATION...

k_exch = k_exch/60./60./24. #m/s
lamma = 3.8235/60./60./24. #s-1
C_atm = 0.
error_perc = 0.05

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([0.119145552,0.170615371,0.20523669,0.113551722,0.075396566,0.115146444,0.05385325,0.053387195]) #bq/l
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
#yd = np.array([14.38]) #average of groundwater data
yd = np.array([7.62]) #min of groundwater data
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc

##################################
#tracer - 2
tracer = 'SF6'
k_exch = 2. #m/d 
k_exch = k_exch/60./60./24. #m/s
lamma = 0. #s-1
C_atm = df_cfc[tracer][-1]
has_age = 1
error_perc = 0.15

print('SF6 equlibrium conc. is ', '%0.2g'%C_atm ,' fMol/g')

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([2.436300826,2.461125903,2.668872477,2.607161224,2.471835049,2.53411851,2.616043849,2.335642312]) #fmol/kg
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
yd = np.array([0.36]) #avg of groundwater data
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc
sTran.Tracers[tracer].has_age = has_age

##################################
#tracer - 3
tracer = 'd18O'
k_exch = 0. #m/s
lamma = 0. #s-1
C_atm = 0.
error_perc = 0.05

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([-18.52,-18.4,-18.25,-18.29,-18.07,-18.14,-17.91,-17.4])
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
yd = np.array([-19.71])
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc

##################################
#tracer - 4
tracer = 'dD'
k_exch = 0. #m/s
lamma = 0. #s-1
C_atm = 0.
error_perc = 0.05

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([-147.79,-146.56,-146.65,-147.34,-144.93,-145.54,-145.08,-142.5])
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
yd = np.array([-153.94])
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc

##################################
#tracer - 5
tracer = 'He4'
#tracer = 'CFC12'
k_exch = 8./60./60./24. #m/s
lamma = 0. #s-1
P=ng.lapse_rate(E)  #range 1121 us to 1070 downstream
C_atm = ng.equil_conc('He',T,P=P)
print('He equlibrium conc. is ', '%0.2g'%C_atm ,' ccSTP/g')
#C_atm = df_cfc[tracer][-1]
error_perc = 0.03

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([4.08E-08,4.13E-08,4.19E-08,4.28E-08,4.30E-08,4.23E-08,4.25E-08,4.22E-08])
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
yd = np.array([1.10e-7])
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc

##################################
#tracer - 6
tracer = 'H3'
k_exch = 0. #m/s
t_half = 12.32/356./24./60./60. #s
lamma = np.log(2.)/t_half #s-1
C_atm = 10. #TU
error_perc = 0.03

#tracer field data 

#river data
xd = np.array([0.000,31.000,54.800,65.100,99.300,121.200,149.200,172.000])*1000.
yd = np.array([7.85,7.69,7.82,7.21,7.25,7.43,7.59,7.60])
FieldData[tracer]=(xd,yd)

#groundwater data
xd = np.array([40.])*1000.
yd = np.array([0.41])  #average of groundwater wells
C_gw = (xd,yd)

#add tracer to simulation
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = C_atm
sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
sTran.Tracers[tracer].C_gw = C_gw
sTran.Tracers[tracer].k = k_exch
sTran.Tracers[tracer].error_perc = error_perc

#######################################################################
#model geometry
#######################################################################
L = 172000.  #m
nx = 1.e4   #number of cells
E = 2e-3/60./60./24. #m/s Q_us = 1.0

# geometery visually estimated at sampling locations
w = (np.array([0.,31.,54.8,65.1,99.3,121.2,149.2,172.])*1000.,\
np.array([25.,30.,30.,35.,50.,45.,40.,60.]))
d = (np.array([0.,31.,54.8,65.1,99.3,121.2,149.2,172.])*1000.,\
np.array([1.,1.,1.,1.,1.,1.,1.,1.]))
A = (w[0],w[1]*d[1])

#Calculate the gas exchange coefficients from Raymond 2012 - equation 7.
E_us = 956.
E_ds = 645.
S = (E_us-E_ds)/d[0][-1]
Q_av = Q_us
A_av = A[1].mean()
V_av = Q_av/A_av
D_av = d[1].mean()
k600 = ng.k600(Q_av,V_av,D_av,S,eqn=7)
Sch_He = ng.schmidt('He',T)
Sch_Rn = ng.schmidt('Rn',T)
k_Rn = (Sch_Rn/600.)**-0.5*k600
k_He = (Sch_He/600.)**-0.5*k600
#k_Rn = 3.85 #lower MCMC 95%
#k_Rn = 5.42 #lower MCMC 95%
sTran.Tracers['Rn'].k = k_Rn/24./60./60.
sTran.Tracers['He4'].k = k_He/24./60./60.
print('Rn Gas Exchange = ', sTran.Tracers['Rn'].k*60*60*24, ' m/d')
print('He Gas Exchange = ', sTran.Tracers['He4'].k*60*60*24, ' m/d')


##################################
#Tributary info
#tributary ([x],[disch],{tracer:[C]]})
Q_trib = (np.array([53.8,56.3,151.2,154.2])*1000,[0.1,0.4,1.2,0.6],{'Rn':[0.2,0.2,0.1,0.1],'SF6':[df_cfc['SF6'][-1],df_cfc['SF6'][-1],df_cfc['SF6'][-1],df_cfc['SF6'][-1]],'d18O':[-18.21,-18.21,-18.21,-18.21],'dD':[-146.0,-146.0,-146.0,-146.0],'He4':[4.2e-8,4.2e-8,4.2e-8,4.2e-8],'H3':[7.5,7.5,7.5,7.5]}) 

##################################
#Groundwater inflows
#inflows
#ql = np.array([5.,5.,5.,5.,5.,10.,10.,10.])
#ql = ql/60./60./24./12.  #the eight is a fudge factor for moving between cookie's and my units.  Approximate width - could be better if you used an interpolated width but works for now...
ql = 1.e-7*np.ones(7)

##################################
#end model input

#######################################################################
#Parameterize Simulation
#######################################################################
sTran.L = L
sTran.nx = nx
sTran.wi = w
sTran.Ai = A
sTran.q_lin = ql
sTran.Q_us = Q_us
sTran.Q_trib = Q_trib
sTran.C_t = df_cfc
sTran.tau_flag = 1
sTran.t_i = t_i
sTran.age_dist = 'dispersion'

#######################################################################
#Run Simulation
#######################################################################
#sTran.CalcTran()
print('Calculating groundwater inflow')
tracer_list=['Rn']
if run_tran:
    sTran.CalcGwInflow(tracer_list,FieldData,DischData,DischargeError=0.03)
    f = file('LittleSmokyRiverSimulationResults_v2.pkl','w')
    p.dump(sTran,f)
else:
    f = file('LittleSmokyRiverSimulationResults_v2.pkl','r')
    sTran = p.load(f)
    f.close()

#######################################################################
#visulization the solution
#######################################################################

#######################################################################
#Figure 1

fig1 = figure()

##################################
#discharge
error_perc=0.05
x = np.linspace(0,L,num=nx)
subplot(3,1,1)
plot(x/1000.,sTran.SimRes.Q.value,'b',lw=2,alpha=0.6,label='modeled')

#field data
errorbar(DischData[0]/1000.,DischData[1],yerr=error_perc*DischData[1],fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
ylabel('discharge')
legend()

##################################
#Rn
tracer = 'Rn'
error_perc = 0.15
subplot(3,1,2)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',alpha=0.6,lw=2,label='modeled')
ylabel(tracer)

#plot field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
legend()

##################################
#Groundwater inflow
subplot(3,1,3)
plot(x/1000.,sTran.SimRes.q_lin.value,lw=2,alpha=0.6,label='Groundwater Disch.')
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
grid(b=True,which='major',linestyle='--')
xlabel('distance (km)')
ylabel('lateral discharge')

tight_layout()
fig1.savefig('PhysicalParamsFigure.png',dpi=200,bbox_inches='tight')

#######################################################################
#Figure 2

fig2 = figure()

##################################
#Helium
error_perc=0.10
tracer = 'He4'
x = np.linspace(0,L,num=nx)
subplot(3,1,1)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',lw=2,alpha=0.6,label='modeled')

#field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
ylabel(tracer)
legend()

##################################
#Rn
tracer = 'H3'
error_perc = 0.03
subplot(3,1,2)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',alpha=0.6,lw=2,label='modeled')
ylabel(tracer)

#plot field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
legend()

##################################
#Groundwater inflow
subplot(3,1,3)
plot(x/1000.,sTran.SimRes.q_lin.value,lw=2,alpha=0.6,label='Groundwater Disch.')
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
grid(b=True,which='major',linestyle='--')
xlabel('distance (km)')
ylabel('lateral discharge')

tight_layout()
fig2.savefig('HeH3.png',dpi=200,bbox_inches='tight')

#######################################################################
#Figure 3

fig3 = figure()

##################################
#dD
error_perc=0.03
tracer = 'dD'
x = np.linspace(0,L,num=nx)
subplot(3,1,1)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',lw=2,alpha=0.6,label='modeled')

#field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
ylabel(tracer)
legend()

##################################
#d18O
tracer = 'd18O'
error_perc = 0.03
subplot(3,1,2)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',alpha=0.6,lw=2,label='modeled')
ylabel(tracer)

#plot field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')

##################################
#Groundwater inflow
subplot(3,1,3)
plot(x/1000.,sTran.SimRes.q_lin.value,lw=2,alpha=0.6,label='Groundwater Disch.')
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
grid(b=True,which='major',linestyle='--')
xlabel('distance (km)')
ylabel('lateral discharge')

tight_layout()
fig3.savefig('IsotopeFigure.png',dpi=200,bbox_inches='tight')

#######################################################################
#Figure 4

fig4 = figure()

##################################
#dD
error_perc=0.03
tracer = 'H3'
x = np.linspace(0,L,num=nx)
subplot(2,1,1)
plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',lw=2,alpha=0.6,label='modeled')

#field data
errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
grid(b=True,which='major',linestyle='--')
ylabel(tracer)
legend()

##################################
#Groundwater inflow
subplot(3,1,3)
plot(x/1000.,sTran.SimRes.q_lin.value,lw=2,alpha=0.6,label='Groundwater Disch.')
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
grid(b=True,which='major',linestyle='--')
xlabel('distance (km)')
ylabel('lateral discharge')

fig4.savefig('TritiumFigure.png',dpi=200,bbox_inches='tight')

#geometery
figure()
subplot(2,1,1)
plot(x/1000.,sTran.SimRes.w.value)
ylabel('width')

subplot(2,1,2)
plot(x/1000.,sTran.SimRes.d.value)
ylabel('depth(m)')
xlabel('distance (km)')

#######################################################################
#Estimate Age
#######################################################################
print('Estimating mean age')
tracer_list=['Rn','SF6']
if calc_age:
    sTran.CalcGwAge(tracer_list,FieldData,DischData,DischargeError=0.03)
    f = file('LittleSmokyRiverSimulationResultsAge_v2.pkl','w')
    p.dump(sTran,f)
else:
    f = file('LittleSmokyRiverSimulationResultsAge_v2.pkl','r')
    sTran = p.load(f)
    f.close()

#######################################################################
#visulization the solution
#######################################################################

#######################################################################
#Age Figures
for tracer in sTran.Tracers.keys():
    if sTran.Tracers[tracer].has_age:
        figure()

        ##################################
        #tracer
        error_perc=0.10
        x = np.linspace(0,L,num=nx)
        subplot(3,1,1)
        plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'b',lw=2,alpha=0.6,label='modeled')

        #field data
        errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='rs',mec='k',label='Data')
        grid(b=True,which='major',linestyle='--')
        ylabel(tracer+' '+sTran.Tracers[tracer].unit)
        ylim(2.0,3.2)
        legend()

        ##################################
        #Groundwater mean age
        error_perc = 0.03
        subplot(3,1,2)
        plot(x/1000.,sTran.SimRes.Ctaui/365.,'b',alpha=0.6,lw=2,label='modeled age')
        ylabel('GW Mean Age (yrs)')

        ##################################
        #Groundwater inflow
        subplot(3,1,3)
        plot(x/1000.,sTran.SimRes.q_lin.value,lw=2,alpha=0.6,label='Groundwater Disch.')
        ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        grid(b=True,which='major',linestyle='--')
        xlabel('distance (km)')
        ylabel('lateral discharge')

        tight_layout()
        savefig('AgeFigure.png',dpi=200,bbox_inches='tight')
        show()

show()

sTran.SimRes.dump2xl()
