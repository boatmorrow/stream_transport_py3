#from StreamTran import *
import StreamTran.StreamTran as st
from pylab import *
import numpy as np
import tracer_tools.cfc_tools as cfc
from tracer_tools.get_atm_conc import convert2aqueous
import pdb

#inputs
tracer='CFC12'
L = 30000
nx = 1.e4
#A = ([0.,200.,2000.,10000.,30000.],[1.,2.,5.,2.,4.])
w = (np.array([0.,30000.]),np.array([10.,10.]))
d = (np.array([0.,30000.]),np.array([1.,1.]))
A = (w[0],w[1]*d[1])
ql = np.ones(10)*.1e-10
Q_us = .1
Q_trib = ([5000.,20000.],[.05,.01],{tracer:[.2,1.]})
#Q_trib = ([0],[0],[0])
C_gw = ([5000.],[0.2])
#C_gw = ([0.,5000.,10000.],[0.2,0.5,0.3])
C_tau = ([5000.],[50.*365])

#make tracer input time series (cfc,sf6,He)
df_cfc = cfc.get_gas_conc()
P=cfc.lapse_rate(2000.)
T = 7. #mean temp
df_cfc = convert2aqueous(df_cfc,T,P=P,addHe=True,addAr39=False,addKr81=False) #cfc should pmol/kg after this un

#test interpolate the gw flow vec
#ql = randn(10)
#ql = ql+abs(min(ql))
q_lin = st.CreateLatFlowVec(ql,nx)
x = np.linspace(0,L,num=int(nx))

#try and initialize for a flow only
sTran = st.StreamTranSim('test')
sTran.L = L
sTran.nx = nx
sTran.wi = w
#sTran.d = d
sTran.Ai = A
sTran.q_lin = ql
sTran.Q_us = Q_us
sTran.Q_trib = Q_trib
sTran.C_t = df_cfc
sTran.tau_flag = 0
sTran.CalcFlow()
figure()
plot(x,sTran.SimRes.Q.value)
xlabel('distance')
ylabel('discharge')

#lets do a transport simulation
#add a tracer
sTran.add_tracer(tracer)
sTran.Tracers[tracer].C_atm = df_cfc[tracer][-1]
sTran.Tracers[tracer].C_us = df_cfc[tracer][-1]
sTran.Tracers[tracer].C_gw = C_gw
sTran.CalcTran()
figure()
plot(x,sTran.SimRes.Q.value)
figure()
plot(x,sTran.SimRes.tracers[tracer].C.value)


print('the interpolated q_lin length is ', len(q_lin))

figure()
plot(x,q_lin)
xlabel('distance')
ylabel('lateral discharge')


#test interpolate material value

Ai = st.InterpolateMatrialValue(A[0],A[1],nx,L)
figure()
plot(x,Ai)
xlabel('distance')
ylabel('area')

print('the interpolated mat val length is ', len(q_lin))

#test interpolate tribs, pumping...
Qtr = st.InterpolateTribPump(Q_trib[0],Q_trib[1],nx,L)
figure()
plot(x,Qtr)
xlabel('distance')
ylabel('tributary flow')

#Interpolate trib concen.
Ctr = st.InterpolateTribPump(Q_trib[0],Q_trib[2][tracer],nx,L)
figure()
plot(x,Ctr)
xlabel('distance')
ylabel('tributary Conc')

#test interpolate groundwater concentration - raw conc.
numSteps = len(ql)
Cgw = st.InterpolateGwConc(C_gw[0],C_gw[1],numSteps,nx,L)
figure()
plot(x,Cgw)
xlabel('distance')
ylabel('GW conc.')

#test interpolate groundwater concentration - mean age.
Cgw = st.InterpolateGwConcTau(C_tau[0],C_tau[1],numSteps,nx,L,df_cfc,'CFC11',df_cfc.index[-1],0.,mod_type='exponential')
figure()
plot(x,Cgw)
xlabel('distance')
ylabel('GW conc.')

show()

