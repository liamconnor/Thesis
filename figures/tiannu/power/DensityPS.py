import numpy as np
import matplotlib
matplotlib.use('PDF')
import pylab as py

#Redshift of ps
zs = '0.010'

#Cosmology parameters
h=0.67
Oc = 0.27+0.05 #Oc + Ob
On = 0.05/93.14/h**2
fc = Oc/(Oc+On)
fn = On/(Oc+On)

#Simulation power spectra
suf='.dat'
dir = './ps/dmnu/'
fname = 'ngp_density_delta'+zs+'_'
dd = np.genfromtxt(dir+fname+'dmdm'+suf)
nn = np.genfromtxt(dir+fname+'nunu'+suf)
dn = np.genfromtxt(dir+fname+'dmnu'+suf)
k = dd[:,1]

dir = './ps/dmon/'
do = np.genfromtxt(dir+fname+'dmdm'+suf)

##Compute poisson noise
L = 1200.0
nc = 13824.0
ni = (L/nc)**(3)
pn = k**3 * ni / 2.0 / np.pi**2

##Compute power spectra
ddm = (dd[:,2]-pn)/dd[:,3]
nnm = (nn[:,2]-pn)/nn[:,3]
dnm = (dn[:,2])#/dn[:,3]

mp = (ddm*fc**2+2*fc*fn*dnm+fn**2*nnm)

#Class transfer functions
dir = './tf/'
lmps = np.genfromtxt(dir+'th2_dtf_z1_pk.dat')
lmps0 = np.genfromtxt(dir+'th2_dtf_nonu_z1_pk.dat')
nlmps = np.genfromtxt(dir+'th2_dtf_nl_z1_pk_nl.dat')
nlmps0 = np.genfromtxt(dir+'th2_dtf_nl_nonu_z1_pk_nl.dat')
tf = np.genfromtxt(dir+'th2_dtf_z1_tk.dat')
tf0 = np.genfromtxt(dir+'th2_dtf_nonu_z1_tk.dat')

tcol=6
ncol=5
dcol=1
#Plot individual power spectra
py.figure(1)

py.loglog(k,dd[:,2]/dd[:,3],'--',label='DMxDM',color='black')
py.loglog(k,dn[:,2]/dn[:,3],'--',label='DMxNU',color='green')
py.loglog(k,nn[:,2]/nn[:,3],'--',label='NUxNU',color='red')


py.loglog(k,pn,':',label='PN',color='grey')
py.loglog(k,pn/dd[:,3],'--',label='PN',color='grey')

py.loglog(k,(dd[:,2]-pn)/dd[:,3],'-',label='DMxDM',color='black')
py.loglog(k,(dn[:,2])/dn[:,3],'-',label='DMxNU',color='green')
py.loglog(k,(nn[:,2]-pn)/nn[:,3],'-',label='NUxNU',color='red')

py.xlabel(r'k (h/Mpc)')
py.axes().set_xlim((2*np.pi/1200,20))
py.ylabel(r'$\Delta(k)$')
py.axes().set_ylim((1e-8,1000))
#py.legend(loc='best',fancybox=True,frameon=False)
py.savefig('./ps.pdf')
py.close(1)

#Plot total matter power spectra
py.figure(1)

py.loglog(k,mp,'--',label='Total MP',color='purple') 
py.loglog(k,dd[:,2]/dd[:,3],':',label='DMxDM',color='blue')
py.loglog(k,do[:,2]/do[:,3],':',label='DOxDO',color='red')

py.xlabel(r'k (h/Mpc)')
py.ylabel(r'$\Delta(k)$')
py.savefig('./mps.pdf')
py.close(1)

