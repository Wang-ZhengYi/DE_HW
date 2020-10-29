import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
# from sympy import *




matplotlib.rcParams['mathtext.fontset'] = 'cm'

def main():
	# a2tEuler()
	# Era()
	Universes()
	# tt = cosmo_time(a_t=(Omg.m0/Omg.L0/2)**(1/3))
	# print(tt)
class Omg(object):
	m0 = 0.3
	L0 = 0.7
	r0 = 8.4e-5
	k0 = 0
	h  = 0.68 
	H0 = 1/14.38

def a2tEuler(ascale=np.linspace(1e-8,1,2048),O_m0=0.3,O_L0=0.7,O_r0=8.4e-5,O_k0=0,a_0=1,h=0.68,Ns=2048):
	t    = np.zeros(shape=Ns,dtype='float')
	a    = ascale
	for i in range(1,a.shape[0]):
		t[i] = t[i-1] +  (a[i]-a[i-1])/a[i]/(O_L0+O_k0/a[i]**2+O_m0/a[i]**3+O_r0/a[i]**4)**(1/2)
		# print('loop:{}'.format(i))
	# print('Finished')

	np.savetxt('t_a.dat',t)
	plt.plot(t*9.78/h,a,'b-',lw=0.8,label='a-t')
	plt.grid(b=True, which='major', color='k', linestyle='-.')
	plt.xlabel(r'$t/Gyr$',fontsize=15)
	plt.ylabel(r'$a$',fontsize=15)
	plt.legend(loc='best')
	plt.savefig('t_a.png',dpi=1200)
	# plt.show()
	plt.clf()
	return t



def Era(O_m0=0.3,O_L0=0.7,O_r0=8.4e-5,O_k0=0,a_0=1,h=0.68,Ns=2000):
	T_h  = 9.78/h
	a    = np.linspace(1e-4,1,Ns)
	O_m  = np.zeros(shape=Ns,dtype='float')
	O_r  = np.zeros(shape=Ns,dtype='float')
	O_L  = np.zeros(shape=Ns,dtype='float')
	O_k  = np.zeros(shape=Ns,dtype='float')
	t    = np.loadtxt('t_a.dat')
	t    = T_h*a2tEuler(a,O_m0,O_L0,O_r0,O_k0,a_0,h,Ns)

	O_L[:] = O_L0
	O_k[:] = O_r0/a[:]**2
	O_m[:] = O_m0/a[:]**3
	O_r[:] = O_r0/a[:]**4
	a_eqrm = O_r0/O_m0
	o_eqrm = O_m0/a_eqrm**3
	a_eqmL = (O_m0/O_L0)**(1/3)
	o_eqmL = O_L0
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	fig1   = plt.figure(num=1)
	plt.clf()
	ax1    = fig1.add_subplot(111)
	ax2	   = ax1.twiny()
	ax1.semilogx(a,np.log10(O_m),'r-' ,lw=1.5,label='Density of matter'     )
	ax1.semilogx(a,np.log10(O_r),'b--',lw=1.5,label='Density of radiation'  )
	ax1.semilogx(a,np.log10(O_L),'g-.',lw=1.5,label='Density of Dark Energy')
	ax2.semilogx([t[1],t[-1]],[0,0],'g-.',lw=0)
	ax1.legend(loc='best')
	ax1.vlines(a_eqrm, np.log10(o_eqrm)-2,np.log10(o_eqrm)+2, linestyles ="dashed", colors ="k")
	ax1.vlines(a_eqmL, np.log10(o_eqmL)-2,np.log10(o_eqmL)+2, linestyles ="dashed", colors ="k")
	# ax1.text(a_eqrm,5,'Radiation-dominated',
	# 	fontproperties='Times New Roman',
	# 	fontsize=15,ha="right", va="top",
 #         bbox=dict(boxstyle="square",
 #                   ec=(0., 0., 1),
 #                   fc=(1., 1., 1.),
 #                   ))
	# ax1.text(a_eqmL,5,'Matter-dominated',
	# 	fontproperties='Times New Roman',
	# 	fontsize=15,ha="right", va="top",
 #         bbox=dict(boxstyle="square",
 #                   ec=(1, 0., 0.),
 #                   fc=(1., 1., 1.),
 #                   ))
	# a,t,z
	text0 = r'$t_{rm} \approx 4.68 \times 10^{4} yr$'
	text00 = r'$t_{m \Lambda} \approx 10.1 G yr$'
	ax1.text(2e-4,5,text0,fontsize=15)
	ax1.text(0.1,4,text00,fontsize=15)
	# r,m,L
	textr = r'$\Omega_{r,0}=8.4 \times 10^{-5}$'+'\n'
	textm = r'$\Omega_{m,0}=0.3$'+'\n'
	textL = r'$\Omega_{\Lambda,0}=0.7$'+'\n'
	textO = textL+textm+textr
	ax1.text(2e-4,-5,textO,fontsize=14)
	ax1.set_xlim(a[0],a[-1])
	ax2.set_xlim(1e9*t[1],1e9*t[-1])
	# ax2.set_xticklabels( ['{:5.2e}'.format(x) for x in np.linspace(t[1],t[-1],5,dtype='float')])
	# ax1.set_xticklabels(np.around(np.logspace(np.log10(a[1]),np.log10(a[-1]),7,dtype='float'),decimals=3))


	ax2.set_xlabel(r"$t/yr$",fontsize=15)
	ax1.minorticks_on()
	ax1.grid(b=True, which='major', color='k', linestyle=':')
	# ax1.grid(b=True, which='minor', color='k', linestyle=':', alpha=0.3)
	ax1.set_xlabel(r'$a$',fontsize=15)
	ax1.set_ylabel(r'$\lg \Omega_{I}$',fontsize=15)
	fig1.tight_layout()
	fig1.savefig('Era_rmL.png',dpi=1200)
	plt.show()
	plt.clf()



def Universes(Or0=8.4e-5,Om0=0.3,OL0=0.7,Ok0=0,H0=1/14.38,nside=2048):
	Om_1    = np.linspace(0,0.5,nside)
	OL_1    = Om_1.copy()#m-L in (0,1/2]
	Om_2    = np.linspace(1,2,nside)
	OL_2    = Om_2.copy()#m-L in [1,2]
	Om_now  = np.linspace(0,2,nside)
	OL_now  = Om_now.copy()# m+L=1
	Om_c    = Om_now.copy()
	OL_c    = Om_now.copy()# m/2=L

	OL_1[:]   = 4*Om_1[:]*np.cosh(np.arccosh(1/Om_1[:]-1)/3)**3#m-L in (0,1/2]
	OL_2[:]   = 4*Om_2[:]*np.cos((np.arccos(1/Om_2[:]-1)+4*np.pi)/3)**3#m-L in [1,2]
	OL_now[:] = 1- Om_now[:]
	OL_c[:]   = Om_c[:]/2
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	fig = plt.figure(figsize=(5,8))
	plt.clf()
	ax1 = fig.add_subplot(111)
	ax1.plot(Om_1,OL_1,'b-',lw=1.5)
	ax1.plot(Om_2,OL_2,'r-',lw=1.5)
	ax1.plot(Om_now,OL_now,'c-',lw=1.5)
	ax1.plot(Om_c,OL_c,'m-',lw=1.5)
	ax1.plot([0,1],[0,0],'k-',lw=1.5)
	# ax1.plot([1,2],[0,0],'k--',lw=1.5)
	# plt.twiny()
	plt.xticks([0,0.5,1,1.5,2])
	axx=plt.twinx()
	axy=plt.twiny()
	plt.yticks([-1,-0.5,0,0.5,1,1.5,2])
	plt.xticks([0,0.5,1,1.5,2])

	ax1.set_xlim(0,2)
	ax1.set_ylim(-1,2)
	axx.set_xlim(0,2)
	axy.set_ylim(-1,2)
	ax1.set_xlabel(r'$\Omega_{m,0}$',fontsize=12)
	ax1.set_ylabel(r'$\Omega_{\Lambda,0}$',fontsize=12)
	ax1.grid(False)


	# plt.subplots_adjust(top=1,bottom=0,left=0,right=1,hspace=0,wspace=0)
	ax1.text(0.1,0.1,r'$\Omega_{m,0}  = 2\Omega_{\Lambda,0}$',rotation=30,fontsize=12,color='m')
	ax1.text(1.3,0.6,'accelerating \n \n deccelerating',rotation=30,fontsize=15,color='m',fontproperties='Segoe Print')
	ax1.text(0.07,1.46,'no big bang \n \n big bang',rotation=63,fontsize=15,color='b',fontproperties='Segoe Print')
	ax1.text(1.3,-0.65,'closed \n \n open',rotation=-53,fontsize=18,color='c',fontproperties='Segoe Print')
	ax1.text(1.3,-0.13,'no recollapse \n \n recollapse',fontsize=15,color='r',fontproperties='Segoe Print')

	ax1.text(0.1,0.1,r'$\Omega_{m,0}  = 2\Omega_{\Lambda,0}$',rotation=30,fontsize=12,color='m')
	ax1.text(0.11,0.53,r'$\Omega_{m,0}  + \Omega_{\Lambda,0}=1$',rotation=-47,fontsize=12,color='c')
	ax1.text(0.3,-0.15,r'$\Omega_{\Lambda,0}=0$',fontsize=12,color='k')
	ax1.vlines(1,0.05,-0.05,linestyle='solid',color='k')

	ax1.text(1.6,1.3,r'$\Omega_{\Lambda,0}=4\Omega_{m,0}\cosh^{3}[ \frac{1}{3} \cosh^{-1}(\frac{1-\Omega_{m,0}}{\Omega_{m,0}})]$',
		fontproperties='Times New Roman',
		fontsize=12,ha="right", va="top",
         bbox=dict(boxstyle="square",
                   ec=(0, 0., 1.),
                   fc=(1., 1., 1.),
                   ))
	ax1.arrow(0.01,-0.05, 1-0.01,0,
              length_includes_head=True,
              head_width=0.05, head_length=0.04,#shape="full",
              fc='k', ec='k',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.arrow(1,-0.05, 0.01-1,0,
              length_includes_head=True,
              head_width=0.05, head_length=0.04,#shape="full",
              fc='k', ec='k',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.arrow(1.1,-0.05, 0.9,0.05,
              length_includes_head=True,
              head_width=0.05, head_length=0.04,#shape="full",
              fc='k', ec='k',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.arrow(2,0, -0.9,-0.05,
              length_includes_head=True,
              head_width=0.05, head_length=0.04,#shape="full",
              fc='k', ec='k',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.arrow(1.6,1.2, Om_1[300]-1.6,OL_1[300]-1.2,
              length_includes_head=True,
              head_width=0., head_length=0.,#shape="full",
              fc='b', ec='b',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.arrow(0.8,-0.7, Om_2[1000]-1.0,OL_2[1000]+0.69,
              length_includes_head=True,
              head_width=0., head_length=0.,#shape="full",
              fc='r', ec='r',overhang=0.5,alpha=0.9,
              # color=c1
              )
	ax1.text(1.5,-0.7,r'$\Omega_{\Lambda,0}=4\Omega_{m,0}\cos^{3}[ \frac{1}{3} \cos^{-1}(\frac{1-\Omega_{m,0}}{\Omega_{m,0}})+\frac{4\pi}{3}]$',
		fontproperties='Times New Roman',
		fontsize=12,ha="right", va="top",
         bbox=dict(boxstyle="square",
                   ec=(1, 0., 0.),
                   fc=(1., 1., 1.),
                   ))
	ax1.text(0.8,1.7,r'$\Omega_{k,0}=1-\Omega_{\Lambda,0}-\Omega_{m,0}$',fontsize=15)
	ax1.set_title('FRLW Universes',fontproperties='Segoe Print',fontsize=15)

	fig.savefig('Universes.png',dpi=1200)
	plt.show()
	plt.clf()


def cosmo_time(a_t=1,O_m0=0.3,O_L0=0.7,O_r0=8.4e-5,O_k0=0,a_0=1,h=0.68,nside=2048):
	# z0 = a_0/a-1
	t = 0
	a = np.linspace(1e-10,a_t,nside)
	for i in range(nside):
		t +=  (a[i]-a[i-1])/a[i]/(O_L0+O_k0/a[i]**2+O_m0/a[i]**3+O_r0/a[i]**4)**(1/2)
	return t*9.78/h





if __name__ == '__main__':
	main()