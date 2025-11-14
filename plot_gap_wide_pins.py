import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os

ap = ArgumentParser()
ap.add_argument("--N_el", "-e", type=int, default=6)
ap.add_argument("--N_orb", "-o", type=int, default=6)
ap.add_argument("--pin_size","-k",type=int,nargs="+",default=[2])
ap.add_argument("--Lz",type=int)
ap.add_argument("--asymmetric",action="store_true",default=False)
ap.add_argument("--special", action="store_true",default=False)
ap.add_argument("--both_positive", action="store_true",default=False)
ap.add_argument("--use_krylov", action="store_true",default=False)
ap.add_argument("--equator", action="store_true",default=False)
ap.add_argument("--interaction", "-i",type=str, default="v1.txt")
ap.add_argument("--north_only", action="store_true",default=False)
ap.add_argument("--positive_north",action="store_true",default=False)
ap.add_argument("--scale_lambda",action="store_true",default=False)
ap.add_argument("--save_plot_data",action="store_true",default=False)
ap.add_argument("--normalize_pins",action="store_true",default=False)
ap.add_argument("--Lz_cutoff",type=float,default=5)
ap.add_argument("--n_gaps","-n",type=int,default=1)

aa = ap.parse_args()
Ne = aa.N_el
No = aa.N_orb

markers_available = ["x--","o--","v--","^--","s--", "d--"]
markers_list      = (markers_available*(aa.n_gaps//len(markers_available)+1))[:aa.n_gaps]

round1dp = lambda x: round(10*x)/10.
if abs(aa.Lz_cutoff-round(aa.Lz_cutoff)) < 0.1: # if integer L
	cutoff = int(aa.Lz_cutoff) 
else: # if half-integer L
	cutoff = round1dp(aa.Lz_cutoff)

fig,ax = plt.subplots(len(aa.pin_size),1,figsize=(5,3*len(aa.pin_size)))


if aa.use_krylov:
	extension = "_krylov"
else:
	extension = ""

if aa.north_only:
	if aa.positive_north:
		northtext = "_positive_north"
	else:
		northtext = "_north"
else:
	northtext = ""

def get_spectrum(k,lamb):
	if aa.both_positive:
		if aa.equator:
			dir_name = "equator_pins_positive_samesign"
		else:
			dir_name = f"wide_pin_k_{k}_positive_samesign"
	else:
		if aa.special:
			if aa.asymmetric:
				dir_name = f"asymmetric_wide_pin_k_{k}_special"
			else:
				dir_name = f"wide_pin_k_{k}_special{northtext}"
		else:
			if aa.asymmetric:
				dir_name = f"asymmetric_wide_pin_k_{k}"
			else:
				dir_name = f"wide_pin_k_{k}{northtext}"
	E  = np.array([])
	LZ = np.array([])
	#for Lz in np.arange(-Ne,Ne+1):
	if aa.north_only:
		if aa.positive_north:
			Lzrange = np.arange(aa.Lz_cutoff+1,1)
		else:
			Lzrange = np.arange(-aa.Lz_cutoff,1,1)
	else:
		Lzrange = np.arange(-aa.Lz_cutoff,aa.Lz_cutoff+1,1)
	for Lz_val in Lzrange:
		if abs(Lz_val - round(Lz_val)) < 0.1:
			Lz = round(Lz_val)
		else:
			Lz = round1dp(Lz_val)
		try:
			with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.3f}_out/eigen{extension}_L_{Lz}.txt") as f:
				data = [list(map(float,line.split())) for line in f.readlines()]
				E    = np.concatenate((E,np.array([datum[0] for datum in data])))
				LZ   = np.concatenate((LZ,np.array([datum[1] for datum in data])))
		except FileNotFoundError:
			try:
				with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.2f}_out/eigen{extension}_L_{Lz}.txt") as f:
					data = [list(map(float,line.split())) for line in f.readlines()]
					E    = np.concatenate((E,np.array([datum[0] for datum in data])))
					LZ   = np.concatenate((LZ,np.array([datum[1] for datum in data])))
			except FileNotFoundError:
				try:
					with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.1f}_out/eigen{extension}_L_{Lz}.txt")as f:
						data = [list(map(float,line.split())) for line in f.readlines()]
						E    = np.concatenate((E,np.array([datum[0] for datum in data])))
						LZ   = np.concatenate((LZ,np.array([datum[1] for datum in data])))
				except FileNotFoundError:
					print(f"Data not found for lambda = {lamb}")
					print("Tried:")
					print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.3f}_out/eigen{extension}_L_{Lz}.txt")
					print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.2f}_out/eigen{extension}_L_{Lz}.txt")
					print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.1f}_out/eigen{extension}_L_{Lz}.txt")
					continue
	E_ret  = np.array(sorted(E))
	LZ_ret = np.array([x for _,x in sorted (zip(E,LZ))])
	return E_ret, LZ_ret


if aa.north_only:
	lambda_range = np.concatenate((np.arange(0,1,0.05),np.arange(1,5.5,0.5),[10,20,50,100,1000]))
else:
	lambda_range = np.arange(0,20.1,1)

def save_plot_data(x,y,filename):
	with open(f"plot_data/{filename}.txt","w+") as f:
		for (xx,yy) in zip(x,y):
			f.write(f"{xx}\t{yy}\n")



if aa.scale_lambda:
	filesuffix = "_scale"
else:
	filesuffix = ""
if aa.special:
	fname_append = f"{northtext}_special"
else:
	fname_append = northtext
	
if aa.both_positive:
	if aa.equator:
		if aa.Lz != None:
			savename = f"spectrum_gap_equator_pins_positive_samesign_{Ne}e_{No}o_Lz_{aa.Lz}{extension}{filesuffix}"
		else:
			savename = f"spectrum_gap_equator_pins_positive_samesign_{Ne}e_{No}o{extension}{filesuffix}"
	else:
		if aa.Lz != None:
			savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o_Lz_{aa.Lz}{fname_append}_positive_samesign{extension}{filesuffix}"
		else:
			savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o{fname_append}_positive_samesign{extension}{filesuffix}"
elif aa.asymmetric:
	if aa.Lz != None:
		savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o_Lz_{aa.Lz}{fname_append}_asym{extension}{filesuffix}"
	else:
		savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o{fname_append}_asym{extension}{filesuffix}"
else:
	if aa.Lz != None:
		savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o_Lz_{aa.Lz}{fname_append}{extension}{filesuffix}"
	else:
		savename = f"spectrum_gap_k_{aa.pin_size}_{Ne}e_{No}o{fname_append}{extension}{filesuffix}"

for ki,k in enumerate(aa.pin_size):
	lambda_plot = []
	gap = [[] for _ in range(aa.n_gaps)]
	gsL = []
	for lamb in lambda_range:
		E, LZ = get_spectrum(k,lamb)
		if len(E)>0:
			if aa.scale_lambda:
				lambda_plot.append(1/lamb)
			else:
				lambda_plot.append(lamb)

			if aa.Lz != None:
				E_6 = E[np.abs(LZ+aa.Lz)<1e-5]
				gap.append(E_6[1]-E_6[0])
				gsL.append(LZ[0])
			else:
				for ni in range(aa.n_gaps):
					gap[ni].append(E[ni+1]-E[0])
				gsL.append(LZ[0])
		print(f"k={k}\tλ={lamb:.2f}\t\tgap = {E[1]-E[0]:.5f}", end="\n\t\t")
		print("\n\t\t".join(f"{E[i]:.5f}   {LZ[i]:.1f}" for i in range(3)))
	for ni in range(aa.n_gaps):
		ax[ki].plot(lambda_plot, gap[ni],markers_list[ni],markerfacecolor="none",label=r"$\Delta E$_"+str(ni+1))
	if aa.save_plot_data:
		for ni in range(aa.n_gaps):
			save_plot_data(lambda_plot,gap[ni],savename.replace(f"{aa.pin_size}",f"{k}")+f"_E_{ni+1}")


	if aa.scale_lambda:
		ax[ki].set_xlabel(r"1/$\lambda$")
	else:
		ax[ki].set_xlabel(r"$\lambda$")

	ax[ki].set_ylabel(r"$\Delta E$_i = E_i - E_0")
	ax[ki].set_xlim([0,20])
	ax[ki].legend()
	ax[ki].set_title(f"k={k}")
#plt.suptitle(r"$H=V_1+V_0$")
if aa.both_positive:
	if aa.equator:
		plt.suptitle(r"$H=V_1 + \lambda (V_{eq1}+V_{eq2})$"+f", Ne={Ne}, No={No}")
	else:
		plt.suptitle(r"$H=V_1 + \lambda (V_{k-north}+V_{k-south})$"+f", Ne={Ne}, No={No}")	

else:
	if aa.north_only:
		if aa.positive_north:
			plt.suptitle(r"$H=V_1 + \lambda V_{k-north}$"+f", Ne={Ne}, No={No}")
		else:
			plt.suptitle(r"$H(1-\lambda)V_1 - \lambda V_{k-north}$"+f", Ne={Ne}, No={No}")
	else:
		plt.suptitle(r"$H(1-\lambda)V_1 + \lambda (-V_{k-north}+V_{k-south})$"+f", Ne={Ne}, No={No}")
plt.tight_layout()
plt.legend()


plt.savefig(f"{savename}.pdf")