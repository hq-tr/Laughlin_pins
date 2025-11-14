import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os

ap = ArgumentParser()
ap.add_argument("--N_el", "-e", type=int, default=6)
ap.add_argument("--N_orb", "-o", type=int, default=6)
ap.add_argument("--pin_size","-k",type=int,default=1)
ap.add_argument("--positive","-p",action="store_true",default=False)
ap.add_argument("--use_krylov", action="store_true",default=False)
ap.add_argument("--interaction", "-i",type=str, default="")
ap.add_argument("--scale_lambda",action="store_true",default=False)
ap.add_argument("--n_gaps","-n",type=int,default=1)
ap.add_argument("--save_plot_data",action="store_true",default=False)
aa = ap.parse_args()
Ne = aa.N_el
No = aa.N_orb

fig = plt.subplots(figsize=(6,4))

markers_available = ["x--","o--","v--","^--","s--", "d--"]
markers_list      = (markers_available*(aa.n_gaps//len(markers_available)+1))[:aa.n_gaps]


if aa.use_krylov:
	extension = "_krylov"
else:
	extension = ""

def get_spectrum(k,lamb):
	if aa.positive:
		dir_name = f"wide_pin_k_{k}_4pins_positive"
	else:
		dir_name = f"wide_pin_k_{k}_4pins"
	try:
		with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.3f}_out/eigen{extension}.txt") as f:
			data = [list(map(float,line.split())) for line in f.readlines()]
			E    = np.array([datum[0] for datum in data])
			LZ   = np.array([datum[1] for datum in data])
	except FileNotFoundError:
		try:
			with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.2f}_out/eigen{extension}.txt") as f:
				data = [list(map(float,line.split())) for line in f.readlines()]
				E    = np.array([datum[0] for datum in data])
				LZ   = np.array([datum[1] for datum in data])
		except FileNotFoundError:
			try:
				with open(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.1f}_out/eigen{extension}.txt")as f:
					data = [list(map(float,line.split())) for line in f.readlines()]
					E    = np.array([datum[0] for datum in data])
					LZ   = np.array([datum[1] for datum in data])
			except FileNotFoundError:
				print(f"Data not found for lambda = {lamb}")
				print("Tried:")
				print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.3f}_out/eigen{extension}.txt")
				print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.2f}_out/eigen{extension}.txt")
				print(f"{dir_name}/{Ne}e_{No}o_{aa.interaction}_{lamb:.1f}_out/eigen{extension}.txt")
				E  = np.array([])
				LZ = np.array([])
	return E, LZ

def save_plot_data(x,y,filename):
	with open(f"plot_data/{filename}.txt","w+") as f:
		for (xx,yy) in zip(x,y):
			f.write(f"{xx}\t{yy}\n")

if aa.scale_lambda:
	filesuffix = "_scale"
else:
	filesuffix = ""

savename = f"spectrum_gap_4pins_k_{aa.pin_size}_{Ne}e_{No}o{extension}{filesuffix}"

lambda_range = np.arange(0,20.1,1)
lambda_range = np.arange(0,1.6,0.1)

k = aa.pin_size
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

		for ni in range(aa.n_gaps):
			gap[ni].append(E[ni+1]-E[0])
		gsL.append(LZ[0])

	print(f"k={k}\tλ={lamb:.2f}\t\tgap = {E[1]-E[0]:.5f}", end="\n\t\t")
	print("\n\t\t".join(f"{E[i]:.5f}   {LZ[i]:.1f}" for i in range(3)))
for ni in range(aa.n_gaps):
	plt.plot(lambda_plot, gap[ni],markers_list[ni],markerfacecolor="none",label=r"$\Delta E$_"+str(ni+1))

if aa.save_plot_data:
	for ni in range(aa.n_gaps):
		save_plot_data(lambda_plot,gap[ni],savename + f"_E_{ni+1}")


if aa.scale_lambda:
	plt.xlabel(r"1/$\lambda$")
else:
	plt.xlabel(r"$\lambda$")

plt.ylabel(r"$\Delta E$_i = E_i - E_0")
#plt.xlim([0,20])
plt.legend()

plt.title(f"V_1 + λV_4pins (k={k}) (Ne={Ne},No={No}) ")
plt.tight_layout()
plt.legend()


plt.savefig(f"{savename}.pdf")