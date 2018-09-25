import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf=matplotlib.backends.backend_pdf.PdfPages("Lcer.pdf")

font1={'family':'Impact','weight':'normal','size':70}
data=pd.read_table("OUT/Lcer_ped.dat", delim_whitespace=True)
# plot ratio
fig,axes=plt.subplots(1,1,figsize=(30,30))
color=['k','r','g','b']
axes.errorbar(data['nrun'],data['cer1'],linestyle="",marker="^",color="green",label="LHRS cer1",ms=30)

axes.legend(loc='best',ncol=2,prop=font1)
#plt.ylim(0.7,1.1)
axes.set_ylabel('Pedestal peak',font1)
axes.set_xlabel('run_number',font1)
axes.tick_params(labelsize=60)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer2'],linestyle="",marker="^",color="green",label="LHRS cer2",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer3'],linestyle="",marker="^",color="green",label="LHRS cer3",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer4'],linestyle="",marker="^",color="green",label="LHRS cer4",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer5'],linestyle="",marker="^",color="green",label="LHRS cer5",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer6'],linestyle="",marker="^",color="green",label="LHRS cer6",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer7'],linestyle="",marker="^",color="green",label="LHRS cer7",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer8'],linestyle="",marker="^",color="green",label="LHRS cer8",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer9'],linestyle="",marker="^",color="green",label="LHRS cer9",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['cer10'],linestyle="",marker="^",color="green",label="LHRS cer10",ms=30)
pdf.savefig(fig)

pdf.close()
