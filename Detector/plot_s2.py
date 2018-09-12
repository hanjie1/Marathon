import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf=matplotlib.backends.backend_pdf.PdfPages("LS2L.pdf")

font1={'family':'Impact','weight':'normal','size':70}
data=pd.read_table("OUT/Ls2L_ped.dat", delim_whitespace=True)
# plot ratio
fig,axes=plt.subplots(1,1,figsize=(30,30))
color=['k','r','g','b']
axes.errorbar(data['nrun'],data['S2L1'],linestyle="",marker="^",color="green",label="LHRS S2L1",ms=30)

axes.legend(loc='best',ncol=2,prop=font1)
#plt.ylim(0.7,1.1)
axes.set_ylabel('single photon peak',font1)
axes.set_xlabel('run_number',font1)
axes.tick_params(labelsize=60)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L2'],linestyle="",marker="^",color="green",label="LHRS S2L2",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L3'],linestyle="",marker="^",color="green",label="LHRS S2L3",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L4'],linestyle="",marker="^",color="green",label="LHRS S2L4",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L5'],linestyle="",marker="^",color="green",label="LHRS S2L5",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L6'],linestyle="",marker="^",color="green",label="LHRS S2L6",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L7'],linestyle="",marker="^",color="green",label="LHRS S2L7",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L8'],linestyle="",marker="^",color="green",label="LHRS S2L8",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L9'],linestyle="",marker="^",color="green",label="LHRS S2L9",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L10'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)

axes.errorbar(data['nrun'],data['S2L11'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L12'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L13'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L14'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L15'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)
axes.errorbar(data['nrun'],data['S2L16'],linestyle="",marker="^",color="green",label="LHRS S2L10",ms=30)
pdf.savefig(fig)


pdf.close()
