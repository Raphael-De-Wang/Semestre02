#!env python

import numpy as np
import operator

# from function import loadToDict, load_pfam2go_toDict, GOAnnotation, plot_switch_view, switch_view, domain_function_list, domain_function_dict
from function import plot_switch_view, switch_view, domain_function_list, GOAnnotation, load_all_func_data

pfam2go_dict, goslimmeta_dict, goDict, dom_gcai, dom_abundance = load_all_func_data()
    
da = sorted(dom_abundance.items(), key=operator.itemgetter(1),reverse=True)
dom_list = np.array(da)[:,0]
dfList = domain_function_list(dom_list,pfam2go_dict)

bio_process,molec_func,cellu_comp = switch_view(dfList,goDict,goslimmeta_dict,dom_gcai,dom_abundance)

plot_switch_view('Biological Process',bio_process,'BiologicalProcess.png',(35,80),caiSeuil=0.5)
plot_switch_view('Cellular Component',cellu_comp, 'CellularComponent.png',(35,30),caiSeuil=0.5)
plot_switch_view('Molecular Function',molec_func, 'MolecularFunction.png',(35,90),caiSeuil=0.5)

