
def heatmap(eforge_dataset,savefig,state=None,width,gamma=1.0):
    """This function was used to make figure 1

    Parameters
    -------------------------------------------------------------
    eforge_dataset: string
        dataset used e.g. encode, erc2-DHS, erc.chart, erc2-H3...
    savefig: string
        path/to/folder/file.pdf
    state: string
        chomatin state and histone mark result .tsv files are aggregated, use this to only take a specific value e.g. 'Enh'
    width: int
        width of figure
    gamma: int
        spread of color to highlight more enriched tissues, higher values increases white spread over blue


    Return
    -------------------------------------------------------------
    Seaborn clustermap of eFORGE analyses performed with the unconsolidated DNase I hotspots data from Epigenomics Roadmap consortium.
    """
    import os
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib import cm
    import matplotlib.gridspec
    from matplotlib.ticker import ScalarFormatter, FuncFormatter
    from matplotlib.transforms import Affine2D
    from matplotlib.collections import PathCollection

    #Folder and file paths
    folder_path="0_tsv_files/"
    sample_size_file=f'{folder_path}0_meta_data.xlsx'

    #Read q values
    files = os.listdir(folder_path)
    column_data_list= list()
    if eforge_dataset=='erc.chart':
        cell_types=['CD3','CD3','CD3','CD3','CD4_Mobilised;','CD4_Mobilised;','CD4','CD4','CD4','CD4','CD8_Mobilised;','CD8','CD8','CD8','CD8','CD14','CD14','CD19','CD19','CD19','CD20','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD34_Mobilised;','CD56_Mobilised;','CD56','CD56','Breast_vHMEC','Breast_vHMEC','H1','H1','H9','H9','Mesenchymal_Stem_Cells','Mesenchymal_Stem_Cells','Mesendoderm_Cultured_Cells','Mesendoderm_Cultured_Cells','Neuronal_Progenitor_Cultured_Cells','Neuronal_Progenitor_Cultured_Cells','Trophoblast_Cultured_Cells','Trophoblast_Cultured_Cells','Fetal_Adrenal_Gland','Fetal_Adrenal_Gland','Fetal_Adrenal_Gland','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Brain','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Heart','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Large','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Intestine_Small','Fetal_Kidney_Left','Fetal_Kidney_Left','Fetal_Kidney_Left','Fetal_Kidney_Left','Fetal_Kidney_Left','Fetal_Kidney_Right','Fetal_Kidney_Right','Fetal_Kidney_Right','Fetal_Kidney_Right','Fetal_Kidney_Right','Fetal_Kidney','Fetal_Kidney','Fetal_Kidney','Fetal_Kidney','Fetal_Kidney','Fetal_Kidney','Fetal_Kidney','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Left','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung_Right','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Lung','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Arm','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Back','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Leg','Fetal_Muscle_Lower_Limb_Skeletal','Fetal_Muscle_Trunk','Fetal_Muscle_Trunk','Fetal_Muscle_Trunk','Fetal_Muscle_Upper_Limb_Skeletal','Fetal_Muscle_Upper_Trunk','Fetal_Placenta','Fetal_Placenta','Fetal_Placenta','Fetal_Placenta','Fetal_Renal_Cortex_Left','Fetal_Renal_Cortex_Left','Fetal_Renal_Cortex_Left','Fetal_Renal_Cortex_Left','Fetal_Renal_Cortex_Right','Fetal_Renal_Cortex_Right','Fetal_Renal_Cortex_Right','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Cortex','Fetal_Renal_Pelvis_Left','Fetal_Renal_Pelvis_Left','Fetal_Renal_Pelvis_Left','Fetal_Renal_Pelvis_Left','Fetal_Renal_Pelvis_Right','Fetal_Renal_Pelvis_Right','Fetal_Renal_Pelvis_Right','Fetal_Renal_Pelvis_Right','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Renal_Pelvis','Fetal_Skin','Fetal_Spinal_Cord','Fetal_Spinal_Cord','Fetal_Spinal_Cord','Fetal_Spleen','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Stomach','Fetal_Testes','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fetal_Thymus','Fibroblasts_Fetal_Skin_Abdomen','Fibroblasts_Fetal_Skin_Abdomen','Fibroblasts_Fetal_Skin_Back','Fibroblasts_Fetal_Skin_Biceps_Left','Fibroblasts_Fetal_Skin_Biceps_Left','Fibroblasts_Fetal_Skin_Biceps_Right','Fibroblasts_Fetal_Skin_Biceps_Right','Fibroblasts_Fetal_Skin_Quadriceps_Left','Fibroblasts_Fetal_Skin_Quadriceps_Left','Fibroblasts_Fetal_Skin_Quadriceps_Right','Fibroblasts_Fetal_Skin_Quadriceps_Right','Fibroblasts_Fetal_Skin_Scalp','Fibroblasts_Fetal_Skin_Scalp','Fibroblasts_Fetal_Skin_Upper_Back','Fibroblasts_Fetal_Skin_Upper_Back','Penis_Foreskin_Fibroblast','Penis_Foreskin_Fibroblast','Penis_Foreskin_Fibroblast','Penis_Foreskin_Fibroblast','iPS_DF_4.7','iPS_DF_6.9','iPS_DF_19.7','iPS_DF_19.11','IMR90','IMR90','IMR90','IMR90','Penis_Foreskin_Keratinocyte','Penis_Foreskin_Keratinocyte','Penis_Foreskin_Keratinocyte','Penis_Foreskin_Keratinocyte','Penis_Foreskin_Melanocyte','Penis_Foreskin_Melanocyte']
    elif 'chromatin' in eforge_dataset:
        cell_types=['E063 Adipose Nuclei','E080 Fetal Adrenal Gland','E029 Primary monocytes from peripheral blood','E030 Primary neutrophils from peripheral blood','E031 Primary B cells from cord blood','E032 Primary B cells from peripheral blood','E033 Primary T cells from cord blood','E034 Primary T cells from peripheral blood','E035 Primary hematopoietic stem cells','E036 Primary hematopoietic stem cells short term culture','E037 Primary T helper memory cells from peripheral blood 2','E038 Primary T helper naive cells from peripheral blood','E039 Primary T helper naive cells from peripheral blood','E040 Primary T helper memory cells from peripheral blood 1','E041 Primary T helper cells PMA-I stimulated','E042 Primary T helper 17 cells PMA-I stimulated','E043 Primary T helper cells from peripheral blood','E044 Primary T regulatory cells from peripheral blood','E045 Primary T cells effector/memory enriched from peripheral blood','E046 Primary Natural Killer cells from peripheral blood','E047 Primary T CD8+ naive cells from peripheral blood','E048 Primary T CD8+ memory cells from peripheral blood','E050 Primary hematopoietic stem cells G-CSF-mobilized Female','E051 Primary hematopoietic stem cells G-CSF-mobilized Male','E062 Primary mononuclear cells from peripheral blood','E115 Dnd41 TCell Leukemia','E116 GM12878 Lymphoblastoid','E123 K562 Leukemia','E124 Monocytes-CD14+ RO01746 Primary Cells','E129 Osteoblast Primary Cells','E067 Brain Angular Gyrus','E068 Brain Anterior Caudate','E069 Brain Cingulate Gyrus','E070 Brain Germinal Matrix','E071 Brain Hippocampus Middle','E072 Brain Inferior Temporal Lobe','E073 Brain_Dorsolateral_Prefrontal_Cortex','E074 Brain Substantia Nigra','E081 Fetal Brain Male','E082 Fetal Brain Female','E125 NH-A Astrocytes Primary Cells','E119 HMEC Mammary Epithelial Primary Cells','E117 HeLa-S3 Cervical Carcinoma','E075 Colonic Mucosa','E077 Duodenum Mucosa','E079 Esophagus','E084 Fetal Intestine Large','E085 Fetal Intestine Small','E092 Fetal Stomach','E094 Gastric','E101 Rectal Mucosa Donor 29','E102 Rectal Mucosa Donor 31','E106 Sigmoid Colon','E109 Small Intestine','E110 Stomach Mucosa','E122 HUVEC Umbilical Vein Endothelial Primary Cells','E027 Breast Myoepithelial Primary Cells','E028 Breast variant Human Mammary Epithelial Cells (vHMEC)','E055 Foreskin Fibroblast Primary Cells skin01','E056 Foreskin Fibroblast Primary Cells skin02','E057 Foreskin Keratinocyte Primary Cells skin02','E058 Foreskin Keratinocyte Primary Cells skin03','E059 Foreskin Melanocyte Primary Cells skin01','E061 Foreskin Melanocyte Primary Cells skin03','E126 NHDF-Ad Adult Dermal Fibroblast Primary Cells','E127 NHEK-Epidermal Keratinocyte Primary Cells','E001 ES-I3','E002 ES-WA7','E003 H1','E008 H9','E014 HUES48','E015 HUES6','E016 HUES64','E024 ES-UCSF4','E004 H1 BMP4 Derived Mesendoderm Cultured','E005 H1 BMP4 Derived Trophoblast Cultured','E006 H1 Derived Mesenchymal Stem Cells','E007 H1 Derived Neuronal Progenitor Cultured','E009 H9 Derived Neuronal Progenitor Cultured','E010 H9 Derived Neuron Cultured','E011 hESC Derived CD184+ Endoderm Cultured','E012 hESC Derived CD56+ Ectoderm Cultured','E013 hESC Derived CD56+ Mesoderm Cultured','E065 Aorta','E083 Fetal Heart','E095 Left Ventricle','E104 Right Atrium','E105 Right Ventricle','E018 iPS-15b','E019 iPS-18','E020 iPS-20b','E021 iPS DF 6.9','E022 iPS DF 19.11','E086 Fetal Kidney','E066 Liver','E118 HepG2 Hepatocellular Carcinoma','E017 IMR90 fetal lung fibroblasts','E088 Fetal Lung','E096 Lung','E114 A549 EtOH 0.02pct Lung Carcinoma','E128 NHLF Lung Fibroblast Primary Cells','E023 Mesenchymal Stem Cell Derived Adipocyte Cultured','E025 Adipose Derived Mesenchymal Stem Cell Cultured','E026 Bone Marrow Derived Cultured Mesenchymal Stem Cells','E049 Mesenchymal Stem Cell Derived Chondrocyte Cultured','E052 Muscle Satellite Cultured','E089 Fetal Muscle Trunk','E090 Fetal Muscle Leg','E100 Psoas Muscle','E107 Skeletal Muscle Male','E108 Skeletal Muscle Female','E120 HSMM Skeletal Muscle Myoblasts','E121 HSMM cell derived Skeletal Muscle Myotubes','E053 Cortex derived primary cultured neurospheres','E054 Ganglion Eminence derived primary cultured neurospheres','E097 Ovary','E087 Pancreatic Islets','E098 Pancreas','E091 Placenta','E099 Placenta Amnion','E076 Colon Smooth Muscle','E078 Duodenum Smooth Muscle','E103 Rectal Smooth Muscle','E111 Stomach Smooth Muscle','E113 Spleen','E093 Fetal Thymus','E112 Thymus']
    elif 'erc2'in eforge_dataset:
        cell_types=['E029 Primary monocytes from peripheral blood','E032 Primary B cells from peripheral blood','E033 Primary T cells from cord blood','E034 Primary T cells from peripheral blood','E046 Primary Natural Killer cells from peripheral ','E050 Primary hematopoietic stem cells G-CSF-mobili','E051 Primary hematopoietic stem cells G-CSF-mobili','E028 Breast variant Human Mammary Epithelial Cells','E003 H1 Cells','E004 H1 BMP4 Derived Mesendoderm Cultured Cells','E005 H1 BMP4 Derived Trophoblast Cultured Cells','E006 H1 Derived Mesenchymal Stem Cells','E007 H1 Derived Neuronal Progenitor Cultured Cells','E008 H9 Cells','E085 Fetal Intestine Small','E080 Fetal Adrenal Gland','E081 Fetal Brain Male','E082 Fetal Brain Female','E083 Fetal Heart','E084 Fetal Intestine Large','E086 Fetal Kidney','E088 Fetal Lung','E090 Fetal Muscle Leg','E089 Fetal Muscle Trunk','E092 Fetal Stomach','E093 Fetal Thymus','E094 Gastric','E021 iPS DF 6.9 Cells','E022 iPS DF 19.11 Cells','E017 IMR90 fetal lung fibroblasts Cell Line','E097 Ovary','E098 Pancreas','E091 Placenta','E100 Psoas Muscle','E055 Foreskin Fibroblast Primary Cells skin01','E056 Foreskin Fibroblast Primary Cells skin02','E057 Foreskin Keratinocyte Primary Cells skin02','E059 Foreskin Melanocyte Primary Cells skin01','E109 Small Intestine']
    elif 'encode' in eforge_dataset:
        cell_types=['HTR8svn','Adult_CD4+','CD14+','CD20+','CD34+','CLL','CMK','GM06990','GM12864','GM12865','GM12878','GM12891','GM12892','GM18507','GM19238','GM19239','GM19240','HL-60','Jurkat','K562','NB4','Th1','Th2','AoAF','AoSMC','HBMEC','HMVEC-dAd','HMVEC-dBl-Ad','HMVEC-dBl-Neo','HMVEC-dLy-Ad','HMVEC-dLy-Neo','HMVEC-dNeo','HMVEC-LBl','HMVEC-LLy','HPAEC','HPAF','HUVEC','Osteobl','BE2_C','Gliobla','Medullo','SK-N-MC','SK-N-SH','HA-h','HMEC','HMF','MCF-7','MCF-7','T-47D','HAc','HeLa-S3','HeLa-S3','Caco-2','HCT-116','HVMF','WI-38','WI-38','A549','HAEpiC','HCPEpiC','HEEpiC','HIPEpiC','HNPCEpiC','HPdLF','HRCEpiC','HRE','HRPEpiC','pHTE','RPTEC','SAEC','HESC','hESCT0','H9ES','HConF','WERI-Rb-1','Chorion','HFF','HFF-Myc','AG09319','HGF','HCFaa','HCF','HCM','iPS','HRGEC','8988T','Hepatocytes','HepG2','Huh-7','Huh-7.5','Stellate','AG04450','HPF','NHLF','E_myoblast','HSMM','HSMM','SKMC','Myometr','NH-A','PANC-1','PanIsletD','PanIslets','HPDE6-E6E7','LNCaP','LNCaP','PrEC','RWPE1','AG04449','AG09309','AG10803','BJ','Fibrobl','FibroP','Melano','NHDF-Ad','NHDF-neo','NHEK','ProgFib','HA-sp','NT2-D1','Urothelia','Urothelia','Ishikawa','Ishikawa']
    elif 'blueprint' in eforge_dataset:
        cell_types=['Acute_myeloid_leukemia','DG-75_Sporadic_Burkitt_lymphoma','KARPAS-422_Germinal_Center_B-Cell-Like_Diffuse_Large_B-Cell_Lymphoma','SU-DHL-5_Germinal_Center_B-Cell-Like_Diffuse_Large_B-Cell_Lymphoma','U-266_Multiple_myeloma','Z-138_Mantle_cell_lymphoma','CD14_positive_CD16_negative_classical_monocyte','CD34-negative_CD41-positive_CD42-positive_megakaryocyte_cell','CD14-positive_CD16-negative_classical_monocyte','CD14-positive_CD16-negative_classical_monocyte','CD14-positive_CD16-negative_classical_monocyte','CD14-positive_CD16-negative_classical_monocyte','inflammatory_macrophage','inflammatory_macrophage','macrophage_-_T_6days_B-glucan','macrophage_-_T_6days_B-glucan','macrophage_-_T_6days_B-glucan','macrophage_-_T_6days_B-glucan','macrophage_-_T_6days_LPS','macrophage_-_T_6days_LPS','macrophage_-_T_6days_LPS','macrophage_-_T_6days_LPS','macrophage_-_T_6days_untreated','macrophage_-_T_6days_untreated','macrophage_-_T_6days_untreated','macrophage_-_T_6days_untreated','macrophage','macrophage','monocyte_-_T_0days','monocyte_-_T_0days']
    
    column_name = "Qvalue"
    mask_list=list()
    included_files=list()
    for file_name in files:
        if eforge_dataset in file_name:
            file_path = os.path.join(folder_path, file_name)
            
            df = pd.read_csv(file_path, sep="\t")
            if state != None:
                df_filtered=df[df['Datatype']==state]
            else:
                df_filtered=df
            if column_name in df_filtered.columns:
                column_data = df_filtered[column_name].tolist()
                
                mask= [x > 0.05 for x in column_data]
                
                if False in mask:
                    column_data_list.append(column_data)
                    included_files.append(file_name.split('.txt')[0])
                    mask_list.append(mask)
                
            else:
                print(f"Column '{column_name}' not found in file '{file_name}'.")
    
    #Read meta data for corresponding files
    df=pd.read_excel(sample_size_file)
    
    filtered_df = df[df['study'].isin(included_files)]
    sample_list=filtered_df['sample']
    sample2_list=filtered_df['sample2']
    
    cmap = mcolors.LinearSegmentedColormap.from_list("", list(zip([0,1],[(1, 1, 1), (0, 0, 0.5)])),6,gamma=1.5)


    sample_dict = {'Whole blood': 'red', 'Cord blood': [1, 0.5, 0.5], 'Peripheral Blood Leukocyte': 'orange', 'Peripheral Blood Mononuclear Cells (PBMCs)':'orange','Buccal':'green','CD4+ T cells': 'orange','/':'grey', 'Placenta': 'yellow','Lymphoblast':'orange','Lung alveolarmacrophage':'purple'}
    sample_colors= sample_list.map(sample_dict)
    sample2_colors=sample2_list.map(sample_dict)

    log_data=-np.log10(column_data_list)
    
    clustergrid = sns.clustermap(log_data,row_colors=[sample_colors,sample2_colors], metric='euclidean',z_score=0,
                                dendrogram_ratio=0.05,row_cluster=True, col_cluster=False,cmap=cmap, yticklabels=included_files, xticklabels=cell_types, 
                                method='average', figsize=(width,8),linewidths=0.003,linecolor='grey',mask=np.array(mask_list))
    
    clustergrid.tick_params(axis='x', labelsize=4,labeltop=True, top=True, labelbottom=False, bottom=False)
    clustergrid.tick_params(axis='y', labelsize=6)
    plt.setp(clustergrid.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(clustergrid.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    clustergrid.gs.update(left=0.08,right=0.7,bottom=0.05,top=0.75)

    #Move clustermap subplot field and add lineplot
    gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.03, right=0.07)
    ax2 = clustergrid.fig.add_subplot(gs2[0])

    #Function to rotate lineplot
    def lineplot_plusplus(orientation = "horizontal", **kwargs):
        line = sns.lineplot(**kwargs)

        r = Affine2D().scale(sx=1, sy=-1).rotate_deg(90)
        for x in line.images + line.lines + line.collections:
            trans = x.get_transform()
            x.set_transform(r+trans)
            if isinstance(x, PathCollection):
                transoff = x.get_offset_transform()
                x._transOffset = r+transoff

        old = line.axis()
        line.axis(old[2:4] + old[0:2])
        xlabel = line.get_xlabel()
        line.set_xlabel(line.get_ylabel())
        line.set_ylabel(xlabel)
        
        return line

    # Take the highest value for each study
    max_values =[max(inner) for inner in log_data]
    reordered_values=[max_values[value] for value in clustergrid.dendrogram_row.reordered_ind]
    
    max_plot=lineplot_plusplus(data=reordered_values,orientation="vertical",ax=ax2)
    max_plot.set_xlabel('Highest enrichment\n-log(q-value)',size=8)
    max_plot.set_ylabel('Index',size=7)
    max_plot.set_yticks([0, round(0.5*len(included_files))])
    max_plot.xaxis.set_label_position('top')
    max_plot.yaxis.set_label_position('right')
    max_plot.yaxis.tick_right()
    max_plot.xaxis.tick_top()
    max_plot.invert_xaxis()
    max_plot.invert_yaxis()
    max_plot.tick_params(axis='x', labelsize=6)
    max_plot.tick_params(axis='y', labelsize=6)
    max_plot.set_position([0.04,0.03,0.02,0.7])
    max_plot.yaxis.set_label_coords(1.1, 1)
    
    

    plt.savefig(savefig)

# Example loop for creating one clustermap per chromatin state
width_list=[16,12,16,18,14,12,12,22,24,18,14,14,14,14,22]
state_list=['TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh','ZNF-Rpts','Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies']
for i in  range(len(state_list)):
    heatmap('erc2-chromatin15state-all.chart',f'figures/{state_list[i]}_g1.5_6split.pdf',state_list[i],width_list[i])