# lung-cancer-CB
In this repository the necessary code and data to carried out the experiments proposed in Heterogeneous Gene Expression Cross-Evaluation of Robust Biomarkers
Using Machine Learning Techniques Applied to Lung Cancer are available. Next a description of available files is presented.

LUAD_pacientes.csv:file to load counts files to R. TCGA-LUAD.

LUSC_files.csv:file to load counts files to R. TCGA-LUSC.

counts: Counts files.

LUAD.cvs: preprocessed RNA-Seq expression matrix (without batch effect treat) for ACC and control samples from TCGA-LUAD project

LUSC.csv: preprocessed RNA-Seq expression matrix (without batch effect treat) for SCC and control samples from TCGA-LUSC project

tri.csv: preprocessed RNA-Seq expression matrix (without batch effect treat) for SCC and control samples from TCGA-LUAD and TCG-LUAD project

mic_exp1.csv: preprocessed microarray expression matrix (All series integrated)

mic_samples.csv: labels for mic_exp samples.

code.R: code to develop the experiments

def.RData: Workspace with the results obtained. 

Data files can be found in drive repository https://drive.google.com/drive/folders/1bJ98M8QWRb1AZpJWMBo9cplxESYkf0ES?usp=sharing because of their size.


# SHINY CODE #

app.R: code to contruct your local intuitive interface to carry out the experiments. To run the application download data.RData, www folder and app.R from https://drive.google.com/drive/folders/1mmOIbRgZWTgdHOywniGx8MDwuY-_CID1?usp=sharing  and run it on R. (All files should be in the same folder)


