# <span style="color:red"> CORRECTIONS
## LICA-FR sequencing-based expression data
### November 26, 2019

The current sequencing expression data (exp_seq.LICA-FR_corrected.tsv.gz) at [https://dcc.icgc.org/releases/current/Projects/LICA-FR](https://dcc.icgc.org/releases/current/Projects/LICA-FR) contains **incorrect** raw read count values. Please download the corrected "exp_seq.LICA-FR.tsv.gz" file at [https://dcc.icgc.org/releases/Supplementary/LICA-FR/corrected_data](https://dcc.icgc.org/releases/Supplementary/LICA-FR/corrected_data).

# ICGC Data Portal Release 28
This is the Data Portal data release 28 of the International Cancer Genome Consortium (ICGC). Release 28 also contains PCAWG mutation data.

## Projects
This directory contains clinical and analyzed data in the following subdirectories, representing all cancer type datasets available from the ICGC DCC Data Portal for Release 28, processed as of March 27,2019:

 - ALL-US Acute Lymphoblastic Leukemia - TARGET, US
 - AML-US Acute Myeloid Leukemia - TARGET, US
 - BLCA-CN Bladder Cancer - CN
 - BLCA-US Bladder Urothelial Cancer - TCGA, US
 - BOCA-FR Soft Tissue cancer - Ewing sarcoma - FR
 - BOCA-UK Bone Cancer - UK
 - BPLL-FR B-Cell Prolymphocytic Leukemia - FR
 - BRCA-EU Breast ER+ and HER2- Cancer - EU/UK
 - BRCA-FR Breast Cancer - FR
 - BRCA-KR Breast Cancer - Very young women - KR
 - BRCA-UK Breast Triple Negative/Lobular Cancer - UK
 - BRCA-US Breast Cancer - TCGA, US
 - BTCA-JP Biliary Tract Cancer - JP
 - BTCA-SG Biliary Tract Cancer - SG
 - CCSK-US Clear Cell Sarcomas of the Kidney - TARGET, US
 - CESC-US Cervical Squamous Cell Carcinoma - TCGA, US
 - CLLE-ES Chronic Lymphocytic Leukemia - ES
 - CMDI-UK Chronic Myeloid Disorders - UK
 - COAD-US Colon Adenocarcinoma - TCGA, US
 - COCA-CN Colorectal Cancer - CN
 - DLBC-US Lymphoid Neoplasm Diffuse Large B-cell Lymphoma - TCGA, US
 - EOPC-DE Early Onset Prostate Cancer - DE
 - ESAD-UK Esophageal Adenocarcinoma - UK
 - ESCA-CN Esophageal Cancer - CN
 - GACA-CN Gastric Cancer - CN
 - GACA-JP Gastirc Cancer - JP
 - GBM-CN Brain Cancer - Glioblastoma Multiforme - CN
 - GBM-US Brain Glioblastoma Multiforme - TCGA, US
 - HNSC-US Head and Neck Squamous Cell Carcinoma - TCGA, US
 - KICH-US Kidney Chromophobe - TCGA, US
 - KIRC-US Kidney Renal Clear Cell Carcinoma - TCGA, US
 - KIRP-US Kidney Renal Papillary Cell Carcinoma - TCGA, US
 - LAML-CN Leukemia - CN
 - LAML-KR Acute Myeloid Leukemia - KR
 - LAML-US Acute Myeloid Leukemia - TCGA, US
 - LGG-US Brain Lower Grade Glioma - TCGA, US
 - LIAD-FR Benign Liver Tumour - FR
 - LICA-CN Liver Cancer - CN
 - LICA-FR Liver Cancer - FR
 - LIHC-US Liver Hepatocellular carcinoma - TCGA, US
 - LIHM-FR Liver Cancer - Hepatocellular macronodules
 - LINC-JP Liver Cancer - NCC, JP
 - LIRI-JP Liver Cancer - RIKEN, JP
 - LMS-FR Soft tissue cancer - Leiomyosarcoma - FR
 - LUAD-US Lung Adenocarcinoma - TCGA, US
 - LUSC-CN Lung Cancer - CN
 - LUSC-KR Lung Cancer - KR
 - LUSC-US Lung Squamous Cell Carcinoma - TCGA, US
 - MALY-DE Malignant Lymphoma - DE
 - MELA-AU Skin Cancer - AU
 - NACA-CN Nasopharyngeal cancer - CN
 - NBL-US Neuroblastoma - TARGET, US
 - NKTL-SG Blood Cancer - T-cell and NK-cell lymphoma - SG
 - ORCA-IN Oral Cancer - IN
 - OS-US Osteosarcoma - TARGET, US
 - OV-AU Ovarian Cancer - AU
 - OV-CN Ovarian Cancer - CN
 - OV-US Ovarian Serous Cystadenocarcinoma - TCGA, US
 - PAAD-US Pancreatic Cancer - TCGA, US
 - PACA-AU Pancreatic Cancer - AU
 - PACA-CA Pancreatic Cancer - CA
 - PACA-CN Pancreatic Cancer - CN
 - PAEN-AU Pancreatic Cancer Endocrine neoplasms - AU
 - PAEN-IT Pancreatic Endocrine Neoplasms - IT
 - PBCA-DE Pediatric Brain Cancer - DE
 - PBCA-US Pediatric Brain Tumor - Multiple subtypes - CHOP, US
 - PEME-CA Pediatric Medulloblastoma - CA
 - PRAD-CA Prostate Adenocarcinoma - CA
 - PRAD-CN Prostate Cancer - CN
 - PRAD-FR Prostate Cancer - Adenocarcinoma - FR
 - PRAD-UK Prostate Adenocarcinoma - UK
 - PRAD-US Prostate Adenocarcinoma - TCGA, US
 - READ-US Rectum Adenocarcinoma - TCGA, US
 - RECA-CN Renal Cancer - CN
 - RECA-EU Renal Cell Cancer - EU/FR
 - RT-US Rhabdoid Tumors - TARGET, US
 - SARC-US Sarcoma - TCGA, US
 - SKCA-BR Skin Adenocarcinoma - BR
 - SKCM-US Skin Cutaneous melanoma - TCGA, US
 - STAD-US Gastric Adenocarcinoma - TCGA, US
 - THCA-CN Thyroid Cancer - CN
 - THCA-SA Thyroid Cancer - SA
 - THCA-US Head and Neck Thyroid Carcinoma - TCGA, US
 - UCEC-US Uterine Corpus Endometrial Carcinoma- TCGA, US
 - UTCA-FR Uterine Cancer - Carcinosarcoma - FR
 - WT-US Wilms Tumor - TARGET, US

## File Formats

Within each of these directories, you will find data files in tab-delimited (TSV) text format. These files have been compressed using gzip.

The first line of each TSV file contains a header line with the names of each column in the file. 

Release 28 submissions used version 0.19a of the data dictionary:

   http://docs.icgc.org/dictionary/viewer/#?vFrom=0.19a

## File Descriptions

Open-access analyzed data:

* clinical.[ICGC project code].tsv.gz: contains aggregated clinical donor, specimen and sample information
* exp_array.[ICGC project code].tsv.gz: gene expression measured at the transcriptional level (mRNA) using array-based platforms
* exp_seq.[ICGC project code].tsv.gz: gene expression measured at the transcriptional level (mRNA) using sequencing-based platforms
* meth_array.[ICGC project code].tsv.gz: array-based methylation data
* mirna_array.[ICGC project code].tsv.gz: array-based microRNA data
* mirna_seq.[ICGC project code].tsv.gz: sequencing-based microRNA data
* copy_number_somatic_mutation.[ICGC project code].tsv.gz: DNA copy number alterations (ie. gain, losses, LOH) of genes and other loci in tumour tissues relative to normal control samples.
* simple_somatic_mutation.open.[ICGC project code].tsv.gz: open-access simple somatic mutations calls. These include single and multiple base substitutions, and small (<=200bp) insertions and deletions that appear in the tumour tissue, but not in the normal control tissues.
* simple_germline_variations.controlled.[ICGC project code].tsv.gz: contains controlled-access simple germline variations (only available to DACO approved users)
* protein_expression.[ICGC project code].tsv.gz:  translational level expression data
* splice_variant.[ICGC project code].tsv.gz: genomic events that affect the splicing of genes.


## Controlled-Access Analyzed Data

Analyzed germline data is hosted at the ICGC DCC and accessible through the ICGC Data Portal. Analyzed germline data includes variants in cancer donors' normal genomes, such as simple germline variations (SGV), and a very small portion of the simple somatic mutations (SSM) that could leak germline information, ie. alleles that are not in concordance with human reference genome. These SSMs are censored through the germline data masking process (see http://docs.icgc.org/portal/methods/#germline-data-masking for more details):

* Censored SSMs: - made available publicly at DCC Data Portal's Data Repository.

* Original unmasked SSMs: - only available to DACO approved users, who can log in to the Data Portal (top right-hand corner) and download controlled access analyzed data via DCC's Data Repository. These type of files will be marked with "controlled" in the file name. Examples: "simple_somatic_mutation.controlled.[ICGC project code].tsv.gz", "simple_germline_variation.controlled.[ICGC project code].tsv.gz"

## ICGC Publication and Embargo Policy
If you plan to publish using data obtained from this portal please read the [ICGC Publication Policy](http://www.icgc.org/icgc/goals-structure-policies-guidelines/e3-publication-policy).
ICGC Publication guidelines and the current embargo status of each ICGC member project is available at http://docs.icgc.org/portal/publication/#current-moratorium-status-for-icgc-projects . Please also review the publication policy of the individual projects concerned. The publication policy of each project can be found on its project page on the [ICGC website](http://www.icgc.org/icgc/cgp/). For an example, see the [US TCGA Invasive Urothelial Bladder Cancer](http://www.icgc.org/icgc/cgp/88/509/1172) project.

## PCAWG Publication and Embargo Policy
The PanCancer Analysis of Whole Genomes (PCAWG) project is committed to the principle of rapid distribution of primary data, alignments files, and secondary analysis data (such as mutation calls) to all members of the research community that have applied for (and have become approved for) data access. The primary data contributed to the project will fall under any publication/embargo policy of the contributing institute or consortia. The community resource for the PCAWG project comprises newly generated files produced by this collaboration (unified alignments, mutation calls and any secondary analysis data). As is the case for similar "community resource projects", users of data generated and distributed by the PCAWG project are asked to respect the desire of the PCAWG Consortium to publish reports on the generation and analysis of their data. Hence, scientific works with a primary focus on pan-cancer (across entity) analysis that made use of the WGS pan cancer community resource will be under publication embargo until the WGS pan-cancer consortium publishes its marker paper or until July 25, 2019, whichever is earlier. Methodology papers may be published prior to this embargo, with agreement from the full scientific working group. Please write to the ICGC Secretariat (Jennifer Jennings) at [jennifer.jennings@oicr.on.ca](jennifer.jennings@oicr.on.ca) for further information. 


## Access to Raw Data (ie. BAM, FASTQ files)

Raw read files for ICGC non-US projects are accessible at the European Bioinformatics Institute's (EBI) European Genome-phenome Archive (EGA) at https://ega-archive.org/

Note: You must be a DACO approved member to access raw data at EGA.

For ICGC U.S. TCGA projects, raw data will be hosted at the NCI's Genomic Data Commons (GDC) site at https://portal.gdc.cancer.gov/. Application to GDC's DAC is required to gain access to TCGA's controlled data. Please see https://gdc.cancer.gov/access-data/obtaining-access-controlled-data for more details regarding TCGA data.

Please refer to the ICGC DCC documentation page at: http://docs.icgc.org/ for more information.

If you have any questions or concerns, please contact us at [dcc-support@icgc.org](dcc-support@icgc.org)

