# Pre-process Step

## Step1 Cellranger

```bash
#!/bin/bash
samples_and_batches=(
    "HD_01 HD_01_BCR"
    "HD_02 HD_02_BCR"
    "HD_03 HD_03_BCR"
    "HD_04 HD_04_BCR"
    "HD_05 HD_05_BCR"
    "ITP_01 ITP_01_BCR"
    "ITP_02 ITP_02_BCR"
    "ITP_03 ITP_03_BCR"
    "ITP_04 ITP_04_BCR"
    "ITP_05 ITP_05_BCR"
)

scRNA_transcriptome="/home/data/tmp_data/zyh-BCR/cellrange8.0.1/refdata-gex-GRCh38-2024-A"
vdj_reference="/home/data/tmp_data/zyh-BCR/cellrange8.0.1/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/"
output_dir="/home/data/tmp_data/zyh-BCR/01_cellranger_multi_outputs"
localcores=156
localmem=512
cellranger_path="/home/data/tmp_data/zyh-BCR/cellrange8.0.1/cellranger-8.0.1/cellranger"

trap 'echo "Terminating..."; kill 0' SIGINT

run_cellranger() {
    IFS=' ' read -r scRNA_sample vdj_sample <<< "$1"
    
    echo "Processing sample: $scRNA_sample"
    echo "Batch: $batch"
    
    scRNA_fastqs="/home/data/tmp_data/zyh-BCR/rawdata/RNA-seq/${scRNA_sample}"
    vdj_fastqs="/home/data/tmp_data/zyh-BCR/rawdata/BCR/${vdj_sample}"
    
    if [ ! -d "$scRNA_fastqs" ] || [ ! -d "$vdj_fastqs" ]; then
        echo "Directory $scRNA_fastqs or $vdj_fastqs does not exist. Skipping sample $scRNA_sample."
        return
    fi
    
    output_dir_sample=${output_dir}/${scRNA_sample}_run
    mkdir -p ${output_dir}/csv
    
    scRNA_samples=$(ls ${scRNA_fastqs} | grep -oE "${scRNA_sample}-[0-9]+" | sort -u)
    
    echo "scRNA Samples: $scRNA_samples"
    echo "vdj Sample: $vdj_sample"
    
    csv_sample="${output_dir}/csv/${scRNA_sample}.csv"
    echo "[gene-expression]" > ${csv_sample}
    echo "reference,${scRNA_transcriptome}" >> ${csv_sample}
    echo "create-bam,true" >> ${csv_sample}
    echo "" >> ${csv_sample}
    echo "[vdj]" >> ${csv_sample}
    echo "reference,${vdj_reference}" >> ${csv_sample}
    echo "" >> ${csv_sample}
    echo "[libraries]" >> ${csv_sample}
    echo "fastq_id,fastqs,feature_types" >> ${csv_sample}
    
    for sample in ${scRNA_samples}; do
        echo "${sample},${scRNA_fastqs},gene expression" >> ${csv_sample}
    done
    
    echo "${vdj_sample},${vdj_fastqs},vdj" >> ${csv_sample}
    
    cat ${csv_sample}

    if [ -e ${output_dir_sample}/_lock ]; then
        echo "Removing existing lock file for ${output_dir_sample}"
        rm -f ${output_dir_sample}/_lock
    fi

    ${cellranger_path} multi \
        --id=${scRNA_sample}_run \
        --csv=${csv_sample} \
        --localcores=${localcores} \
        --localmem=${localmem}
}

export -f run_cellranger
export scRNA_transcriptome
export vdj_reference
export output_dir
export localcores
export localmem
export cellranger_path

parallel -j 4 run_cellranger ::: "${samples_and_batches[@]}"
```

## Step2 Seurat

```R
setwd("/home/data/tmp_data/zyh-BCR/")

rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library (harmony)
library(tidyverse)
library(patchwork)

assays <- dir("./02_10X_data/")
dir <- paste0("./02_10X_data/", assays)
samples_name = c("H01", "H02", "H03", "H04", "H05", "P01", "P02", "P03", "P04", "P05")

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #Insufficient data values to produce 24 bins.  
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
}

setwd("/home/data/tmp_data/zyh-BCR/CCA/")
dir.create("Integrate")
setwd("./Integrate")

names(scRNAlist) <- samples_name
#system.time(save(scRNAlist, file = "Integrate/scRNAlist0.Rdata")) 
system.time(saveRDS(scRNAlist, file = "scRNAlist0.rds"))
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA

table(scRNA$orig.ident)
#save(scRNA,file = 'scRNA_orig.Rdata')
saveRDS(scRNA,file = 'scRNA_orig.rds')

theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "orig.ident"
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

minGene=200
minUMI=500
pctMT=5
pctHB=1

scRNA <- subset(scRNA, subset = nCount_RNA > minUMI & nFeature_RNA > minGene & 
                percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8) 

#scRNA_initial <- scRNA

all.genes <- rownames(scRNA)
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(features = all.genes)
scRNA <- RunPCA(scRNA, verbose = F)

#JackStrawPlot(scRNA, dims = 1:50)

ElbowPlot(scRNA, ndims = 50)

scRNA <- FindNeighbors(scRNA, dims = 1:30, reduction = "pca")
scRNA <- FindClusters(scRNA, resolution = 1, cluster.name = "unintegrated_clusters")

scRNA <- RunUMAP(scRNA, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scRNA, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
p <- DimPlot(scRNA, group.by = "orig.ident")
ggsave("UMAP_Samples.pdf", p, width = 8, height = 6)

p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("UMAP_Samples_Split.pdf", p, width = 18, height = 12)
saveRDS(scRNA, "scRNA.rds")


scRNA <- IntegrateLayers(
  object = scRNA, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  verbose = FALSE)

scRNA <- FindNeighbors(scRNA, reduction = "integrated.cca", dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.9, cluster.name = "cca_clusters")
scRNA <- RunUMAP(scRNA, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

## DoubletFinder
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#ls("package:DoubletFinder")

scRNA <- JoinLayers(scRNA)
sweep.res.list <- paramSweep(scRNA, PCs = 1:50, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()
 
p <- ggplot(bcmvn, aes(x=pK, y=BCmetric, group=1)) + 
  geom_point() +
  geom_line()
 
ggsave("pk_plot.pdf", p, width = 8, height = 6)
annotations <- scRNA$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
print(homotypic.prop)
nExp_poi <- round(0.07*nrow(scRNA@meta.data))    
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scRNA <- doubletFinder(scRNA, PCs = 1:50, 
                       pN = 0.25, pK = pk_best, nExp = nExp_poi, # nExp = nExp_poi.adj, 
                       reuse.pANN = FALSE, sct = FALSE)
 
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == "pANN_0.25_0.03_6637"] <- "Double_score"
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == "DF.classifications_0.25_0.03_6637"] <- "Is_Double"
head(scRNA@meta.data[, c("Double_score", "Is_Double")])

pdf("doublefinder_result.pdf", width =12, height = 8)
DimPlot(scRNA, reduction = "umap", group.by = "Is_Double")
VlnPlot(scRNA, group.by = "Is_Double", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0, ncol = 2)
dev.off()

scRNA <- subset(scRNA, subset = Is_Double == "Singlet")

#NOTICE: Run the former steps again after processed withDoubletFinder

pdf("marker_t.pdf", width = 8, height = 8)
DotPlot(scRNA, features =  c("CD3G", "CD3E"), group.by ="cca_clusters") + RotatedAxis()
dev.off()

scRNA_subset <- subset(scRNA, cells = WhichCells(scRNA, idents = setdiff(unique(Idents(scRNA)), "21")))
subset_cells <- colnames(scRNA_subset)
scRNA_final <- subset(scRNA_initial, cells = subset_cells)
print(scRNA_final)
scRNA <- scRNA_final
```

## Step3 Dandelion

### Step3.1 Dandelion Pre-process

```bash
singularity pull library://kt16/default/sc-dandelion:latest

sample_list=("H01" "H02" "H03" "H04" "H05" "P01" "P02" "P03" "P04" "P05")

for sample in "${sample_list[@]}"
do
    vdj_dir="/home/data/tmp_data/zyh-BCR/01_cellranger_multi_outputs/${sample}_run/outs/per_sample_outs/${sample}_run/vdj_b/"
    destination_dir="/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/10x_vdj/${sample}"

    mkdir -p "$destination_dir"

    cp "${vdj_dir}filtered_contig_annotations.csv" "$destination_dir"
    cp "${vdj_dir}filtered_contig.fasta" "$destination_dir"
done

cd /home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/10x_vdj/
singularity run -B $PWD /home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/sc-dandelion_latest.sif dandelion-preprocess --meta meta.csv --file_prefix filtered 
```

### Step3.2 Dandelion Process

```python
import os
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
import functools
import seaborn as sns
import scipy.stats
import anndata
import dandelion as ddl
ddl.logging.print_header()

os.chdir(os.path.expanduser('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/process/'))
# I'm importing scanpy here to make use of its logging module.
sc.settings.verbosity = 3
warnings.filterwarnings('ignore')
sc.logging.print_header()

samples = ["H01", "H02", "H03", "H04", "H05", "P01", "P02", "P03", "P04", "P05"]
vdj_list = []
for sample in samples:
    vdj = ddl.read_10x_airr('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/10x_vdj/'+sample+'/dandelion/filtered_contig_dandelion.tsv')
    #the dandelion output already has the sample ID prepended at the start of each contig
    vdj_list.append(vdj)
vdj = ddl.concat(vdj_list)

#vdj.write('vdj_concat.h5ddl')
#vdj = ddl.read_h5ddl("/home/data/tmp_data/zyh-BCR/CCA/01_dandelion/process/vdj-processed.h5ddl")
adata = sc.read("/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/object.h5ad", gex_only=True)
adata = sc.read("/home/data/tmp_data/zyh-BCR/CCA/01_dandelion/process/gex-processed.h5ad", gex_only=True)
old_rownames = adata.obs.index
new_rownames = [name.split('-', 1)[0] for name in old_rownames]
adata.obs.index = new_rownames
print(adata.obs.head())
adata.obs.to_csv('adata_obs.csv')
#vdj.data['cell_id'] = vdj.data['cell_id'].apply(lambda x: x.rsplit('_', 1)[0])
vdj, adata = ddl.pp.filter_contigs(vdj, adata)
vdj, adata = ddl.pp.check_contigs(vdj, adata)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.1)
#sc.settings.set_figure_params(dpi=80)

ddl.tl.find_clones(vdj)
#ddl.tl.generate_network(vdj, key = "junction", num_cores = 64)
ddl.tl.generate_network(vdj, num_cores = 64)

ddl.tl.clone_size(vdj)
ddl.tl.clone_size(vdj, max_size=3)

ddl.tl.transfer(adata, vdj)

adata.write("gex-processed.h5ad")
vdj.write("vdj-processed.h5ddl")
```
