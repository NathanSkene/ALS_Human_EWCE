#############################################################
## GSE18920 --- Splicing array data from human spinal cord ##
#############################################################
library(limma)
library(One2One)
library(tidyr)
library(EWCE)
library(R.utils)
library(drake)
library(ggplot2)
library(reshape)
sourceDirectory("R")
#source("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/generate.bootstrap.plots.for.transcriptome_SOD1.r")

mkdir("Results")
mkdir("Results/Tables")
mkdir("Results/Figures")

load(file="Data/Tidy/ctd_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")

plan <- drake_plan(
  path = "/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/",
  data_spinal = load_als_data(path),
  annot_spinal = load_als_annot(path,data_spinal),
  tt_spinal = run_als_diffExp_analysis(annot_spinal,data_spinal),
  ewce_spinal = run_ewce(ctd,tt_spinal),
  catchOut = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=tt_spinal,thresh=250,annotLevel=1,reps=1000,full_results=ewce_spinal,listFileName="Human ALS Spinal")
)

config <- drake_config(plan)
vis_drake_graph(config)

make(plan)

# loadd(ctd,ewce_spinal,tt_spinal)