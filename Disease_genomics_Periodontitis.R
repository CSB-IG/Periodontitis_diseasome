
# Installation instructions
#library(devtools)
#install_bitbucket("ibi_group/disgenet2r", force=TRUE)

library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "ehernandez@inmegen.gob.mx", 
  password = "********" )

Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

setwd("/Users/enrique/Dropbox/Periodontitis_Genomics")

POD <- disease2gene(disease = "C0031099", 
                     database = "ALL",
                     score = c(0.3,1) 
)

POD

plot(POD,
     class = "Network",
     prop = 10)

plot(POD,
     class="ProteinClass")

#??"DataGeNET.DGN"

# PDA <-extract(POD)
# 
# PDA

PDD <-POD@qresult

## PDD <- POD@qresult  ## alternative to extract for a DataGeNET.DGN object


PDD

## Diseases_sharing_genes_with Periodontitis

PGSD <- c("C0031099","C0948089","C0023465","C0023479",
          "C0001418","C0001973","C0002395","C0002726",
          "C0002871","C0038013","C0003467","C0004096",
          "C0004153","C0004238","C0006840","C4048328",
          "C0024117","C1561643","C1956346","C0010068",
          "C0010346","C0011581","C0242339","C0376618",
          "C0014859","C0031069","C0017661","C0018213",
          "C0018520","C0278996","C0018965","C0019080",
          "C0019158","C1857728","C0019348","C0019829",
          "C0020429","C0020456","C0085580","C0021390",
          "C0021655","C0684249","C0024796","C0023418",
          "C0026640","C0026764","C0027051","C0400966",
          "C0028754","C0585362","C0029456","C0235974",
          "C0030297","C0030567","C1704436","C0001623",
          "C0032285","C0032460","C0033687","C0033860",
          "C0003873","C0036690","C1527336","C0520679",
          "C0038454","C0024141","C0040034","C0011854",
          "C0009324","C0346153","C0006142","C0678222",
          "C0011860","C1852092","C0011849","C0019693",
          "C0376358","C2931456","C4722327","C0015499")

PODNet <- disease2gene(
  disease = PGSD,
  database = "CURATED",
  score =c(0.4,1),
  verbose  = TRUE )


plot(PODNet,
     class = "Network",
     prop = 2)

plot(PODNet,
     class="Heatmap",
     limit =30,
     cutoff=0.2)

plot(PODNet,
     class="ProteinClass")

# PODNet

POT <-extract(PODNet)

#POT

write.csv(POT, file="Diseases_Genes_Periodontitis_extract.csv")

PDRS<-PODNet@qresult

#PDRS

write.csv(PDRS, file="Diseases_Genes_Periodontitis.csv")

#### Diseasome networks

DisPer <- disease2disease_by_gene(
  disease  = "C0031099", 
  database = "CURATED", ndiseases = 100 )

DPeri <- extract(DisPer)

plot(DisPer, 
     class = "Network",
     layout="layout.lgl" ,
     prop = 0.2 )

plot( DisPer, 
      class = "Diseasome",
      layout="layout.lgl" ,
      prop = 0.1 )

plot( DisPer,
      class="Barplot")

write.csv(DPeri, file="Top_Diseases_shared_genes_Periodontitis.csv")

PerDis <-DPeri[c("disease1_name", "disease2_name","jaccard_genes","ngenes", "pvalue")]

write.csv(PerDis, file="Top_Diseases_shared_genes_Periodontitis_abridged.csv")




