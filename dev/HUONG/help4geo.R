setwd("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq/")
library(tidyverse)

dirs = list.dirs("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq",full.names = T)[2:108]
number_samples = lapply(dirs, basename) %>% unlist()

fastqs = lapply(dirs, function(x) list.files(x,full.names=T))
names(fastqs) = number_samples

gz_final=list()
for (s in names(fastqs)) {
  
  sn=str_replace(s,"-",".")
  
  tmp = data.frame(library_name=paste0("X",sn),
                  raw = fastqs[[s]])
  r1=tmp$raw[grepl(pattern = "R1",x=tmp$raw)]
  r2=tmp$raw[grepl(pattern = "R2",x=tmp$raw)]
  r3=tmp$raw[grepl(pattern = "R3",x=tmp$raw)]
  
  tmp = tmp %>% 
    mutate(
      r = case_when(
        str_detect(pattern = "R1",string = raw) ~ "R1",
        str_detect(pattern = "R2",string = raw) ~ "R2",
        str_detect(pattern = "R3",string = raw) ~ "R3")
    )
  R1 = filter(tmp, r == "R1") %>% pull(raw)
  R2 = filter(tmp, r == "R2") %>% pull(raw)
  R3 = filter(tmp, r == "R3") %>% pull(raw)
  
  gz_final[[s]] = data.frame(library_name=paste0("X",sn),
                             r1=R1,r2=R2,r3=R3)
}
gz_final2 = do.call("rbind",gz_final)
samples_order = readr::read_csv("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/samples.txt",col_names = F) %>% unique() %>% unlist() %>% as.character()

gz_final2 = gz_final2 %>% mutate(library_name=factor(library_name,levels=samples_order)) %>% arrange(library_name) 

readr::write_csv(gz_final2,file = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/samples_paths.csv")
