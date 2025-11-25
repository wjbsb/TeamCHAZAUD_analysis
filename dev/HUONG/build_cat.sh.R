data = readr::read_csv("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/sra.csv")

cols = colnames(data)[c(15:17)]


library(tidyverse)
data_f = data[,c("*library name",cols)]  %>% 
  pivot_longer(cols=cols,names_to = "raw_file",values_to = "values") %>%
  mutate(lane = case_when(
    
    str_detect(string = values,pattern = "L001") ~ "L1",
    str_detect(string = values,pattern = "L002") ~ "L2",
    str_detect(string = values,pattern = "L003") ~ "L3",
    str_detect(string = values,pattern = "L004") ~ "L4",
    str_detect(string = values,pattern = "L005") ~ "L5",
    str_detect(string = values,pattern = "L006") ~ "L6",
    str_detect(string = values,pattern = "L007") ~ "L7",
    str_detect(string = values,pattern = "L008") ~ "L8"
  ),R=case_when(
    str_detect(string=values,pattern="R1") ~ "R1",
    str_detect(string=values,pattern="R2") ~ "R2",
    str_detect(string=values,pattern="R3") ~ "R3"
  ),name_file = paste0("raw_file_",R)) %>%
  pivot_wider(id_cols=c("*library name","R"),names_from = "lane",values_from = "values") %>%
  mutate(
    namefile=paste0(`*library name`,"_",R,".fastq.gz"),
    cat=paste("cat",L1,L2,L3,L4,L5,L6,L7,L8,">",namefile,sep=" "))
 

cat = data_f$cat
write.table(cat, file="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/cat.sh",row.names = F)

cols = colnames(data)[c(1:14)]
data_ff =  data[,cols] %>% unique() %>%
  inner_join(dplyr::select(data_f,`*library name`,namefile,R) %>%
               pivot_wider(id_cols="*library name",names_from="R",values_from="namefile"))

write.table(data_ff, file="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/metadata.txt",row.names = F)
``