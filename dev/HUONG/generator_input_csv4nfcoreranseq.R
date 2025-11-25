
path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq"
samples = "2-428077|2-428080|2-428081|2-428082|2-428083|2-428084|2-428085|2-428086|2-428087|2-428088|2-428089|2-428090|2-428091|2-428093|2-428094|2-428095|2-428098|2-428099|2-428152|2-428153|2-428154|2-428155|2-428156|2-428158"
dirs = lapply(list.dirs(path), function(x) list.files(x,full.names = T))[2:129]
names(dirs) = str_replace_all(list.dirs(path)[2:129],paste0(path,"/"),"")

dirs_f = dirs[str_split(samples,"\\|")[[1]]]

dirs_ff = data.frame(do.call("rbind",dirs_f)) %>% 
  rownames_to_column(var="dataset") %>%
  pivot_longer(cols=paste0("X",seq(1,24)),names_to = "files",values_to="path_files") %>%
  mutate(
    files = case_when(
      str_detect(path_files,"_R1_001") ~ "fastq_1",
      str_detect(path_files,"_R2_001") ~ "index",
      str_detect(path_files,"_R3_001") ~ "fastq_2"),
    lanes = case_when(
      str_detect(path_files,"_L001_") ~ "L1",
      str_detect(path_files,"_L002_") ~ "L2",
      str_detect(path_files,"_L003_") ~ "L3",
      str_detect(path_files,"_L004_") ~ "L4",
      str_detect(path_files,"_L005_") ~ "L5",
      str_detect(path_files,"_L006_") ~ "L6",
      str_detect(path_files,"_L007_") ~ "L7",
      str_detect(path_files,"_L008_") ~ "L8")
  ) %>% filter(files != "index") %>%
  pivot_wider(id_cols=c("dataset","lanes"),names_from="files",values_from = "path_files") %>%
  mutate(strandedness="auto") %>%
  dplyr::select(-lanes) %>%
  rename("sample"="dataset")

write.csv(dirs_ff,file = "/media/mna_bioinfo/DD_linux/SCRIPTS/16042024/input_nfcorernaseq_3.14.0.16042024.csv",row.names = F)
