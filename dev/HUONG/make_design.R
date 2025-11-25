csv = readxl::read_xlsx("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/seq_template CHAZAUD.xlsx",skip = 32)[2:890,]

library(tidyverse)

csv_like = data.frame(library.name=csv[1:856,1]) %>%
  rename("library.name"="X.library.name")
#   separate(col="library.name",into=c("first","time","age","cells","second","lane"),sep = "_")

design = read_csv("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/design.csv")[,2:5]
l = length(list.dirs(path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq"))
samples=setNames(list.dirs(path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq")[2:l],
                 paste0("X",basename(list.dirs(path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/fastq")[2:l])) %>% str_replace_all("-",".")
)
files=list()
for( s in names(samples)) {
  files[[s]]=data.frame(
    fullpath=list.files(path=samples[[s]],full.names=T)) %>%
    mutate(filename=basename(fullpath) %>% str_replace_all("\\.fastq\\.gz",""),
           sample=s) %>%
    separate(col=filename,into=c("first","time","age","type","second","lane","R","track"),sep="_") %>%
    mutate(library.name=paste(first,time,age,type,second,lane,sep="_")) %>%
    select(sample,library.name) %>%
    unique()
}
final = do.call("rbind",files)
final_f = do.call("rbind",files) %>% filter(library.name %in% csv_like$library.name) %>%
  mutate(library.name=factor(library.name,levels=csv_like$library.name)) %>%
  arrange(library.name)

write_csv(final_f,file="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/HUONG/design_simplified.csv")
