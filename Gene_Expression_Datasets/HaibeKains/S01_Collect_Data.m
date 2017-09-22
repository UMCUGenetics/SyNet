clc;
clear;

%% Download data
% ! wget http://compbio.dfci.harvard.edu/pubs/sbtpaper/data.zip
% ! whet http://compbio.dfci.harvard.edu/pubs/sbtpaper/sbtrobustnessres.zip
% ! unzip data.zip

%% Convert to CSV (run this in R)
%{ 
rm(list=ls())
dir.create("csv_data")
file_list <- list.files("./data/")
for (file_item in file_list) {
    print(sprintf("Reading [%s]", file_item))
    study_name <- strsplit(file_item, "[.]")
    load(sprintf("data/%s", file_item), envir = study <- new.env())
    for (item in ls(study)) {
        n_item <- dim(study[[item]])
        print(sprintf("    Writing variable [%s] with dim [%d, %d, %d]", item, n_item[1], n_item[2], n_item[3]))
        Sys.sleep(0.3)
        write.table(study[[item]], sprintf("./csv_data/Data_%s_%s.csv", study_name[[1]][1], item), sep = "\t", row.names=TRUE, col.names =TRUE, quote = FALSE)
        Sys.sleep(0.3)
    }
    rm(study)
}
%}

