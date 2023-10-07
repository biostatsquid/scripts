########################################
#### PART 1: CHEMBRIDGE DATASET ####
# Objective: preprocessing compound library datasets with chemical features for each compound
# 
# 


########################################


##### 1. IMPORTING RELEVANT LIBRARIES ###################
require(dplyr)
require(tidyverse)
require(readxl) #readxl allows to work with excel workbook with several worksheets
require(tools)
require(ggplot2)
#require(ggVennDiagram)

#path
path <- "C:/Users/laura/Documents/Biostatsquid/Scripts/Compound_dataset"
setwd(path)

######2. Allocating more memory for next steps ###################

memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size = 56000)


######3.DATA IMPORT AND INITIAL CLEAN-UP 
############3.1. Importing datasets #################
#'Description of match_names
#'function to import the features generated in OCHEM and match the names (compound ids)
#'file == name of OCHEM generated .xlsx file (from original .sdf). Note. before saving
#'file as.xlsx, change col name with compound IDs to "id"
#'names == name of .csv exported from MarvinView (from original .sdf)
#'library == which library the files are coming from, e.g. "ChemBridge"

match_names <- function(file, names, library){
  # 1. Extract names from sdf file of the given library (i.e. FDA)
  names_file <- read.csv(names)
  cID <- names_file$id #compound ID, dataframe with 1 column with the compound names
  
  # 2. Load descriptors (.xlsx) 
  sheetnames <- excel_sheets(file) #create list of all of the sheets
  all_features <- sheetnames %>% map(~read_xlsx(file, .)) #join all tables by row (with map_dfc original column headers were not kept)
  all_features_df <- bind_cols(lapply(all_features, as.data.frame)) #convert to dataframe
  
  colnames(all_features_df) <- tolower(colnames(all_features_df)) #all column names in lowercase
  
  #3. Bind names to features
  final <- cbind(cID, all_features_df) #add names dataframe to first column
  final$cID = as.factor(as.character(final$cID)) # cID as factors
  final <- add_column(final, library = as.factor(library), .after = 2) #add library columnn as factors in 3rd column
  
  return(final)
}


# fda <- match_names("descriptors_FDA_xlsx.xlsx", "FDA_sdf.csv", "FDA") # 1280 x 13477
# chembridge <- match_names("descriptors_ChemBridge_xlsx.xlsx", "ChemBridge_sdf.csv", "ChemBridge") # 5000 x 16251
# ppi <- match_names("descriptors_PPI_xlsx.xlsx", "PPI_sdf.csv", "PPI") # 5000 x 14079
# chemdiv0 <- match_names("descriptors_ChemDiv_0_xlsx.xlsx", "ChemDiv_0_sdf.csv", "ChemDiv") # 5001 x 17563
# chemdiv1 <- match_names("descriptors_ChemDiv_1_xlsx.xlsx", "ChemDiv_1_sdf.csv", "ChemDiv") # 4999 x 16278
# enamine0 <- match_names("descriptors_Enamine_0_xlsx.xlsx", "Enamine_0_sdf.csv", "Enamine") # 5001 x 15808
# enamine1 <- match_names("descriptors_Enamine_1_xlsx.xlsx", "Enamine_1_sdf.csv", "Enamine") # 4997 x 15414
chembridge <- match_names("descriptors_ChemBridge_xlsx1.xlsx", "ChemBridge_sdf1.csv", "ChemBridge") # 199 x 16251

# dfs <- list(chembridge, ppi, fda, chemdiv0, chemdiv1, enamine0, enamine1) # list of databases
# alldata <- bind_rows(dfs) # bind all data together 
alldata <- chembridge


# 
# alldata$library[alldata$cID == "Prestw-919"] <- "FDA"
# alldata <- alldata[alldata$library != "FDA_194", ]
# #alldata$`mw:(alvadesc)`[alldata$cID == "Spiramycin"]
# alldata <- alldata[-11280, ] #remove spyramicin wrong entry
# #summary(as.factor(alldata$library)) 
# print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) #31279 x 28803

#Export if needed
#write.csv(alldata, file = "alldata_withspir_withprestw919.csv")
#alldata <- read.csv(file = "alldata_withspir_withprestw919.csv")

############3.2. Data filtering ########################

alerts_changename <- function(df){
  i <- 0
  n <- ncol(df)
  for (j in (1:n)){
    if (grepl("^alert", colnames(df)[j])){
      colnames(df)[j] <-  paste0("alert_", i) #most columns names are "Alert0_... very long name" - reduce to "Alert0", "Alert1"...
      i <- i + 1
    }
  }
  return(df)
} #most columns names are "Alert0_... very long name" - reduce to "Alert0", "Alert1"...
alldata <- alerts_changename(alldata)


#Dealing with missing values: if na for all compounds (rows), drop columns
alldata <- alldata[, !apply(alldata, 2, function(x) all(is.na(x)))] # no all na features 
print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) #31279 obs x 28795 variables

## Remove columns with more than 90% NA
alldata <- alldata[, which(colMeans(!is.na(alldata)) > 0.9)] 
print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) #31279 obs x 10683 variables


#if the same value for all compounds (rows), drop columns
alldata <- alldata[vapply(alldata, function(x) length(unique(x)) > 1, logical(1L))] 
print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) # 31279 obs x 10679 variables

#remove rows with all nas
alldata <- alldata[!!rowSums(!is.na(alldata[,-c(which(colnames(alldata) == "cID"), 
                                                which(colnames(alldata) == "library"),
                                                which(colnames(alldata) == "smiles"))])), ] 
print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) #31279 obs x 10679 variables


#Remove duplicated rows (without taking into account library column)
#alldata <-  distinct(.data = alldata[, -c(which(colnames(alldata) == "library"))], .keep_all = TRUE)
#duprows <- alldata[duplicated(alldata[,-which(colnames(alldata) == "library")]) | duplicated(alldata[,-which(colnames(alldata) == "library")], fromLast = TRUE),]
alldata <- alldata %>% distinct(cID, smiles, .keep_all = TRUE)
print(paste("Nrows = ", nrow(alldata), "Ncols = ", ncol(alldata))) #31052 x 10679
row.names(alldata) <- 1:nrow(alldata) #reset rownames after deleting 227

#MWE dataset: 199 x 5656

############3.3. Checking for duplicates ####
#Checking colnames in common and cID in common
#'NOTE. Chemdiv0 and ppi have 24 identical observations (cID and smiles)
#'Chemdiv1 and ppi have 203 identical observations (cID and smiles) and 2
#'observations with = cID but different smiles

#all_colnames <- list(fda = colnames(fda), chembridge = colnames(chembridge), 
#                     chemdiv0 = colnames(chemdiv0), chemdiv1 = colnames(chemdiv1),
#                     enamine0 = colnames(enamine0), enamine1 = colnames(enamine1),
#                     ppi = colnames(ppi))

#all_cID <- list(fda = as.character(fda$cID), chembridge = as.character(chembridge$cID), 
#                chemdiv0 = as.character(chemdiv0$cID), chemdiv1 = as.character(chemdiv1$cID),
#                enamine0 = as.character(enamine0$cID), enamine1 = as.character(enamine1$cID),
#                ppi = as.character(ppi$cID))

#installation
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")

#library("ggVennDiagram")
#ggVennDiagram(all_colnames, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "light blue", high = "yellow")

#ggVennDiagram(all_cID, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "light blue", high = "yellow")

#cID_0 <- intersect(chemdiv0$cID, ppi$cID) #observations in common (cID)
#cID_1 <- intersect(chemdiv1$cID, ppi$cID)

#smiles_0 <- intersect(chemdiv0$smiles, ppi$smiles) #observations in common (smiles)
#smiles_1 <- intersect(chemdiv1$smiles, ppi$smiles)

#Checking chemdiv0 and ppi
#testing0 <- data.frame(cID_0)
#for (i in (1:length(cID_0))){
#  testing0[i, 2] <- ppi$smiles[ppi$cID == cID_0[i]]
#  testing0[i, 3] <- chemdiv0$smiles[chemdiv0$cID == cID_0[i]]
#}
#testing0
#testing0$V2 == testing0$V3

#Checking chemdiv1 and ppi
#testing1 <- data.frame(cID_1)
#for (i in (1:length(cID_1))){
#  testing1[i, 2] <- ppi$smiles[ppi$cID == cID_1[i]]
#  testing1[i, 3] <- chemdiv1$smiles[chemdiv1$cID == cID_1[i]]
#}


#testing1
#testing1$cID[testing1$V2 != testing1$V3]
#chemdiv1$smiles[chemdiv1$cID == "C466-0123"] == ppi$smiles[ppi$cID == "C466-0123"]
#chemdiv1$smiles[chemdiv1$cID == "C466-0146"] == ppi$smiles[ppi$cID == "C466-0146"]

#intersect(alldata$cID[alldata$library == "PPI"], alldata$cID[alldata$library == "ChemDiv"])
#alldata[alldata$cID == "C466-0123", 1:2]
#alldata[alldata$cID == "C466-0146", 1:2]

# levels(alldata$cID) <- c(levels(alldata$cID), "C466-0123_dup", "C466-0146_dup")
# alldata$cID[alldata$cID == "C466-0123" & alldata$library == "ChemDiv"] <- "C466-0123_dup"
# alldata$cID[alldata$cID == "C466-0146" & alldata$library == "ChemDiv"] <- "C466-0146_dup"

#checking
#intersect(alldata$cID[alldata$library == "PPI"], alldata$cID[alldata$library == "ChemDiv"])

#alldata$smiles[alldata$cID == "C466-0123"]
#alldata$smiles[alldata$cID == "C466-0123_dup"]

#alldata$smiles[alldata$cID == "C466-0146"]
#alldata$smiles[alldata$cID == "C466-0146_dup"]


############3.4. Removing features with Near Zero Variance (not applied) ######################
#see nearZeroVar for details
#removes too many features, won't use it yet

#library(caret)
# get indices of data.frame columns (pixels) with low variance
#badCols <- nearZeroVar(alldata) #takes around 15 min
#print(paste("Fraction of nearZeroVar columns:", round(length(badCols)/length(alldata),4))) # 0.882
# remove those "bad" columns from the dataframe
#alldata_nzv <- alldata[, -badCols]

#alldata2 <- alldata #saving point






############3.5. Duplicated columns ##############################

#if all the values for all compounds are identical between two or more descriptors/columns, 
#drop duplicate columns while keeping representative column with joined column name

#'Description
#'finds duplicated columns (all values in the column), merges the names of duplicated 
#'columns into one (sep = ;), keeps only column with merged names and discards the rest
#'input dataframe
#'output list of duplicated and non-duplicated df with merged col names

dupcols <- function(df = testframe){
  
  require(digest)
  
  #separate into 2 dataframes; with and without duplicated columns
  nondups <- df[!duplicated(lapply(df, digest))]
  dups <- df[duplicated(lapply(df, digest))]
  
  for(j in 1:ncol(dups)){
    if (j %% 100 == 0){
      print(paste(j, "out of ", ncol(dups)))
    }
    for(i in 1:ncol(nondups)){
      test <- nondups[,i] == dups[,j]
      if(FALSE %in% test | all(is.na(test))){
        NULL
      }
      else if(all(test == TRUE, na.rm = TRUE)){
        names(nondups)[i] <- paste(as.character(names(nondups[i])), as.character(names(dups[j])), sep = "; ")
      }
    }
  }
  
  return(list(df1 = nondups, df2 = dups))
}
Sys.time()
dups_nodups <- dupcols(alldata) #takes aprox 1h (1 min for MWE)
Sys.time()
nodu <- dups_nodups$df1
#dups <- dups_nodups$df2
print(paste("Nrows no dups = ", nrow(nodu), "Ncols no dups = ", ncol(nodu))) #31052 x 7574

longnames_nodu <- colnames(nodu)
#write.csv(longnames_nodu, "~/MWE\\longnames_nodu.csv") #export merged column names to csv
colnames(nodu) <-  sub("\\; .*", "_dup", colnames(nodu)) #columns that are duplicated will have the name of the first column + _dup to shorten names

#saveRDS(nodu, file = "nodu.Rds")
#saveRDS(alldata, file = "alldata_afterpreprocessing.Rds")
##### RESULTING DATAFRAME nodu:
#Dataframe with no duplicated, no all na rows, no na cols, no >90% na cols, no cols
#with all = value, no duplicated columns, no duplicated rows
#Colnames in lowercase letters, combination of all duplicated column names if 
#applicable, also shortened names ("Alert0_" instead of "Alert0_kjnjasdnfkasdn")
#columns that were duplicated were renamed to "..._dup" - all names are stored in longnames_nodups


#rm(list = c("longnames_nodu", "chembridge", "alldata", "chemdiv1", "enamine0", "enamine1", "fda", "ppi", "chemdiv0", "dfs", "dups_nodups"))

#save.image("~/MHB/ICB/Data/HTS_transfer/data/Files_preprocessing/features_clean20.RData")

#write.csv(nodu, "~/MHB/ICB/Data/HTS_transfer/data/Files_preprocessing\\nodu.csv") #save nodu as csv


# change feature columns to numeric for PCA
#nodu <- readRDS(file = "nodu.Rds") #31052 x 7574
row.names(nodu) <- 1:nrow(nodu)
#nodu[c(134504, 134597, 134707, 135001, 135151), 1:5] #there are 5 rows with just NAs
#nodu <- nodu[-c(10300, 10393, 10503, 10797, 10947), ] #31047 x 7574

#colnames(nodu)[colSums(is.na(nodu)) > 0]
#2 compounds missing smiles, enter manually
#nodu$smiles[10797] <- "CCC(C)C1C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NCCCCC(C(=O)NC(C(=O)N1)CCCN)NC(=O)C(C(C)CC)NC(=O)C(CCC(=O)O)NC(=O)C(CC(C)C)NC(=O)C2CSC(=N2)C(C(C)CC)N)CC(=O)N)CC(=O)O)CC3=CN=CN3)CC4=CC=CC=C4"
#nodu$smiles[10947] <- "CC(C)(C)CC(C)(C)C1=CC=C(C=C1)O.C=O.C1CO1"

#nodu <- na.omit(nodu) # some rows with NAs #29771 x 7558 #will not do this since it deletes all FDA rows
#nodu <- nodu[ ,colSums(is.na(nodu)) == 0] #31047 x 5509

#Checking na columns
# listnas <- list(nodu = colnames(nodu[ ,colSums(is.na(nodu)) != 0]), 
#      spir = colnames(nodu[which(is.na(nodu[11275,]))]),
#      prestw = colnames(nodu[which(is.na(nodu[11276,]))]))
# ggVennDiagram(listnas, label_alpha = 0) + 
#   ggplot2::scale_fill_gradient(low = "light blue", high = "yellow")


#nodu <- nodu[, !apply(nodu, 2, function(x) all(is.na(x)))]
#nodu2 <- nodu

columns <- colnames(nodu)[-c(1:4, 4473)] #select columns except first 4 and last ("cID", "smiles","library","moleculeid", "inchi")
nodu[columns] <- sapply(nodu[columns],as.numeric)

#Checking for NAs before PCA 
nasperrow <- apply(nodu, MARGIN = 1, function(x) sum(is.na(x)))
summary(nasperrow) #should all be 0
which(is.na(nodu))




############3.6. Checking the LIPINSKI RULE ############################
library(ggVennDiagram)
all_lipinski <- list(donors = as.character(nodu$cID[nodu$`nhdon:(alvadesc)` > 5]), 
                     acceptors = as.character(nodu$cID[nodu$`nhacc:(alvadesc)` > 10]), 
                     mw = as.character(nodu$cID[nodu$`mw:(alvadesc)` > 500]),
                     mlogp = as.character(nodu$cID[nodu$`mlogp:(alvadesc)` > 4.15]),
                     lipinski = as.character(nodu$cID[nodu$`ro5:(alvadesc)` == 1]))
#library("ggVennDiagram")
ggVennDiagram(all_lipinski, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "white", high = "yellow")

# lipinski_df <- nodu[nodu[, "ro5:(alvadesc)"] == 1, c("cID", "ro5:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)", 
#                                       "nhacc:(alvadesc)", "nhdon:(alvadesc)")]
# 
# lipinski_df[11,] <- nodu[nodu[, "cID"] == "C200-5954", c("cID", "ro5:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)", 
#                                             "nhacc:(alvadesc)", "nhdon:(alvadesc)")]
# lipinski_df[12,] <- nodu[nodu[, "cID"] == "C694-0016", c("cID", "ro5:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)", 
#                                                          "nhacc:(alvadesc)", "nhdon:(alvadesc)")]
# lipinski_df[13,] <- nodu[nodu[, "cID"] == "D268-0425", c("cID", "ro5:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)", 
#                                                          "nhacc:(alvadesc)", "nhdon:(alvadesc)")]



# #Finding false positives of lipinski failure 
# nodu$cID[nodu$`ro5:(alvadesc)` == 1 & nodu$`nhdon:(alvadesc)` > 5 & 
#            nodu$`nhacc:(alvadesc)` <= 10 & nodu$`mw:(alvadesc)` <= 500] #Pemetrexed disodium
# nodu$cID[nodu$`ro5:(alvadesc)` == 1 & nodu$`nhdon:(alvadesc)` <= 5 & 
#            nodu$`nhacc:(alvadesc)` <= 10 & nodu$`mw:(alvadesc)` > 500 &
#            nodu$`mlogp:(alvadesc)` <= 4.15] #Reserpine and Verteporfin
# nodu$cID[nodu$`ro5:(alvadesc)` == 1 & nodu$`nhdon:(alvadesc)` <= 5 & 
#            nodu$`nhacc:(alvadesc)` <= 10 & nodu$`mw:(alvadesc)` <= 500 &
#            nodu$`mlogp:(alvadesc)` > 4.15] # Halofantrine hydrochloride
# 
# #Finding false negatives of lipinski failure 
# nodu$cID[nodu$`ro5:(alvadesc)` == 0 & nodu$`nhdon:(alvadesc)` <= 5 & 
#            nodu$`nhacc:(alvadesc)` > 10 & nodu$`mw:(alvadesc)` > 500 &
#            nodu$`mlogp:(alvadesc)` <= 4.15] #C200-5954 C694-0016
# nodu$cID[nodu$`ro5:(alvadesc)` == 0 & nodu$`nhdon:(alvadesc)` <= 5 & 
#            nodu$`nhacc:(alvadesc)` > 10 & nodu$`mw:(alvadesc)` <= 500 &
#            nodu$`mlogp:(alvadesc)` > 4.15] #D268-0425
# 
# 
# library(tableHTML)
# lipinski_conditional <- tableHTML(lipinski_df) %>%
#   add_css_conditional_column(conditional = "==",
#                              value = 1, 
#                              css = list('background-color', "red"),
#                              columns = c("ro5:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 4.15, 
#                              css = list('background-color', "red"),
#                              columns = c("mlogp:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 500, 
#                              css = list('background-color', "red"),
#                              columns = c("mw:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 10, 
#                              css = list('background-color', "red"),
#                              columns = c("nhacc:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 5, 
#                              css = list('background-color', "red"),
#                              columns = c("nhdon:(alvadesc)"))
# lipinski_conditional
# #Summary table
# lipinski_df[,-c(1, 3, 4)] <- lapply(lipinski_df[,-c(1, 3, 4)], as.factor)
# summary(lipinski_df[-c(1)])
# all_lipinski_2 <- list(donors = as.character(nodu$cID[nodu$`nhdon:(alvadesc)` > 5]), 
#                      acceptors = as.character(nodu$cID[nodu$`nhacc:(alvadesc)` > 10]), 
#                      mw = as.character(nodu$cID[nodu$`mw:(alvadesc)` > 500]),
#                      mlogp = as.character(nodu$cID[nodu$`mlogp:(alvadesc)` > 4.15]))
# #library("ggVennDiagram")
# ggVennDiagram(all_lipinski_2, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "white", high = "yellow")


#nodu$cID[nodu$`mw:(alvadesc)` > 500 & nodu$`nhacc:(alvadesc)` > 10]
#nodu$cID[nodu$`mlogp:(alvadesc)` > 4.15 & nodu$`nhacc:(alvadesc)` > 10]
#nodu$cID[nodu$`ro5:(alvadesc)` == 1]



############3.7. Checking the RULE OF 3 (`lls_01:(alvadesc)`) ################################
#all_ro3 <- list(mw = as.character(nodu$cID[nodu$`mw:(alvadesc)` > 300]), 
#                donors = as.character(nodu$cID[nodu$`nhdon:(alvadesc)` <= 3]), 
#                acceptors = as.character(nodu$cID[nodu$`nhacc:(alvadesc)` <= 3]),
#                mlogp = as.character(nodu$cID[nodu$`mlogp:(alvadesc)` <= 3]),
#                rbn = as.character(nodu$cID[nodu$`rbn:(alvadesc)` <= 3]),
#                tpsa_no = as.character(nodu$cID[nodu$`tpsa(no):(alvadesc)` <= 60]),
#                ro3 = as.character(nodu$cID[nodu$`lls_01:(alvadesc)` != 0]))

all_ro3 <- list(mw = as.character(nodu$cID[nodu$`mw:(alvadesc)` > 300]), 
                donors = as.character(nodu$cID[nodu$`nhdon:(alvadesc)` <= 3]), 
                acceptors = as.character(nodu$cID[nodu$`nhacc:(alvadesc)` <= 3]),
                mlogp = as.character(nodu$cID[nodu$`mlogp:(alvadesc)` <= 3]),
                rbn = as.character(nodu$cID[nodu$`rbn:(alvadesc)` <= 3]),
                tpsa_tot = as.character(nodu$cID[nodu$`tpsa(tot):(alvadesc)` <= 60]),
                ro3 = as.character(nodu$cID[nodu$`lls_01:(alvadesc)` == 0]))

ggVennDiagram(all_ro3, label_alpha = 0) + 
  ggplot2::scale_fill_gradient(low = "white", high = "yellow")


# ro3_df <- nodu[, c("cID", "lls_01:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)",
#                    "nhacc:(alvadesc)", "nhdon:(alvadesc)", "rbn:(alvadesc)", 
#                    "tpsa(tot):(alvadesc)")]
# ro3_df$`lls_01:(alvadesc)` <- as.factor(ro3_df$`lls_01:(alvadesc)`)
# ro3_df$`nhdon:(alvadesc)` <- as.factor(ro3_df$`nhdon:(alvadesc)`)
# ro3_df$`nhacc:(alvadesc)` <- as.factor(ro3_df$`nhacc:(alvadesc)`)
# ro3_df$`rbn:(alvadesc)` <- as.factor(ro3_df$`rbn:(alvadesc)`)
# summary(ro3_df)
# 
# ro3_df <- ro3_df %>% group_by(`lls_01:(alvadesc)`) %>% sample_n(5)
# ro3_df
# 
# ro3_df <- nodu[, c("cID", "lls_01:(alvadesc)", "mlogp:(alvadesc)", "mw:(alvadesc)",
#                    "nhacc:(alvadesc)", "nhdon:(alvadesc)", "rbn:(alvadesc)", 
#                    "tpsa(tot):(alvadesc)")]
# ro3_df$`lls_01:(alvadesc)` <- as.factor(ro3_df$`lls_01:(alvadesc)`)
# ro3_df <- ro3_df %>% group_by(`lls_01:(alvadesc)`) %>% sample_n(5)
# 
# ro3_conditional <- tableHTML(ro3_df) %>%
#   add_css_conditional_column(color_rank_theme = "RAG", 
#                              columns = c("lls_01:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 3, 
#                              css = list('background-color', "red"),
#                              columns = c("mlogp:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 300, 
#                              css = list('background-color', "red"),
#                              columns = c("mw:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 3, 
#                              css = list('background-color', "red"),
#                              columns = c("nhacc:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 3, 
#                              css = list('background-color', "red"),
#                              columns = c("nhdon:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 3, 
#                              css = list('background-color', "red"),
#                              columns = c("rbn:(alvadesc)")) %>%
#   add_css_conditional_column(conditional = ">",
#                              value = 60, 
#                              css = list('background-color', "red"),
#                              columns = c("tpsa(tot):(alvadesc)")) 
# 
# ro3_conditional

#write.csv(ro3_df, "~/MHB/ICB/Data/HTS_transfer/data/Files_preprocessing\\ro3_df_lls.csv") 

#NOte: rule of 3 vs rule of 6
#lls_02 = lls_02 (rule of 8) 
#llls_01 (rule of 3) <-------- !!!!!!!!! use lls.1!!



############3.8. Tested compounds (NOT USED FOR MWE DATASET) #####
####################3.8.1. Checking tested compounds ####
# tested <- read_xlsx("testedcompounds.xlsx", col_types = "text")
# tested$`Library Id` <- as.factor(tested$`Library Id`)
# library(ggVennDiagram)
# 
# check_tested <- list(tested = as.character(tested$`Compound Id`), 
#                      all = as.character(nodu$cID))
# ggVennDiagram(check_tested, label_alpha = 0) + 
#   ggplot2::scale_fill_gradient(low = "light blue", high = "yellow")
# #194 compounds do not have matching cID 
# 
# cid_tested <- list(FDA = as.character(tested$`Compound Id`[tested$`Library Id` == "FDA"]), 
#                      ChemBridge = as.character(tested$`Compound Id`[tested$`Library Id` == "ChemBridge"]),
#                      ChemDiv = as.character(tested$`Compound Id`[tested$`Library Id` == "ChemDiv"]),
#                      CNS = as.character(tested$`Compound Id`[tested$`Library Id` == "CNS"]),
#                      PPI = as.character(tested$`Compound Id`[tested$`Library Id` == "PPI"]),
#                      Enamine = as.character(tested$`Compound Id`[tested$`Library Id` == "Enamine"]))
# ggVennDiagram(cid_tested, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "white ", high = "yellow")

####################3.8.2. Finding out the 194 unmatched compounds are all from FDA library #####
# FDA_sdf <- read.csv(file = "FDA_sdf.csv")
# cID_FDA <- FDA_sdf$id #compound ID, dataframe with 1 column with the compound names
# cID_FDA <- as.data.frame(cID_FDA)
#cID_FDA$cID_FDA <- as.factor(cID_FDA$cID_FDA)
# check_tested3 <- list(tested = as.character(tested$`Compound Id`), 
#                       fda = as.character(cID_FDA))
# ggVennDiagram(check_tested3, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "light blue", high = "yellow")
# 
# filter(tested_unique, grepl("^Prest*", `Compound Id`))
# filter(cID_FDA, grepl("^Acetaz*", cID_FDA))
# filter(cID_FDA, grepl("^Prest*", cID_FDA))

####################3.8.3. Finding out CNS = PPI compounds ####
# CNSwho <- list(CNS = as.character(tested$`Compound Id`[tested$`Library Id` == "CNS"]),
#                FDA = as.character(pca_data$cID[pca_data$Library == "FDA"]), 
#                ChemBridge = as.character(pca_data$cID[pca_data$Library == "ChemBridge"]),
#                ChemDiv = as.character(pca_data$cID[pca_data$Library == "ChemDiv"]),
#                PPI = as.character(pca_data$cID[pca_data$Library == "PPI"]),
#                Enamine = as.character(pca_data$cID[pca_data$Library == "Enamine"]))
# CNS_PPI <- list(CNS = as.character(tested$`Compound Id`[tested$`Library Id` == "CNS"]),
#                 PPI = as.character(pca_data$cID[pca_data$Library == "PPI"]))
# ggVennDiagram(CNS_PPI, label_alpha = 0) + ggplot2::scale_fill_gradient(low = "white ", high = "yellow")

####################3.8.4. FDA tested compounds cID mismatch ####
###### FDA_tested - tested compounds from FDA library had different cID + smiles, 
#ran OCHEM to identify them through features
# tested_FDA <- read.csv("tested_FDA.csv")
# 
# FDA194 <- excel_sheets("descriptors_FDA194.xlsx") %>% map(~read_xlsx("descriptors_FDA194.xlsx", .)) #join all tables by row (with map_dfc original column headers were not kept)
# FDA194 <- bind_cols(lapply(FDA194, as.data.frame)) #convert to dataframe 194 x 10750
# FDA194 <- cbind(tested_FDA$Compound.Id, FDA194) #194 x 10751
# colnames(FDA194) <- tolower(colnames(FDA194))
# FDA194 <- alerts_changename(FDA194)
# names(FDA194)[names(FDA194) == "tested_fda$compound.id"] <- "cID"
# 
# 
# FDA194 <- FDA194[ ,c("cID", "smiles", "mw:(alvadesc)", "alogps_logp", "alogps_logs", "nhacc:(alvadesc)", "nhdon:(alvadesc)")]
# FDA194[,-c(1:2)] <- sapply(FDA194[,-c(1:2)],as.numeric)
# 
# nodu_FDA <- nodu[nodu$library == "FDA", c("cID", "smiles", "mw:(alvadesc)", "alogps_logp", "alogps_logs", "nhacc:(alvadesc)", "nhdon:(alvadesc)")]
# nodu_FDA[,-c(1:2)] <- sapply(nodu_FDA[,-c(1:2)],as.numeric)
# 
# 
# equal <- data.frame()
# for (i in 1:194){
#   for (j in 1:1275){
#     if (FDA194$`mw:(alvadesc)`[i] == nodu_FDA$`mw:(alvadesc)`[j]){
#       if (FDA194$alogps_logp[i] == nodu_FDA$alogps_logp[j]){
#         if (FDA194$alogps_logs[i] == nodu_FDA$alogps_logs[j]){
#           if (FDA194$`nhacc:(alvadesc)`[i] == nodu_FDA$`nhacc:(alvadesc)`[j]){
#             if (FDA194$`nhdon:(alvadesc)`[i] == nodu_FDA$`nhdon:(alvadesc)`[j]){
#               equal[i, 1] <- FDA194$cID[i]
#               equal[i, 2] <- nodu_FDA$cID[j]
#             }
#           }
#         }
#       }
#     }
#   }
# } #find their cID in the nodu (FDA library) dataframe, by comparing mw, alogps and nhacc, nhdon 
# equal
# 
# #Checking smiles are equal
# equal[6,]
# FDA194$smiles[FDA194$cID == "Prestw-102"] == nodu_FDA$smiles[nodu_FDA$cID == "Amoxapine"]
# for (i in 1:194){
#   print((FDA194$smiles[FDA194$cID == as.character(equal[i,1])]) == 
#     (nodu_FDA$smiles[nodu_FDA$cID == as.character(equal[i,2])]))
# }
# equal[16,]
# FDA194$smiles[FDA194$cID == "Prestw-1105"] == nodu_FDA$smiles[nodu_FDA$cID == "Verteporfin"]
# FDA194$smiles[FDA194$cID == "Prestw-1105"] 
# nodu_FDA$smiles[nodu_FDA$cID == "Verteporfin"]
# equal[26,]
# FDA194$smiles[FDA194$cID == "Prestw-1188"] == nodu_FDA$smiles[nodu_FDA$cID == "Ziprasidone  Hydrochloride"]
# FDA194$smiles[FDA194$cID == "Prestw-1188"]
# nodu_FDA$smiles[nodu_FDA$cID == "Ziprasidone  Hydrochloride"]
# #though smiles are different they seem to refer to the same structure

#no match for #144 and #179
#Prestw-745 - is Spiramycin (incorrect smiles)
#Prestw-919 - unidentified, large compound

#> FDA194[179, ]
# smiles
# 179 CCC(C)C(N)C1=NC(CS1)C(=O)NC(CC(C)C)C(=O)NC(CCC(O)=O)C(=O)NC(C(C)CC)C(=O)NC1CCCCNC(=O)C(CC(N)=O)NC(=O)C(CC(O)=O)NC(=O)C(CC2=CN=CN2)NC(=O)C(CC2=CC=CC=C2)NC(=O)C(NC(=O)C(CCCN)NC1=O)C(C)CC
# mw:(alvadesc) alogps_logp alogps_logs nhacc:(alvadesc) nhdon:(alvadesc)
# 179          1420       -1.07       -4.77               32               20

# > FDA194[144, ]
# cID smiles mw:(alvadesc) alogps_logp alogps_logs nhacc:(alvadesc) nhdon:(alvadesc)
# 144 Prestw-745      C            16        0.92        -3.2                0                0



####################3.8.5. Conclusion for tested values ######
#' 1. Save equivalences
#write.csv(equal, "~/MHB/ICB/Data/HTS_transfer/data/Files_preprocessing\\equal_testedFDAandnoduFDA.csv") 
#' 2. Substitute FDA names in tested df (also change Prestw-745 to Spiramycin) ***new plan
#' 3. Substitute Spiramycin in nodu with new (correct) features after runnning OCHEM with real smiles
#' 4. Add Prestw-919 as new compound in nodu df *** new plan
#' 5. Create dummy variable with 0 = not tested and 1 = tested in nodu dataframe (dont forget to add Prest-919 to tested list first)
#' 6. Check how many 1s - should be 7464 = # tested compounds
#' 7. Save new nodu
#' 8. Run PCAs & UMAPS again


##### 4. Using features to classify the compounds ####

#############4.1. PCA ####
####################4.1.1. Intro ####
#' prcomp returns x, sdev and rotation
#' 
#' x = contains the PCs for drawing a graph. there are as many PCs as samples, in this case 29771
#' the first PC accounts for the most variation in the original data, 2nd PC 
#' for the second most variation. to plot 2d pca graph we usually use the first 2 PCs
#' 
#' sdev = standard deviation, we use the square of sdev to calculate how much 
#' variation in the original data each principal component accounts for, 
#' pca_var <- pca$sdev^2 
#' however usually turn it into a percentage by using 
#' percentages <- round(pca_var/sum(pca_var) * 100, 1)
#' we can plot percentages with a barplot.
#' barplot(percentages, main = "Scree Plot", xlab = "Principal Component", ylab = "Percent Variation")
#' 
#' Scree plot with ggplot2
#' pca_data <- data.frame(Sample = rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])
#' this created a dataframe with 1 column with sample ids, and two columns for the 
#' X and Y coordinates for each sample
#' ggplot(data = pca_data, aes(x = X, y = Y, label = Sample)) + geom_text() + 
#' xlab(paste("PC1 - ", percentages[1], %, sep = ")) + 
#' ylab(paste("PC2 - ", percentages[2], %, sep = ")) + 
#' theme_bw()
#' 
#' check which features have the largest effect on where samples are plotted 
#' in the PCA plot by using loadings
#' rotation = loading scores for each PC. Eg for PC1, if it accounts for most 
#' variation, will have loading scores. Features or variables that push the 
#' samples to the right side of the graph will have large positive values, 
#' varuables that push samples to the left side of the graph will have negative
#' since we are interested in both sets of variables, we can use the abs() 
#' function to sort based on the numbers magnitude independent of negative or positive
#' we can get the names of the top 10 features with the largest loading scores
#' loading_scores <- pca$rotation[,1]
#' feature_scores <- abs(loading_scores)
#' feature_scores_ranked <- sort(feature_scores, decreasing = TRUE)
#' top_10_features <- names(feature_scores_ranked[1:10])
#' pca$rotation[top_10_features, 1] #show the scores with - or + sign

require(ggplot2)
#nodu[, -c(1:4)]

#nodu[, -c(1:4, 5509, 5510)] <- sapply(nodu[, -c(1:4, 5509, 5510)],as.numeric)
nodu[, -c(1:4, 4473)] <- sapply(nodu[, -c(1:4, 4473)],as.numeric)
test <- nodu[, -c(1:4, 4473)] #removing character/factor columns: cID, smiles, library, moleculeid, inchi_key, tested
removeZeroVar <- function(df){
  df[, !sapply(df, function(x) min(x, na.rm = TRUE) == max(x, na.rm = TRUE))]
} #remove zero variance columns 
test <- removeZeroVar(test)  #31045 X 5504

#test2 <- test
####################4.1.2. PCA unscaled ####
##PCA UNSCALED
# nodu_pca <- prcomp(test) #PCA unscaled
# rotations <- as.data.frame(nodu_pca$rotation[, 1:2])
# rotations$absPC1 <- abs(rotations$PC1)
# rotations$absPC2 <- abs(rotations$PC2)
# 
# df_nodu_pca <- as.data.frame(nodu_pca$x)
# df_nodu_pca$group <- nodu$library
# df_nodu_pca$lipinski <- factor(nodu$`ro5:(alvadesc)`, levels = unique(nodu$`ro5:(alvadesc)`)) #number of lipinski failures
# df_nodu_pca$RO3 <- factor(nodu$`lls_01:(alvadesc)`, levels = unique(nodu$`lls_01:(alvadesc)`)) #rule of three
# df_nodu_pca$size <- nodu$`mw:(alvadesc)`
# 
# # labels for x and y axis
# percentage_unscaled <- round(nodu_pca$sdev / sum(nodu_pca$sdev) * 100, 2)
# percentage_unscaled <- paste(colnames(df_nodu_pca), "(", paste(as.character(percentage), "%", ")", sep=""))
# 
# ggplot(df_nodu_pca, aes(x = PC1, y = PC2, colour = "size")) +
#   geom_point() +
#   xlab(percentage_unscaled[1]) +
#   ylab(percentage_unscaled[2]) +
#   theme_minimal()
# #geom_text(data = subset(df_nodu_pca_scaled, PC1 > 100), aes(label=rownames(subset(df_nodu_pca_scaled, PC1 > 100))))
# #geom_count() # adds count data where points are identical
# ggplot(df_nodu_pca, aes(x = PC1, y = PC2)) +
#   stat_binhex()+
#   xlab(percentage[1])+
#   ylab(percentage[2])+
#   theme_minimal()



####################4.1.3. PCA scaled ####
########## PCA SCALED 
Sys.time()
nodu_pca_scaled <- prcomp(test, scale. = TRUE) #PCA scaled #takes about 1h 30 min
df_nodu_pca_scaled <- as.data.frame(nodu_pca_scaled$x) #dataframe with PCs 31046 x 6933
Sys.time()
#saveRDS(nodu_pca_scaled, file = "nodu_pca_scaled.Rds")
#nodu_pca_scaled <- readRDS(file = "nodu_pca_scaled.Rds")

#Screeplot
pca_var <- nodu_pca_scaled$sdev^2 
percentages <- round(pca_var/sum(pca_var) * 100, 1)
x <- barplot(percentages[1:20], main = "Scree Plot", xlab = "Principal Component", ylab = "Percent Variation", ylim=c(0,30))
text(x, percentages[1:20] + 2, labels = as.character(percentages[1:20]))



#PCA plot
# pca_data <- data.frame(cID = nodu$cID,
#                        Library = nodu$library, 
#                        MW = nodu$`mw:(alvadesc)`,
#                        PCA1 = nodu_pca_scaled$x[,1], 
#                        PCA2 = nodu_pca_scaled$x[,2],
#                        TWC = nodu$`twc:(alvadesc)`,
#                        RO3 = nodu$`lls_01:(alvadesc)`,
#                        Lipinski = nodu$`ro5:(alvadesc)`,
#                        tested = nodu$tested)

pca_data <- data.frame(cID = nodu$cID,
                       MW = nodu$`mw:(alvadesc)`,
                       PCA1 = nodu_pca_scaled$x[,1],
                       PCA2 = nodu_pca_scaled$x[,2],
                       TWC = nodu$`twc:(alvadesc)`)
rownames(pca_data) <- paste0(nodu$cID, "_", nodu$library) #change rownames to cID_library
#write.csv(pca_data, file = "~/MHB/ICB/Data/HTS_transfer/data/Files_preprocessing\\pca_data.csv")

ggplot(data = pca_data, aes(x = PCA1, y = PCA2)) + 
  geom_point(alpha = 0.1) + 
  xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
  ylab(paste("PC2 - ", percentages[2], "%", sep = "")) 

loading_scores <- nodu_pca_scaled$rotation[,1]
feature_scores <- abs(loading_scores)
feature_scores_ranked <- sort(feature_scores, decreasing = TRUE)
top_20_features <- names(feature_scores_ranked[1:20])
top_20_features
as.data.frame(nodu_pca_scaled$rotation[top_20_features, 1]) #show the scores with - or + sign



#Check for clustering
ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = MW > 250)) + 
  geom_point(alpha = 0.1) + 
  stat_ellipse(type = "t") +
  xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
  ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
  theme_bw()


# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = Library)) + 
#   stat_binhex() + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw()

#In aes, Change color to MW, TWC, Lipinski, RO3...
# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = RO3)) + 
#   geom_point(alpha = 0.5) + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw() +
#   scale_color_gradient(low = "red", high = "blue")
# 
# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = as.factor(RO3) == 1)) + 
#   geom_point(alpha = 0.5) + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw() +
#   scale_colour_manual(values = c("gray", "black")) 
# 
# lipinski_fail <- pca_data %>% filter(Lipinski == 1)
# false_positives <- pca_data %>% 
#   filter(cID %in% c("Pemetrexed disodium", "Reserpine", "Verteporfin", "Halofantrine hydrochloride"))
# false_negatives <- pca_data %>% 
#   filter(cID %in% c("C200-5954", "C694-0016", "D268-0425"))
# ggplot(data = pca_data, aes(x = PCA1, y = PCA2)) + 
#   geom_point(alpha = 0.5, color = "light gray") +
#   geom_point(data = lipinski_fail, 
#              aes(x = PCA1, y = PCA2), color = "red") +
#   geom_point(data = false_positives, 
#              aes(x = PCA1, y = PCA2), color = "yellow") +
#   geom_point(data = false_negatives, 
#              aes(x = PCA1, y = PCA2), color = "black") +
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw() 
# 
# 
# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = MW > 500)) + 
#   geom_point(alpha = 0.5) + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw()
# 
# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = TWC > 14)) + 
#   geom_point(alpha = 0.5) + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw()


# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = Lipinski)) + 
#   stat_binhex() + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw() 

# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID, color = as.factor(RO3))) + 
#   geom_point(alpha = 0.5) + 
#   stat_ellipse(type = "t") +
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw() +
#   scale_color_manual(values = c("blue", "green", "red", "orange", "purple", "pink", "cyan"))


# ggplot(data = pca_data, aes(x = PCA1, y = PCA2, label = cID)) + 
#   geom_point(alpha = 0.5, aes(colour = as.factor(MW > 500))) + 
#   xlab(paste("PC1 - ", percentages[1], "%", sep = "")) + 
#   ylab(paste("PC2 - ", percentages[2], "%", sep = "")) +
#   theme_bw()



####################4.1.4. PCA scaled top 500 features ####

#############4.2. UMAP ##################
# require(umap)
# Sys.time()
# nodu_umap <- umap(test) #use dataset with no zero variance and no character columns (see pca)
# Sys.time()
# df_nodu_umap <- data.frame(x = nodu_umap$layout[,1], #layout holds matrix with coordinates used to visualize the dataset
#                          y = nodu_umap$layout[,2], 
#                          cID = nodu$cID,
#                          Library = nodu$library,
#                          MW = nodu$`mw:(alvadesc)`,
#                          RO3 = nodu$`lls_01:(alvadesc)`,
#                          lipinski = factor(nodu$`ro5:(alvadesc)`, levels = unique(nodu$`ro5:(alvadesc)`)))
# 
# 
# ggplot(df_nodu_umap, aes(x, y, colour = factor(lipinski))) +
#   geom_point(alpha = 0.2) +
#   theme_minimal()
# 
# ggplot(df_nodu_umap, aes(x, y, colour = factor(Library))) +
#   geom_point(alpha = 0.5) +
#   theme_minimal()
# 
# ggplot(df_nodu_umap, aes(x, y, colour = factor(RO3))) +
#   geom_point(alpha = 0.5) +
#   theme_minimal()
# 
# 
# ggplot(df_nodu_umap, aes(x, y, colour = as.factor(RO3) == 1)) +
#   geom_point(alpha = 0.5) +
#   theme_minimal() +
#   scale_colour_manual(values = c("gray", "black")) 
# 
# 
# 
# high_mw <- df_nodu_umap %>% filter(MW > 500)
# ggplot(df_nodu_umap, aes(x, y, colour = "black")) +
#   geom_point(alpha = 0.5) +
#   theme_minimal() + 
#   geom_point(data = high_mw, aes(x, y, colour = "gray"))
# 
# 
# 
# # k-mean clustering - needs defined number of clusters - use elbow plot to determine - biased 
# 








########################################
#### PART 2: PRIMARY SCREEN DATASET ####
# Objective: build a model that is able to predict anticancer properties of compounds by using chemical features as input
# 
# 


########################################


#3. NORMALISED DATASET ####
#### 3.1. Import data ####
# screen_cc <- read.csv("~/MWE/PrimaryScreen_PI_CC_withNEN.csv",
#                       stringsAsFactors = F, na.strings=c("", " "))
screen_cc <- read.csv("~/MWE/PrimaryScreen_PI_CC_withNEN_MWE1.csv",
                      stringsAsFactors = F, na.strings=c("", " "))

screen_cc$Readout.PI <- as.integer(screen_cc$Readout.PI) #38400 x 7
screen_cc$Readout..cell.count <- as.integer(screen_cc$Readout..cell.count) 
screen_cc <- screen_cc[!is.na(screen_cc$Readout..cell.count), ] #36765 x 7

#### 3.2. Adding controls ####
## For each of the controls (DMSO = negative control, NEN + dilazep, NEN + DMSO - pos control) we
## Plot average cell count, PI cell viability and ratio PI/CC per wellplate to look for inter-plate variability/bias in screening assay)

############## 3.2.1. DMSO ####
DMSO <- screen_cc %>% group_by(Plate.Id) %>% filter(Well.Literal == "POP_1") %>% 
  summarise(DMSO_CC_mean = mean(Readout..cell.count), 
            DMSO_CC_var = var(Readout..cell.count),
            DMSO_PI_mean = mean(Readout.PI), 
            DMSO_PI_var = var(Readout.PI),
            ratio = mean(Readout.PI/Readout..cell.count),
            ratio_var = var(Readout.PI/Readout..cell.count))

# ggplot(data = DMSO, aes(x = Plate.Id, y = DMSO_CC_mean)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = DMSO_CC_mean - sqrt(DMSO_CC_var), ymax = DMSO_CC_mean + sqrt(DMSO_CC_var)))
# 
# ggplot(data = DMSO, aes(x = Plate.Id, y = DMSO_PI_mean)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = DMSO_PI_mean - sqrt(DMSO_PI_var), ymax = DMSO_PI_mean + sqrt(DMSO_PI_var)))
# 
ggplot(data = DMSO, aes(x = Plate.Id, y = ratio)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ratio - sqrt(ratio_var), ymax = ratio + sqrt(ratio_var)))

############## 3.2.3. NEN + DMSO ####
nen_dmso <- screen_cc %>% group_by(Plate.Id) %>% filter(Well.Literal == "POP_3") %>% 
  summarise(nen_dmso_CC_mean = mean(Readout..cell.count), 
            nen_dmso_CC_var = var(Readout..cell.count),
            nen_dmso_PI_mean = mean(Readout.PI), 
            nen_dmso_PI_var = var(Readout.PI),
            ratio = mean(Readout.PI/Readout..cell.count),
            ratio_var = var(Readout.PI/Readout..cell.count))

# ggplot(data = nen_dmso, aes(x = Plate.Id, y = nen_dmso_CC_mean)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = nen_dmso_CC_mean - sqrt(nen_dmso_CC_var), ymax = nen_dmso_CC_mean + sqrt(nen_dmso_CC_var)))
# 
# ggplot(data = nen_dmso, aes(x = Plate.Id, y = nen_dmso_PI_mean)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = nen_dmso_PI_mean - sqrt(nen_dmso_PI_var), ymax = nen_dmso_PI_mean + sqrt(nen_dmso_PI_var)))

ggplot(data = nen_dmso, aes(x = Plate.Id, y = ratio)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ratio - sqrt(ratio_var), ymax = ratio + sqrt(ratio_var)))


#### 3.3. Preprocessing ####
############## 3.3.1. Remove rows of controls and rows with missing values ####
screen_cc <- screen_cc[complete.cases(screen_cc),] #30360 x 7

############## 3.3.2. Density distributions of PI, CC and ratios  ####
##Compound density distribution of PI (fluorescence = citotoxicity assay), cell count and ratio PI/CC
ggplot(data = screen_cc, aes(x = Readout..cell.count, color = as.factor(Plate.Id))) +
  geom_density()
ggplot(data = screen_cc, aes(x = Readout.PI, color = as.factor(Plate.Id))) +
  geom_density()
ggplot(data = screen_cc, aes(x = Readout.PI/Readout..cell.count, color = as.factor(Plate.Id))) +
  geom_density()

#### Obtain a table with % inhibition, p_inhibition, a measurement of cell viability/citotoxicity of a compound defined as 
### (ratio compound - ratio neg control)/(ratio pos control - ratio neg control)
screen_cc$PI_CC <- screen_cc$Readout.PI/screen_cc$Readout..cell.count
screen_cc$PI_normalised <- NA
screen_cc$CC_normalised <- NA
screen_cc$PI_CC_normalised <- NA
for (i in 1:nrow(screen_cc)){
  j <- as.numeric(which(as.character(DMSO$Plate.Id) %in% as.character(screen_cc$Plate.Id[i])))
  screen_cc$PI_normalised[i] <- screen_cc$Readout.PI[i]/DMSO$DMSO_PI_mean[j]
  screen_cc$CC_normalised[i] <- screen_cc$Readout..cell.count[i]/DMSO$DMSO_CC_mean[j]
  screen_cc$PI_CC_normalised[i] <- (screen_cc$Readout.PI[i]/screen_cc$Readout..cell.count[i])/DMSO$ratio[j]
  screen_cc$p_inhibition[i] <- (screen_cc$PI_CC[i] - DMSO$ratio[j])/(nen_dmso$ratio[j] - DMSO$ratio[j])
}



# ggplot(data = screen_cc, aes(x = CC_normalised, color = as.factor(Plate.Id))) +
#   geom_boxplot()
# ggplot(data = screen_cc, aes(x = PI_normalised, color = as.factor(Plate.Id))) +
#   geom_boxplot()
ggplot(data = screen_cc, aes(x = PI_CC_normalised, color = as.factor(Plate.Id))) +
  geom_boxplot()

ggplot(data = screen_cc, aes(x = p_inhibition, color = as.factor(Plate.Id))) +
  geom_density()

#Note. We will remove outliers from regression by keeping values within mean +- 3SD 
length(which(screen_cc$p_inhibition > (mean(screen_cc$p_inhibition, na.rm = TRUE) + 3 * sd(screen_cc$p_inhibition, na.rm = TRUE))))
length(which(screen_cc$p_inhibition < (mean(screen_cc$p_inhibition, na.rm = TRUE) - 3 * sd(screen_cc$p_inhibition, na.rm = TRUE))))
#removes 170 compounds in total

############## 3.3.3. Investigating duplicates  
############## 3.3.4. Changing "Prestw" names 
############## 3.3.5.Checking compound IDs between screen_cc and features 
############## 3.3.6. Final clean-up ####
#' We will :
#' 1. delete compounds in GHEN12 plate
#' 2. keep only compounds within mean +- 3SD %inhibition value (p_inhibition)
#' 3. keep duplicates! they will be assigned the same features
#' 4. remove compounds without features
############## 3.3.7. Formatting input ####






