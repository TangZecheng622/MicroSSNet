devtools::load_all("C:\\Users\\86135\\Desktop\\MicroSSNet")
library(MicroSSNet) 
library(dplyr)



load(file = "C:\\Users\\86135\\Documents\\WeChat Files\\wxid_o0ez2x4irqn722\\FileStorage\\File\\2024-11\\cov_depth.RData")
group_pig <- read.csv("C:\\Users\\86135\\Desktop\\Eurp_meta_use.csv",row.names = 1) 

depth <- depth_filter[rownames(cov_filter),colnames(cov_filter)]
cover <- cov_filter

# filtered_depth <- ifelse(cover > 0.4, depth, NA)
# filtered_depth <- as.data.frame(filtered_depth)

filtered_depth <- depth %>%
  mutate(across(everything(), ~ ifelse(cover[[cur_column()]] > 0.4, ., NA)))
filtered_depth[is.na(filtered_depth)] <- 0


europe_pig <- filtered_depth[,group_pig$sample[group_pig$Group == "Europe"]]

PCA_draw(europe_pig, group_pig, "Region", save = FALSE, showType = "centroid", mode = "PCOA", offset = FALSE)
PCA_draw(europe_pig, group_pig, "Region", save = FALSE, showType = "centroid", mode = "PCA", offset = FALSE)
# europe_pig1 <- read.csv("C:\\Users\\86135\\Desktop\\Eurp_depth_filter.csv",row.names = 1)


# group_pig_t <- group_pig[93:107,]

setwd("C:\\Users\\86135\\Desktop\\zangpig\\all")
ssn_pipeline(table1 = europe_pig,group_df = group_pig,vscol1 = "Region",ssn_method = "ssPCC",pvl_threshold = 0.5,r_threshold = 0,save = TRUE,property = FALSE,scale = TRUE,control = "all",log = TRUE)
setwd("C:\\Users\\86135\\Desktop\\zangpig\\pergroup/")
ssn_pipeline(europe_pig,group_df = group_pig,vscol1 = "Region",ssn_method = "ssPCC",pvl_threshold = 0.5,r_threshold = 0,save = TRUE,property = FALSE,scale = TRUE,control = "pergroup")
cor <- data.table::fread("C:\\Users\\86135\\Desktop\\zangpig\\pergroup\\SSN\\sspcc\\Scaled_Total_-1_1scaled.tsv")

table1 <- europe_pig
group_df <- group_pig
vscol1 = "Region"
vscol2 = vscol1
ssn_method = "ssPCC"
pvl_threshold = 0.5
r_threshold = 0
save = TRUE
property = FALSE
scale = TRUE
control = "all"
log = TRUE

######test1 
setwd("C:\\Users\\86135\\Desktop\\ssn-test")
table1 <- read.csv("otu_table_09_14.csv",row.names = 1)
group_df <- read.csv("GROUP2.csv")
group_df <- group_df[-(1:24),]
table2 <- table1[,-(1:24)]
check_result <- process_group_and_table(table = table2, group_df = group_df, vscol = "Warm")
table3 <- check_result$table
group_df2 <- check_result$group_df
group_list <- check_result$group_list
group_list2 <- unique(group_df$Warm)
table_i_list <- list()
vscol1 <- "Warm"
vscol2 <- vscol1
pvl_threshold <- 0.5
for (i in group_list) {
  i <- as.character(i)
  # i = group_list[[1]]
  
  # 根据group_df过滤table1中的列
  selected_table <- table2[, group_df2[,1][group_df2[[vscol1]] == i], drop = FALSE]
  
  for (j in group_list2) {
    # j <- group_list2[[2]]
    j <- as.character(j)
    cols_to_select <- group_df2[, 1][group_df2[[vscol2]] == j]
    valid_cols <- cols_to_select[cols_to_select %in% colnames(selected_table)]
    if(length(valid_cols) == 0){
      next
    }else{
      table_i <- selected_table[, valid_cols]
      table_i <- filter_OTU2(table_i, Pre = pvl_threshold)
      n <- paste(i,"_",j)
      table_i_list[[n]] <- table_i
    }
  }
  merged_table <- table_i_list %>%
    purrr::map(~ .x %>% tibble::rownames_to_column(var = "row_name")) %>% # 将行名转换为一列
    purrr::reduce(dplyr::full_join, by = "row_name") %>%                  # 按"row_name"列合并
    replace(is.na(.), 0) %>%                                              # 填充缺失值为0
    tibble::column_to_rownames(var = "row_name") 
}
table3 <- table2[rownames(merged_table),]

table1p <- log1p(table2)
aggregation_netpipeline(table1 = table1p,group_df = group_df2,vscol1 = "Warm",pvl.threshold = 0.5,r.threshold = 0.75,p.threshold = 0.05,method = "spearman",rm.p.list = seq(0, 0.6, by = 0.1),
                        output_dir = "./test1r8p01_spearman/",ncpus = 10,step = 50)
######test2
table1i <- read.csv("table2iseed8.csv",row.names = 1)
sow1i <- read.csv("sow2iseed8.csv",row.names = 1)
group_dfi <- group_df[group_df$Sample_name %in% colnames(table1i),]

check_result <- process_group_and_table(table = table1, group_df = group_dfi, vscol = "Warm")
table1i <- check_result$table
group_dfi <- check_result$group_df
group_list <- check_result$group_list
group_list2 <- unique(group_dfi$Warm)
table_i_list <- list()
for (i in group_list) {
  i <- as.character(i)
  # i = group_list[[1]]
  
  # 根据group_df过滤table1中的列
  selected_table <- table2[, group_dfi[,1][group_dfi[[vscol1]] == i], drop = FALSE]
  
  for (j in group_list2) {
    # j <- group_list2[[2]]
    j <- as.character(j)
    cols_to_select <- group_dfi[, 1][group_dfi[[vscol2]] == j]
    valid_cols <- cols_to_select[cols_to_select %in% colnames(selected_table)]
    if(length(valid_cols) == 0){
      next
    }else{
      table_i <- selected_table[, valid_cols]
      table_i <- filter_OTU2(table_i, Pre = pvl_threshold)
      n <- paste(i,"_",j)
      table_i_list[[n]] <- table_i
    }
  }
  merged_table <- table_i_list %>%
    purrr::map(~ .x %>% tibble::rownames_to_column(var = "row_name")) %>% # 将行名转换为一列
    purrr::reduce(dplyr::full_join, by = "row_name") %>%                  # 按"row_name"列合并
    replace(is.na(.), 0) %>%                                              # 填充缺失值为0
    tibble::column_to_rownames(var = "row_name") 
}
table4_16 <- merged_table

PCA_draw(table4_16, group_dfi, "Warm", save = FALSE, showType = "centroid", mode = "PCA", offset = FALSE)
PCA_draw(sow1i, group_dfi, "Warm", save = FALSE, showType = "centroid", mode = "PCA", offset = FALSE)

PCA_draw(table4_16, group_dfi, "Warm", save = FALSE, showType = "centroid", mode = "PCOA", offset = FALSE)
PCA_draw(sow1i, group_dfi, "Warm", save = FALSE, showType = "centroid", mode = "PCOA", offset = FALSE)


######test3
cor <- data.table::fread("C:\\Users\\86135\\Desktop\\test3\\SSN\\sspcc\\Scaled_Total_-1_1scaled.tsv")
cor <- as.data.frame(cor)
cor[, -c(1, 2)][abs(cor[, -c(1, 2)]) < 0.3] <- 0
cor_top <- cor[rowSums(abs(cor[, -c(1, 2)])) > 0, ]
SOW <- Calculate_sum_of_weights(network = cor_top)
SOW1 <- SOW[,-c(1:24)]

table_abu_limma <- Differential_nodes(network = table3,sample_metadata = group_df2,offset = FALSE,vscol = "Warm",save = FALSE,plot = TRUE)
table_sow_limma <- Differential_nodes(network = SOW1,sample_metadata = group_df2,offset = FALSE,vscol = "Warm",save = FALSE,plot = TRUE)

abu_limma_result <- table_abu_limma$LimmaResult
abu_limma_result$group <- "Abundance"
sow_limma_result <- table_sow_limma$LimmaResult
sow_limma_result$group <- "Node_Weight"

compare_abu_sow <- rbind(abu_limma_result,sow_limma_result)



########test4
##auc roc

