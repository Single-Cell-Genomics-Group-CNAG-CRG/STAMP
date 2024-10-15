wh <- function(w,h){
    options(repr.plot.width = w, repr.plot.height = h)
}
###############################################################################
# 0. Palette general

library(Polychrome)

palette_general <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(20, c("#fc6060", "#74f774", "#7c7cfc"))
  return(unname(color_palette))
}
###############################################################################

###############################################################################
# 1. Pal Lvl2 T
library(Polychrome)
  set.seed(1234) # for reproducibility
  pal_lvl2_t <- Polychrome::createPalette(10, c("#fc6060", "#74f774", "#7c7cfc"))
  names(pal_lvl2_t) <- c("CD4","CD8","NK","LowQ")
###############################################################################


pal_clines <- c('SKBR3' = '#8DD3C7', 'MCF7' = '#BEBADA', 'LnCAP' = '#FB8072')
###############################################################################
# 1. Palette for seurat cell_types_pred

library(Polychrome)


palette_cell_types <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(200, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- c("NKT cells" = color_palette[1], 
               "Naive CD4+ T cells" = color_palette[2],
               "Cytotoxic CD8+ T cells" = color_palette[3],
               "Central Memory CD4+ T cells" = color_palette[4],
               "Tr1" = color_palette[5],
               "Th" = color_palette[6],
               "Proliferative CD4+ T cells" = color_palette[7],
               "Treg" = color_palette[8],
               "Itgb1+ Memory CD4+ T cells" = color_palette[9],
               "Naive CD8+ T cells" = color_palette[10],
               "Activated Effector Memory CD8+ T cells" = color_palette[11],
               "Central Memory CD8+ T cells" = color_palette[12],
               "IFN+ Effector CD4+ T cells" = color_palette[13],
               "Cytotoxic Effector Memory CD4+ T cells" = color_palette[14],
               "IFNgR+ Memory CD8+ T cells" = color_palette[15],
               "dn T cells" = color_palette[16],
               "IFN+ Effector Memory CD8+ T cells" = color_palette[17],
               "Proliferative CD8+ T cells" = color_palette[18],
               "MAIT cells" = color_palette[19],
               "Il6st+ Naive CD4+ T cells" = color_palette[20])
  return(palette)
}
###############################################################################


###############################################################################
# 1.2 Full palette for cell types, including non T cells

library(Polychrome)


palette_cell_types_full <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(200, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- c("NKT cells" = color_palette[1], 
               "Naive CD4+ T cells" = color_palette[2],
               "Cytotoxic CD8+ T cells" = color_palette[3],
               "Central Memory CD4+ T cells" = color_palette[4],
               "Tr1" = color_palette[5],
               "Th" = color_palette[6],
               "Treg" = color_palette[8],
               "Itgb1+ Memory CD4+ T cells" = color_palette[9],
               "Naive CD8+ T cells" = color_palette[10],
               "Activated Effector Memory CD8+ T cells" = color_palette[11],
               "Central Memory CD8+ T cells" = color_palette[12],
               "Cytotoxic Effector Memory CD4+ T cells" = color_palette[14],
               "dn T cells" = color_palette[16],
               "IFN+ Effector Memory CD8+ T cells" = color_palette[17],
               "MAIT cells" = color_palette[19],
               "Il6st+ Naive CD4+ T cells" = color_palette[20],
               "Classical Monocytes Il1b+" = color_palette[21],
               "Classical Monocytes" = color_palette[22],
               "CD4+ T cell-Platelet complexes" = color_palette[23],
               "ILC" = color_palette[24],
               "Switched CD27- Memory B cells" = color_palette[25],
               "Non-Classical Monocytes" = color_palette[26],
               "IFNgR+ Memory CD8+ B cells" = color_palette[27],
               "Classical Monocytes CD163+" = color_palette[28],
               "Naive B cells" = color_palette[29],
               "Monocytes Csfr3+" = color_palette[30],
               "CD16+ high RB NK cells" = color_palette[31],
               "CD56dim CD16+ CD57- NK cells" = color_palette[32],
               "Classical Monocyte-platelet complexes" = color_palette[33],
               "Transitional B cells" = color_palette[32],
               "CD56dim CD16+ CD57+ NK cells" = color_palette[33],
               "DC4" = color_palette[34],
               "Unswitched Memory B cells" = color_palette[35],
               "XISThi B cells" = color_palette[36],
               "Intermediate Monocytes" = color_palette[37],
               "CD19hiCD11c+ Memory B cells" = color_palette[38],
               "DC2 Cd1c-" = color_palette[39],
               "Progenitor cells" = color_palette[40],
               "CD8+ T cell-Platelet complexes" = color_palette[41],
               "CD56bright NK cells" = color_palette[42]
               
  )
  return(palette)
}
###############################################################################
###############################################################################
# 1. Palette for clonal homeostasis plot

palette_homeostasis <- function(){
  
  library(Polychrome)
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(32, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- color_palette[10:15]
  
  palette <- c("Expanded (180 < X <= 300)" = "#B2182B",
               "Medium (5 < X <= 20)" = "#90EE90" ,
               "Large (20 < X <= 80)" = "#FD8D3C",
               "Large_expanded (80 < X <= 180)" = "#D6604D",
               "Small (1 < X <= 5)" = "grey")
  return(palette)
}
###############################################################################

###############################################################################
# 1.B Continuous palette

palette_continuous <- function(){
  library(Polychrome)
  set.seed(1234) # for reproducibility
  color_palette <- colorRampPalette(c("#322783","#fbcc1c"))
  return(color_palette)
}
###############################################################################

###############################################################################
# 4. Function to calculate diversity within cell subtypes in scRepertoire

subtype_diversity = function(combined,subtype){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 2. Calculate diversity in NR-antiIL6-W0 
  antiIL6_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W0 <- clonalDiversity(antiIL6_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W0 <- bind_rows(diversity_antiIL6_NR_W0,diversity_antiIL6_R_W0)
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_W0,index, value, colnames(diversity_antiIL6_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 4. Calculate diversity in NR-JAKi-W0 
  JAKi_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W0 <- clonalDiversity(JAKi_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W0 <- bind_rows(diversity_JAKi_NR_W0,diversity_JAKi_R_W0)
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_W0,index, value, colnames(diversity_JAKi_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 6. Calculate diversity in NR-antiIL6-W12 
  antiIL6_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W12 <- clonalDiversity(antiIL6_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W12 <- bind_rows(diversity_antiIL6_NR_W12,diversity_antiIL6_R_W12)
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_W12,index, value, colnames(diversity_antiIL6_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 8. Calculate diversity in NR-JAKi-W12 
  JAKi_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W12 <- clonalDiversity(JAKi_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W12 <- bind_rows(diversity_JAKi_NR_W12,diversity_JAKi_R_W12)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_W12,index, value, colnames(diversity_JAKi_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  
  
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  # Calculate p values
  statistics <- full_diversity %>%
    group_by(therapy, index, week) %>%
    wilcox_test(data =.,
                formula = value~response,
                paired = FALSE) %>%
    adjust_pvalue(method = "fdr")
  
  
  # add xy_position and y_position to correctly put the p values on the x and y axis
  statistics <- statistics %>% 
    add_x_position(x = "response", dodge = 0.8)%>%
    add_y_position(scales = "free")
  
  #adjust the number of significant digits
  statistics$p.adj <- round(statistics$p.adj, digits = 5)
  
  
  function_output <- list(full_diversity, statistics, subtype)
  names(function_output) <- c("full_diversity","statistics","subtype")
  
  return(function_output)
  
}
###############################################################################

# 4. Function to calculate diversity within cell lineages in scRepertoire

cell_lineages_diversity = function(combined,cell_lineage){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_lineages_pred == cell_lineage & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 2. Calculate diversity in NR-antiIL6-W0 
  antiIL6_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W0 <- clonalDiversity(antiIL6_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W0 <- bind_rows(diversity_antiIL6_NR_W0,diversity_antiIL6_R_W0)
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_W0,index, value, colnames(diversity_antiIL6_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 4. Calculate diversity in NR-JAKi-W0 
  JAKi_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W0 <- clonalDiversity(JAKi_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W0 <- bind_rows(diversity_JAKi_NR_W0,diversity_JAKi_R_W0)
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_W0,index, value, colnames(diversity_JAKi_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_lineages_pred == cell_lineage & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 6. Calculate diversity in NR-antiIL6-W12 
  antiIL6_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W12 <- clonalDiversity(antiIL6_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W12 <- bind_rows(diversity_antiIL6_NR_W12,diversity_antiIL6_R_W12)
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_W12,index, value, colnames(diversity_antiIL6_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 8. Calculate diversity in NR-JAKi-W12 
  JAKi_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W12 <- clonalDiversity(JAKi_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W12 <- bind_rows(diversity_JAKi_NR_W12,diversity_JAKi_R_W12)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_W12,index, value, colnames(diversity_JAKi_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  
  
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  function_output <- list(full_diversity, cell_lineage)
  names(function_output) <- c("full_diversity","cell_lineage")
  
  return(function_output)
  
}

###############################################################################
#function to calculate celltype specific diversity in Responders

subtype_R_diversity = function(combined,subtype){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  #filter for Responders
  subtyped_combined <- lapply(combined, function(df) subset(df, response == "R"))
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_R_W0,index, value, colnames(diversity_antiIL6_R_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_R_W0,index, value, colnames(diversity_JAKi_R_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_R_W12,index, value, colnames(diversity_antiIL6_R_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_R_W12,index, value, colnames(diversity_JAKi_R_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  # Calculate p values
  statistics <- full_diversity %>%
    group_by(therapy, index, response) %>%
    wilcox_test(data =.,
                formula = value~week,
                paired = FALSE) %>%
    adjust_pvalue(method = "fdr")
  
  
  # add xy_position and y_position to correctly put the p values on the x and y axis
  statistics <- statistics %>% 
    add_x_position(x = "week", dodge = 0.8)%>%
    add_y_position(scales = "free")
  
  #adjust the number of significant digits
  statistics$p.adj <- round(statistics$p.adj, digits = 5)
  
  
  function_output <- list(full_diversity, statistics, subtype)
  names(function_output) <- c("full_diversity","statistics","subtype")
  
  return(function_output)
  
}

###############################################################################
# 5. Function to plot cell type specific diversity 

plot_celltypediversity_function <- function(diversity_output){
  
  therapies <- c("JAKi","antiIL6")
  
  diversity <- diversity_output[[1]]
  statistics <- diversity_output[[2]]
  subtype <- diversity_output[[3]]
  
  plot <- list()
  for (t in therapies) {
    
    df <- diversity %>% filter(therapy == t)
    df <- df %>% filter(index == "Shannon")
    stat <- statistics %>% filter(therapy == t)
    
    
    plot[[t]] <- ggplot(df, aes(x = response, y = value, fill = week)) +
      facet_wrap(~therapy, scales = "free") +
      geom_boxplot(alpha = 0, outlier.alpha = 0, aes(color = week), lwd = 1, position = "dodge") +
      scale_color_manual(values = c("lightblue","orchid4")) +
      geom_jitter(alpha= 0.8, size = 1.5, position = position_jitterdodge()) +
      theme(axis.text = element_text(size = 12),
            axis.text.x=element_text(angle=-45, vjust = 1),
            axis.text.x.bottom = element_text(size = 30),
            text = element_text(size = 30))
  }
  
  plot_wrapped <- patchwork::wrap_plots(plot)
  return(plot_wrapped)
}

###############################################################################

###############################################################################
# 5. Function to plot cell type specific diversity in responders across timepoints

plot_celltypediversity_R_function <- function(diversity_output){
  
  therapies <- c("JAKi","antiIL6")
  
  diversity <- diversity_output[[1]]
  statistics <- diversity_output[[2]]
  subtype <- diversity_output[[3]]
  
  plot <- list()
  for (t in therapies) {
    
    df <- diversity %>% filter(therapy == t)
    stat <- statistics %>% filter(therapy == t)
    
    
    plot[[t]] <- ggplot(df, aes(x = week, y = value)) +
      facet_wrap(~therapy~index, scales = "free") +
      geom_poi(aes(color = week), position = position_jitterdodge(0.1)) +
      geom_boxplot(alpha = 0.5, aes(fill = week), position = "dodge") +
      scale_color_manual(values = c("darkred","darkorchid4")) +
      scale_fill_manual(values = c("darkred", "darkorchid4")) +
      theme_pubr() +
      theme(axis.text = element_text(size = 12),
            axis.text.x=element_text(angle=-45, vjust = 1)) +
      stat_pvalue_manual(stat,label = "p.adj") +
      labs(title = t, subtitle = subtype) +
      theme(text = element_text(size = 20))  +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
  }
  
  plot_wrapped <- patchwork::wrap_plots(plot)
  return(plot_wrapped)
}

###############################################################################
# 6. Function to plot Symphony results

plotBasic = function(umap_labels,                  # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = NULL,                 # Plot title
                     color.by = NULL,              # metadata column name for coloring
                     facet.by = NULL,              # (optional) metadata column name for faceting
                     color.mapping = NULL,         # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  
  p = umap_labels %>%
    dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 1, stroke = 0.4, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

###############################################################################


###############################################################################
# 9. Function to export table summary for quality control

table_summary_qc <- function(seurat_obj, title){
  df <- list()
  
  libraries <- unique(seurat_obj$gem_id)
  for(library in libraries){
    seurat_fun <- subset(seurat_obj, gem_id == library)
    
    df[[library]] <- data.frame(
      "Number of cells" = length(Cells(seurat_fun)),
      "median umi per cell" = median(seurat_fun$nCount_RNA),
      "median gene per cell" = median(seurat_fun$nFeature_RNA),
      "mitochondrial fraction" = mean(seurat_fun$mt_pct))
    
    df[[library]]$gem_id <- library
    
  }
  
  df <- bind_rows(df)
  df <- cbind(df[,5],df[,-5])
  df <- rename(df, "gem_id" = "df[, 5]")
  tab <- gt(df, groupname_col = group_vars(df)) 
  tab <- tab_header(tab,title = title)
  
  return(tab)
}

# 11. Function to plot homeostasis plot with scRepertoire function 

homeostasis_plot <- function(seurat_obj){
  
  seurat <- seurat_obj
  
  #subset by therapy & response
  seurat_R_antiIL6 <-  subset(seurat, response == "R" & therapy == "antiIL6")
  seurat_NR_antiIL6 <- subset(seurat, response == "NR" & therapy == "antiIL6")
  
  seurat_R_JAKi <-  subset(seurat, response == "R" & therapy == "JAKi")
  seurat_NR_JAKi <- subset(seurat, response == "NR" & therapy == "JAKi")
  
  # Put in a list
  homeostasis_list <- list(seurat_NR_antiIL6,seurat_R_antiIL6,seurat_NR_JAKi,seurat_R_JAKi)
  names(homeostasis_list) <- c("seurat_NR_antiIL6","seurat_R_antiIL6","seurat_NR_JAKi","seurat_R_JAKi")
  
  
  #plot
  homeostasis_plot_list <- list()
  homeostasis_tables <- list()
  for(name in names(homeostasis_list)){
    df <- homeostasis_list[[name]]
    p <- occupiedscRepertoire(df, x.axis = "ident", 
                              facet.by = "week",
                              label = FALSE) +
      labs(title = name, x = "Cell types", y = "single cells") +
      coord_flip() +
      theme(text = element_text(size = 15)) +
      scale_fill_manual(values = palette)
    
    
    homeostasis_plot_list[[name]] <- p
    
    t <- occupiedscRepertoire(df, x.axis = "ident", 
                              facet.by = "week",
                              label = FALSE, 
                              exportTable = TRUE)
    homeostasis_tables[[name]] <- t
  }
  
  # remove redundant legend in plot[1]
  homeostasis_plot_list[[1]] <- homeostasis_plot_list[[1]] +
    guides(fill = "none")
  #wrap plots
  antiIL6 <- wrap_plots(homeostasis_plot_list[1:2], ncol = 2) +
    plot_annotation(title = "Clonal Homeostasis") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", text= element_text(size = 15)) 
  
  JAKi <- wrap_plots(homeostasis_plot_list[3:4], ncol = 2) +
    plot_annotation(title = "Clonal Homeostasis") &
    guides(fill = "none")
  
  
  # Plot
  list_plot <- list(antiIL6, JAKi)
  plot <- patchwork::wrap_plots(list_plot, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  output <- list()
  output[[1]] <- plot
  output[[2]] <- homeostasis_tables
  return(output)
}


##### FUNCTION TO PRE PROCESS A SUBSET OF DATA REMOVING THE TCR GENES FROM THE CLUSTERING ######

process_sub <- function(t_seu, resolution){

## Select TCR genes to exclude from downstream analysis

# For TCR Alpha Chain (TRA) Genes
TRAV <- grep("^TRAV", rownames(t_seu), value = TRUE)
TRAJ <- grep("^TRAJ", rownames(t_seu), value = TRUE)
TRAC <- grep("^TRAC", rownames(t_seu), value = TRUE)

# For TCR Beta Chain (TRB) Genes
TRBV <- grep("^TRBV", rownames(t_seu), value = TRUE)
TRBD <- grep("^TRBD", rownames(t_seu), value = TRUE)
TRBJ <- grep("^TRBJ", rownames(t_seu), value = TRUE)
TRBC <- grep("^TRBC", rownames(t_seu), value = TRUE)

# For TCR Gamma Chain (TRG) Genes
TRGV <- grep("^TRGV", rownames(t_seu), value = TRUE)
TRGJ <- grep("^TRGJ", rownames(t_seu), value = TRUE)
TRGC <- grep("^TRGC", rownames(t_seu), value = TRUE)

# For TCR Delta Chain (TRD) Genes
TRDV <- grep("^TRDV", rownames(t_seu), value = TRUE)
TRDD <- grep("^TRDD", rownames(t_seu), value = TRUE)
TRDJ <- grep("^TRDJ", rownames(t_seu), value = TRUE)
TRDC <- grep("^TRDC", rownames(t_seu), value = TRUE)

TCR_genes <- c(TRAV, TRAJ, TRAC, TRBV, TRBD, TRBJ, TRBC, TRGV, TRGJ, TRGC, TRDV, TRDD, TRDJ, TRDC)

# Normalize data and identify highly variable genes
t_seu <- t_seu %>%
  NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures(nfeatures = 2000)

## Filter the HVGs
hvg_names <- VariableFeatures(t_seu)
hvg_names <- setdiff(hvg_names, TCR_genes)
VariableFeatures(t_seu) <- hvg_names


t_seu <- t_seu %>%
  ScaleData() %>%
  RunPCA()

## Integrate

t_seu <- RunHarmony(t_seu, group.by.vars = "gem_id")
t_seu <- RunUMAP(t_seu, dims = 1:35, reduction = "harmony")

t_seu <- t_seu %>%
    FindNeighbors(dims = 1:35, reduction = "harmony") %>%
    FindClusters(resolution = resolution)


return(t_seu)
}
###############################################################################

###############################################################################
'%!in%' <- function(x, table) {
  !(x %in% table)
}
###############################################################################
					   
###############################################################################
generate_plots <- function(tmp_seurat, assay_name, feature_name) {
	# Ensure the DefaultAssay is set to the one you're interested in
	DefaultAssay(tmp_seurat) <- assay_name
	
	# Adjusted to directly use the provided feature_name for plotting
	plot_list <- list()
	plot_list[['Pre']] <- FeaturePlot(tmp_seurat[, tmp_seurat$timepoint == "Pre"], features = feature_name, split.by = "response",
			   min.cutoff = "q05", max.cutoff = "q95", keep.scale = "all") &
	scale_colour_viridis_c() &
	theme(legend.position = "right") &
	plot_annotation(title = "Pre",
	theme = theme(plot.title = element_text(size = 30)))
	
	plot_list[['Post']] <- FeaturePlot(tmp_seurat[, tmp_seurat$timepoint == "Post"], features = feature_name, split.by = "response",
				min.cutoff = "q05", max.cutoff = "q95", keep.scale = "all") &
	scale_colour_viridis_c() &
	theme(legend.position = "right") &
	plot_annotation(title = "Post",
	theme = theme(plot.title = element_text(size = 30)))
	
	return(plot_list)
}
###############################################################################

###############################################################################
# Adjusted function to generate unique identifiers ensuring all are unique
generate_unique_ids <- function(n, seed = 123) {
library(stringr)
  potential_ids <- expand.grid(letters, letters, sprintf("%02d", 0:99), stringsAsFactors = FALSE)
  potential_ids <- apply(potential_ids, 1, paste, collapse = "")
  potential_ids <- sample(potential_ids, n) # Randomly sample without replacement
  
  if (length(unique(potential_ids)) < n) {
    stop("Not enough unique IDs can be generated with the current method.")
  }
  
  return(potential_ids)
}
###############################################################################
					   
###############################################################################
addAssayDataToMetadata <- function(seu, assayName) {
  # Extract assay data
  assayData <- GetAssayData(seu, assay = assayName, slot = "data")
  
  # Convert to dataframe and transpose
  assayData <- as.data.frame(t(assayData))
  
  # Remove the now redundant barcode column
  assayData$barcode <- NULL
  
  # Add the modified assay data to the Seurat object's metadata
  seu <- AddMetaData(object = seu, metadata = assayData)
  
  return(seu)
}
###############################################################################

###############################################################################
library(ggplot2)
library(dplyr)
library(glue)
library(patchwork)

create_plots <- function(spe, variable, markers, pal) {
  # Plot reduced dimensions
  gg_um <- plotReducedDim(spe, "UMAP", colour_by = variable, text_by = variable, point_size = 1, scattermore = T, rasterise = F) +
    scale_color_manual(values = pal) +
    theme(legend.position = "none")
  
  # Create a data frame for the bar plots
  cells <- as.data.frame(table(spe[[variable]])) %>%
    dplyr::rename(Cluster = Var1) %>%
    mutate(pct = round(Freq / sum(Freq), 3))
  
  # Bar plot 1
  gg_barplot1 <- ggplot(cells, aes(x = Cluster, y = Freq, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    theme_minimal() +
    theme(text = element_text(size = 20), axis.ticks.x = element_blank(), legend.position = "none") +
    coord_flip()
  
  # Bar plot 2
  gg_barplot2 <- ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    geom_text(aes(label = paste0(pct * 100, "%")), position = position_stack(vjust = 0.5), size = 3) +
    theme_minimal() +
    scale_y_reverse() +
    labs(x = "", y = "Percentage",
         title = glue("N = {sum(cells$Freq)} cells")) +
    theme(axis.text.x = element_blank(), text = element_text(size = 20), axis.ticks.x = element_blank(), legend.position = "right")
  
  # Dotplot
  gg_dots <- plotDots(spe, features = intersect(markers, rownames(spe)), group = variable, scale =TRUE, center = TRUE) +
 coord_flip() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

  # Combine plots
  A <- wrap_plots(gg_um, gg_barplot1, gg_barplot2, ncol = 3) + plot_layout(widths = c(2, 2, 1))
  plot <- wrap_plots(A, gg_dots, nrow = 2) + plot_annotation(tag_levels = "A")
  
  return(plot)
}

# Example usage:
# spe <- your_spe_object
# pal <- your_color_palette
# variable <- "anno_lvl1"
# result_plot <- create_plots(spe, variable, pal)
# print(result_plot)
###############################################################################


###############################################################################

###############################################################################
library(ggplot2)
library(ggrepel)
library(glue)

# Define the function to create the plot
plot_hvg <- function(dec.var, sub) {
  # Convert dec.var to a data frame for ggplot2
  df <- as.data.frame(dec.var@listData)
  df$gene <- dec.var@rownames
  df <- df[sample(rownames(df)),]
  
  # Ensure hvg is a factor
  df$hvg <- as.factor(df$hvg)
  n_hvg <- n_distinct(df$gene[df$hvg == 'yes'])
  pct_hvg <- round(n_hvg / nrow(df) * 100, 2)
  
  # Map colors
  colors <- setNames(c("grey", "red"), levels(df$hvg)) # Ensure levels match "yes" and "no"
  
  # Plot with ggplot2
  p <- ggplot(df, aes(x = mean, y = total)) +
    geom_point(aes(color = hvg), shape = 16, size = 0.5) +
    scale_color_manual(values = colors) +
    labs(x = "Mean expression", y = "Variance", subtitle = paste0("Subset: ", sub),
         title = glue("HVGs, n = {n_hvg} / {pct_hvg}%")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(color = "black", size = 18),
      plot.subtitle = element_text(color = "black", size = 14),
      axis.title = element_text(color = "black", size = 14),
      axis.text = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      legend.title = element_text(color = "black", size = 14)
    )
  
  # Add the curve (assuming trend function is stored in metadata slot)
  trend_function <- metadata(dec.var)$trend
  p <- p + stat_function(fun = trend_function, color = "dodgerblue", size = 0.5)
  
  # Identify and label the top 20 points on the y-axis
  top10 <- df %>% arrange(desc(total)) %>% head(20)
  p <- p + geom_text_repel(data = top10, aes(label = gene), 
                           vjust = -1, color = "black", size = 3,
                           nudge_y = 0.05, # Adjust to nudge labels upwards
                           box.padding = 0.5, # Increase the padding around each label
                           point.padding = 0.5, # Increase the padding around each point
                           force = 2, 
                           segment.size = 0.2, segment.color = "darkgrey"
                           )
  
  return(p)
}
###############################################################################
library(parallel)
num_cores <- detectCores()
print(num_cores)
bp <- BiocParallel::MulticoreParam(num_cores)
###############################################################################

###############################################################################
library(ggplot2)
library(glue)

plot_density <- function(df, var, group_var, pal, title, thr) {
  # Ensure the variables exist in the data frame
  if (!var %in% names(df)) {
    stop(glue("Variable '{var}' not found in the data frame."))
  }
  if (!group_var %in% names(df)) {
    stop(glue("Group variable '{group_var}' not found in the data frame."))
  }
  
  # Filter the data frame
  df_filtered <- df[df[[var]] < thr, ]
  
  # Calculate mean values
  mean_values <- sapply(unique(df[[group_var]]), function(x) {
    round(mean(df[[var]][df[[group_var]] == x], na.rm = TRUE), 0)
  })
  
  # Create the plot
  gg_feats <- ggplot(df_filtered, aes_string(x = var)) + 
    geom_density(aes_string(color = group_var), key_glyph = draw_key_abline) +
    theme_minimal() +
    labs(x = var, y = "Density", title = title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = pal,
                       labels = setNames(glue("{unique(df[[group_var]])} (Mean={mean_values})"), unique(df[[group_var]])))
  
  return(gg_feats)
}

###############################################################################

##############################################################################
pal_lvl1 <- c('B' = '#FB5F5C', 'T' = '#77F873', 'Myeloid' = '#7978F7')
###############################################################################

###############################################################################
feats_cosmx <- c("CD3E","CD3D","CD247","CD4","CD8A","CD8B","NCAM1","NKG7",'KLRB1','KLRC1','KLRD1','KLRF1','KLRK1',"CD19","CD79A","CD79B","MS4A1","MZB1","IL6","CD1C","CD14","FCGR3A/B","ITGAX")
feats_xenium <- c("CD3E","CD3D","CD247","CD4","CD8A","CD8B","NCAM1","NKG7",'CMKLR1','KLRB1','KLRC1','KLRD1','KLRF1','KLRK1',"CD19","CD79A","CD79B","MS4A1","MZB1","IL6","CD1C","CD14","FCGR3A")


###############################################################################
library(ggplot2)
library(dplyr)
library(glue)
library(patchwork)
library(scater)
create_plots2 <- function(spe, variable, markers) {
  # Plot UMAP 
  df <- as.data.frame(colData(spe))
  um <- reducedDim(spe,"UMAP")
  df <- cbind(df,um)
  df <- df[sample(rownames(df)),]
  
  pal <- palette_general()
  length(pal) <- length(unique(df[[variable]]))
  names(pal) <- unique(df[[variable]])
  
  gg_um <- ggplot(df, aes(x = UMAP1, y= UMAP2, color = !!sym(variable))) +
    geom_point(shape = 16, size = 0.1) +
    scale_color_manual(values = pal) +
    theme_bw() +
    theme(panel.background = element_blank(), panel.grid = element_blank())
  
  # Create a data frame for the bar plots
  cells <- as.data.frame(table(spe[[variable]])) %>%
    dplyr::rename(Cluster = Var1) %>%
    mutate(pct = round(Freq / sum(Freq), 3))
  
  # Bar plot 2
  gg_barplot2 <- ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    geom_text(aes(label = paste0(pct * 100, "%")), position = position_stack(vjust = 0.5), size = 3) +
    theme_minimal() +
    scale_y_reverse() +
    coord_flip() +
    theme(legend.position = "bottom")
  labs(x = "", y = "Percentage",
       title = glue("N = {sum(cells$Freq)} cells")) +
    theme(axis.text.x = element_blank(), text = element_text(size = 20), axis.ticks.x = element_blank(), legend.position = "right")
  
  # Dotplot
  gg_dots <- plotDots(spe, features = intersect(markers, rownames(spe)), group = variable, scale =TRUE, center = TRUE) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  # Combine plots
  A <- wrap_plots(gg_um, gg_barplot2, ncol = 2) + plot_layout(widths = c(1,2.5))
  plot <- wrap_plots(A, gg_dots, nrow = 2) + plot_annotation(tag_levels = "A")
  
  return(plot)
}

###############################################################################

###############################################################################
library(ggplot2)
library(dplyr)
library(glue)
library(patchwork)

create_plots3 <- function(spe, variable, markers, pal) {
  # Plot reduced dimensions
  gg_um <- plotReducedDim(spe, "UMAP", colour_by = variable, point_size = 1, scattermore = T, rasterise = F) +
    scale_color_manual(values = pal) +
    theme(legend.position = "none")
  
  # Create a data frame for the bar plots
  cells <- as.data.frame(table(spe[[variable]])) %>%
    dplyr::rename(Cluster = Var1) %>%
    mutate(pct = round(Freq / sum(Freq), 3))
  
  # Bar plot 2
  gg_barplot2 <- ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    geom_text(aes(label = paste0(pct * 100, "%")), position = position_stack(vjust = 0.5), size = 3) +
    theme_minimal() +
    scale_y_reverse() +
    coord_flip() +
    theme(legend.position = "bottom")
  labs(x = "", y = "Percentage",
       title = glue("N = {sum(cells$Freq)} cells")) +
    theme(axis.text.x = element_blank(), text = element_text(size = 20), axis.ticks.x = element_blank(), legend.position = "right")
  
  
  # Combine plots
  A <- wrap_plots(gg_um, gg_barplot2, ncol = 2) + plot_layout(widths = c(1,2.5)) + plot_annotation(tag_levels = "A")

  return(A)
}
###############################################################################

mcf7_markers = c("TFF1","SPTSSB","AREG","MDK","PSMD6","ACKR3","EEF1A2","USP32","GATA3")
skbr3_markers = c("KRT7","DHRS2","ERBB2","LGALS3BP","SUSD2","S100A9","S100A6","CLU","VTCN1")
lncap_markers = c("UGT2B17", "UGT2B28", "UGT2B11", "KLK4", "KLK3", "KLK2", "ALDH3B2", "FAAH", "DDC")

###############################################################################
scientific_10 <- function(x) {
  parse(text=gsub("e\\+?0?", " %*% 10^", scales::scientific_format()(x)))
}
###############################################################################
extr_thr <- function(df){
  return <- list()
  lower <- round(attr(df, 'thresholds')['lower'],0)
  higher <- round(attr(df, 'thresholds')['higher'],0)
  return(c(lower,higher))
}

# Same as extr_thr, but to use when setting batch in isOutlier
extr_thr_batch <- function(df,batch){
  return <- list()
  lower <- round(attr(df, 'thresholds')['lower',batch],0)
  higher <- round(attr(df, 'thresholds')['higher',batch],0)
  return(c(lower,higher))
}

###############################################################################
pal_s13a <- c(
  ctrl_4h_r2        = "#1F78B4",
  LPS_4h_r2         = "#33A02C",
  aCD3aCD28_4h_r2   = "#E31A1C",
  ctrl_4h_r1        = "#A6CEE3",
  LPS_4h_r1         = "#B2DF8A",
  aCD3aCD28_4h_r1   = "#FB9A99",
  ctrl_24h_r2       = "#FF7F00",
  LPS_24h_r2        = "#6A3D9A",
  aCD3aCD28_24h_r2  = "#B15928",
  ctrl_24h_r1       = "#FDBF6F",
  LPS_24h_r1        = "#CAB2D6",
  aCD3aCD28_24h_r1  = "#FFFF99"
)
###############################################################################