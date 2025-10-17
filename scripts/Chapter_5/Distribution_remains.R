#1. Load libraries and data ----
library(vcd)
library(tidyverse)
library(ggplot2)
library(here)
library(RColorBrewer)

phy_dataset <- read.csv(
  here("data", "Chapter 5","Database_Vicofertile_Nogarole.csv"),
  na.strings = c("", "NA", "N/A")
) #general data set 

phy_distribution <- read.csv(
  here("data", "Chapter 5","phy_distribution.csv"),
  na.strings = c("", "NA", "N/A")
) #subset of the general data set inly including information about the position of the remains in the urns 

#customised palette
physical_palette <- c("#e14444","#3cf065","#5050f3","#a0e8b2","#3c7f49","#f0f050","#5cb85c","darkgrey","lightgrey")
cb_palette <-c("#e14444","darkgrey","navyblue","#E0E0E0")
lb_palette <-c("#3cf065","darkgrey","navyblue","#E0E0E0")
sb_palette <-c("#5050f3","darkgrey","navyblue","#E0E0E0")
ul_palette <-c("#a0e8b2","darkgrey","navyblue","#E0E0E0")
ll_palette <-c("#3c7f49","darkgrey","navyblue","#E0E0E0")
t_palette <-c("#f0f050","darkgrey","navyblue","#E0E0E0")

#2. Transforming characters in factors ----
phy_dataset <- phy_dataset %>%
  mutate_if(is.character, as.factor) %>%
  mutate(element_1= factor(element_1, levels=c("CB","LB","SB","UL","LL","T","L","C","I")))%>%
  mutate(section = factor(section, levels = c("upper", "central", "lower")))

#3. Calculating percentages ----
phy_dataset <- phy_dataset %>%
  mutate(total_percentage = (total_element / total_urn) * 100)

phy_dataset <- phy_dataset%>%
  mutate(section_percentage = (total_element/total_section)*100)  

phy_dataset <- phy_dataset%>%
  mutate(half_A_B_percentage = (total_element/total_half_A_B)*100) 

phy_dataset <- phy_dataset%>%
  mutate(half_C_D_percentage = (total_element/total_half_C_D)*100)

#4. Plotting distribution of elements per urn by sites ----
#VICOFERTILE
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Vicofertiele_Phy_Dist_Tot.jpeg"), #save as jpeg
     width = 700, height = 500)

ggplot(data = phy_dataset %>%
         filter(site=="VFT") %>%
         filter(burial %in% c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32")) %>%
         filter(!is.na(element_1)) %>%
         mutate(burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32"))),
       aes(x = burial, y = total_percentage, fill = element_1)) +
  geom_col() +
  ggtitle("Vicofertile Physical elements total distribution") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2)) +
  labs(x = "Burial", y = "%") +
  scale_fill_manual(name = NULL, 
                    values = physical_palette,
                    breaks = c("CB", "LB", "SB", "UL", "LL", "T", "L", "C", "I"),
                    labels = c("(CB) Cranial Bones", "(LB) Long Bones", "(SB) Spongy Bone", 
                               "(UL) Upper Limbs", "(LL) Lower Limbs", "(T) Trunk", "(L) Limbs", 
                               "(C) Cortical bone", "(I) Indeterminate")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  coord_cartesian(ylim = c(0, 100))
dev.off()

#NOGAROLE
jpeg(file = here("figures","Chapter 5","Distribution_remains", "Nogarole_Phy_Dist_Tot.jpeg"), #save as jpeg
     width = 700, height = 500)

ggplot(data = phy_dataset %>%
         filter(site=="NOGR") %>%
         filter(burial %in% c("26n", "33(US359)", "33(US360)", "36", "56")) %>%
         filter(!is.na(element_1)) %>%
         mutate(burial = factor(burial, levels = c("26n", "33(US359)", "33(US360)", "36", "56"))),
       aes(x = burial, y = total_percentage, fill = element_1)) +
  geom_col() +
  ggtitle("Nogarole Physical elements total distribution") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2)) +
  labs(x = "Burial", y = "%") +
  scale_fill_manual(name = NULL, 
                    values = physical_palette,
                    breaks = c("CB", "LB", "SB", "UL", "LL", "T", "L", "C", "I"),
                    labels = c("(CB) Cranial Bones", "(LB) Long Bones", "(SB) Spongy Bone", 
                               "(UL) Upper Limbs", "(LL) Lower Limbs", "(T) Trunk", "(L) Limbs", 
                               "(C) Cortical bone", "(I) Indeterminate")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  coord_cartesian(ylim = c(0, 100))
dev.off()

#5. Plotting anatomical categories by sides and sites ----
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Elements_side_VFT_NOGR.jpeg"), #save as jpeg
     width = 600, height = 400)

ggplot(
  phy_dataset %>%
    filter(
      !is.na(side),
      !is.infinite(total_element),
      side != "I"
    ) %>%
    group_by(element_1, side, site) %>%
    summarise(
      summed_percentage = sum(total_element, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(summed_percentage > 0),
  aes(x = side, y = summed_percentage, fill = site)  
) +
  geom_boxplot(color = "black", alpha = 0.7, position = position_dodge(width = 0.8)) +  
  stat_summary(fun = median, geom = "point", aes(group = site), shape = 22, size = 3, color = "black") +
  facet_wrap(~element_1, scales = "fixed") +  
  scale_fill_manual(values = c("VFT" = "#FF5733", "NOGR" = "#3498DB")) + 
  scale_y_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, by = 25)
  ) +
  labs(
    title = "Distribution of Anatomical Categories by Side and Site",
    x = "Side",
    y = "Count",
    fill = "Site"  
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  )
dev.off()


#6. Plotting vertical distribution of elements by site ----
#VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Sections_total_VFT.jpeg"), #save as jpeg
     width = 700, height = 500)

  ggplot(phy_dataset %>% 
           filter(!is.na(element_1), site != "NOGR") %>%
           mutate(burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32"))),
         aes(x = section, y = section_percentage, fill = element_1)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Vicofertile",
    y = "%",
    x = "Section"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))
  dev.off()


#NOGR   
  jpeg(file = here("figures", "Chapter 5","Distribution_remains","Sections_total_NOGR.jpeg"), #save as jpeg
       width = 700, height = 500)
  
ggplot(phy_dataset %>% filter(!is.na(element_1), site != "VFT"), 
       aes(x = section, y = section_percentage, fill = element_1)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Nogarole Rocca",  
    y = "%",                
    x = "Section"           
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))  
dev.off()


#7. Mosaic Plot - Burial ~ element_section ----

#Organise data
v_phy_dist_merged <- phy_distribution %>%
  mutate(section_abbr = case_when(
    section == "upper" ~ "u",
    section == "central" ~ "c",
    section == "lower" ~ "l",
    TRUE ~ section  
  )) %>%
  mutate(element_section = factor(paste(element_1, section_abbr, sep = "."))) %>%
  dplyr::select(-section_abbr)  

#Replace some labels
v_phy_dist_merged <- v_phy_dist_merged %>%
  mutate(burial = gsub("\\(US359\\)", "(359)", burial),
         burial = gsub("\\(US360\\)", "(360)", burial))

#Order factors
v_phy_dist_merged <- v_phy_dist_merged %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.u", "CB.c", "CB.l", 
    "LB.u", "LB.c", "LB.l", 
    "SB.u", "SB.c", "SB.l", 
    "UL.u", "UL.c", "UL.l", 
    "LL.u", "LL.c", "LL.l", 
    "T.u", "T.c", "T.l"
  )))

# Convert 'burial' to a factor with specified order directly
v_phy_dist_merged$burial <- factor(v_phy_dist_merged$burial, 
                                   levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                              "26v", "29", "30", "32", "26n", "33(359)", "33(360)", "36", "56"))


#Create subsets for VFT and NOGR
v_phy_dist_vft <- subset(v_phy_dist_merged, site == "VFT")
v_phy_dist_nogr <- subset(v_phy_dist_merged, site == "NOGR")

#Order the 'burial' variable in the desired order for each subset
#VFT subset
v_phy_dist_vft$burial <- factor(v_phy_dist_vft$burial, 
                                levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                           "26v", "29", "30", "32"))

v_phy_dist_vft <- v_phy_dist_vft %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.u", "CB.c", "CB.l", 
    "LB.u", "LB.c", "LB.l", 
    "SB.u", "SB.c", "SB.l", 
    "UL.u", "UL.c", "UL.l", 
    "LL.u", "LL.c", "LL.l", 
    "T.u", "T.c", "T.l"
  )))

#NOGR subset
v_phy_dist_nogr$burial <- factor(v_phy_dist_nogr$burial, 
                                 levels = c( "26n", "33(359)", "33(360)", "36", "56"))

v_phy_dist_nogr <- v_phy_dist_nogr %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.u", "CB.c", "CB.l", 
    "LB.u", "LB.c", "LB.l", 
    "SB.u", "SB.c", "SB.l", 
    "UL.u", "UL.c", "UL.l", 
    "LL.u", "LL.c", "LL.l", 
    "T.u", "T.c", "T.l"
  )))%>% 
  filter(!element_section %in% c("LB.u", "LB.c")) %>%
  droplevels()


#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_vertical_VFT.jpeg"), #save as jpeg
     width = 700, height = 500)

vcd::mosaic(~burial + element_section, 
            data = v_phy_dist_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_vertical_NOGR.jpeg"), #save as jpeg
     width = 700, height = 500)

vcd::mosaic(~burial + element_section, 
            data = v_phy_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0.2, 0.9),  
              rot_labels = c(60, 0),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Convert data frame to contingency table for Chi-Square extraction
v_phy_table_vft <- table(v_phy_dist_vft$burial, v_phy_dist_vft$element_section)
v_phy_table_nogr <- table(v_phy_dist_nogr$burial, v_phy_dist_nogr$element_section)

#Perform chi-square test
v_chi_res_nogr <- chisq.test(v_phy_table_nogr, simulate.p.value = TRUE, B = 10000) 
v_chi_res_vft <- chisq.test(v_phy_table_vft, simulate.p.value = TRUE, B = 10000) 

#Extract standardized Pearson residuals
v_residuals_df_vft <- as.data.frame(v_chi_res_vft$stdres)  
v_residuals_df_nogr <- as.data.frame(v_chi_res_nogr$stdres)  

#Heatmap of extracted residuals
#VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","element_vertical_VFT_heat.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(v_residuals_df_vft, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  labs(
    title = "Vicofertile",
    x = "Burial",
    y = "Element_section",
    fill = "Residuals"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)  
  )
dev.off()

#heatmap of extracted residuals
jpeg(file = here("figures", "Chapter 5","Distribution_remains","element_vertical_NOGR_heat.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(v_residuals_df_nogr, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  labs(title = "Nogarole Rocca",
       x = "Burial",
       y = "Element_section",
       fill = "Residuals") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)  # Center the title
  )
dev.off()

#8. Plotting horizontal distribution of elements by site - Half A/B ----
#VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Elements_halfAB_VFT.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(
  phy_dataset %>%
    filter(
      !is.na(half_A_B),
      !is.na(half_A_B_percentage),
      site != "NOGR"
    ) %>%
    mutate(burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32"))),
  aes(x = half_A_B, y = half_A_B_percentage, fill = element_1)
) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Vicofertile",
    y = "%",
    x = "Half"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))
dev.off()

#NOGR
jpeg(file = here("figures","Chapter 5", "Distribution_remains","Elements_halfAB_NOGR.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(
  phy_dataset %>%
    filter(
      !is.na(half_A_B),
      !is.na(half_A_B_percentage),
      site != "VFT"
    ),
  aes(x = half_A_B, y = half_A_B_percentage, fill = element_1)
) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Nogarole Rocca",
    y = "%",
    x = "Half"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))
dev.off()

#9. Plotting horizontal distribution of elements by site - Half C/D ----
# VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Elements_halfCD_VFT.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(
  phy_dataset %>%
    filter(
      !is.na(half_C_D),
      !is.na(half_C_D_percentage),
      site != "NOGR"
    ) %>%
    mutate(burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32"))),
  aes(x = half_C_D, y = half_C_D_percentage, fill = element_1)
) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Vicofertile",
    y = "%",
    x = "Half"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))
dev.off()

#NOGR
jpeg(file = here("figures", "Chapter 5","Distribution_remains","Elements_halfCD_NOGR.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(
  phy_dataset %>%
    filter(
      !is.na(half_C_D),
      !is.na(half_C_D_percentage),
      site != "VFT"
    ),
  aes(x = half_C_D, y = half_C_D_percentage, fill = element_1)
) +
  geom_bar(stat = "identity") +
  facet_wrap(~ burial) +
  theme_minimal() +
  scale_fill_manual(values = physical_palette) +
  labs(
    title = "Nogarole Rocca",
    y = "%",
    x = "Half"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL))
dev.off()


#10. Mosaic Plot - Burial ~ element_section (HALF A/B) ----
#remove NA
h_phy_dist_merged <- v_phy_dist_merged %>%
  filter(!is.na(half_A_B) & !is.na(half_C_D))

#mix element_1 and halves-ha/b
h_phy_dist_merged <- h_phy_dist_merged %>%
  mutate(section_abbr = case_when(
    half_A_B == "half_A" ~ "ha",
    half_A_B == "half_B" ~ "hb",
    TRUE ~ half_A_B  
  )) %>%
  mutate(element_section = factor(paste(element_1, section_abbr, sep = "."))) %>%
  dplyr::select(-section_abbr)  

#Create subsets for VFT and NOGR
hab_phy_dist_vft <- subset(h_phy_dist_merged, site == "VFT")
hab_phy_dist_nogr <- subset(h_phy_dist_merged, site == "NOGR")

#Order the 'burial' variable in the desired order for each subset
# VFT subset
hab_phy_dist_vft$burial <- factor(hab_phy_dist_vft$burial, 
                                  levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                             "26v", "29", "30", "32"))

hab_phy_dist_vft <- hab_phy_dist_vft %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.ha", "CB.hb",
    "LB.ha", "LB.hb",  
    "SB.ha", "SB.hb",  
    "UL.ha", "UL.hb", 
    "LL.ha", "LL.hb", 
    "T.ha", "T.hb"
  )))
# NOGR subset
hab_phy_dist_nogr$burial <- factor(hab_phy_dist_nogr$burial, 
                                   levels = c( "26n", "33(359)", "33(360)", "36", "56"))

hab_phy_dist_nogr <- hab_phy_dist_nogr %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.ha", "CB.hb",
    "LB.ha", "LB.hb",  
    "SB.ha", "SB.hb",  
    "UL.ha", "UL.hb", 
    "LL.ha", "LL.hb", 
    "T.ha", "T.hb"
  )))%>% 
  filter(!element_section %in% c("LB.ha")) %>%
  droplevels()


#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_halfAB_VFT_mos.jpeg"), #save as jpeg
     width = 500, height = 400)

vcd::mosaic(~burial + element_section, 
            data = hab_phy_dist_vft, 
            shade = TRUE,
            main= "Vicofertile",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_halfAB_NOGR_mos.jpeg"), #save as jpeg
     width = 500, height = 400)

vcd::mosaic(~burial + element_section, 
            data = hab_phy_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#11. Mosaic Plot - Burial ~ element_section (HALF C/D) ----
#remove NA
h_phy_dist_merged <- v_phy_dist_merged %>%
  filter(!is.na(half_A_B) & !is.na(half_C_D))

#mix element_1 and halves-hc/d
h_phy_dist_merged <- h_phy_dist_merged %>%
  mutate(section_abbr = case_when(
    half_C_D == "half_C" ~ "hc",
    half_C_D == "half_D" ~ "hd",
    TRUE ~ half_C_D  
  )) %>%
  mutate(element_section = factor(paste(element_1, section_abbr, sep = "."))) %>%
  dplyr::select(-section_abbr)  

#Create subsets for VFT and NOGR
hcd_phy_dist_vft <- subset(h_phy_dist_merged, site == "VFT")
hcd_phy_dist_nogr <- subset(h_phy_dist_merged, site == "NOGR")

#Order the 'burial' variable in the desired order for each subset
# VFT subset
hcd_phy_dist_vft$burial <- factor(hcd_phy_dist_vft$burial, 
                                  levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                             "26v", "29", "30", "32"))

hcd_phy_dist_vft <- hcd_phy_dist_vft %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.hc", "CB.hd",
    "LB.hc", "LB.hd",  
    "SB.hc", "SB.hd",  
    "UL.hc", "UL.hd", 
    "LL.hc", "LL.hd", 
    "T.hc", "T.hd"
  )))


# NOGR subset
hcd_phy_dist_nogr$burial <- factor(hcd_phy_dist_nogr$burial, 
                                   levels = c( "26n", "33(359)", "33(360)", "36", "56"))

hcd_phy_dist_nogr <- hcd_phy_dist_nogr %>%
  mutate(element_section = factor(element_section, levels = c(
    "CB.hc", "CB.hd",
    "LB.hc", "LB.hd",  
    "SB.hc", "SB.hd",  
    "UL.hc", "UL.hd", 
    "LL.hc", "LL.hd", 
    "T.hc", "T.hd"
  )))%>% 
  filter(!element_section %in% c("LB.hd")) %>%
  droplevels()

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_halfCD_VFT_mos.jpeg"), #save as jpeg
     width = 500, height = 400)

vcd::mosaic(~burial + element_section, 
            data = hcd_phy_dist_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Distribution_remains","elements_halfCD_NOGR_mos.jpeg"), #save as jpeg
     width = 500, height = 400)

vcd::mosaic(~burial + element_section, 
            data = hcd_phy_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()