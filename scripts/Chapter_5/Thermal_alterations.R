#1. Load libraries, data, and customise palette ----
library(vcd)
library(tidyverse)
library(ggplot2)
library(here)
library(RColorBrewer)

phy_dataset <- read.csv(
  here("data", "Chapter 5","Database_Vicofertile_Nogarole.csv"),
  na.strings = c("", "NA", "N/A")
) #general data set 

phases_fdsl <- read.csv(
  here("data", "Chapter 5","Phases_fract_stai_dim.csv"),
  na.strings = c("", "NA", "N/A")
) #subset specific for phases, fractures, staining, and dimensional categories

#customised palette
physical_palette <- c("#e14444","#3cf065","#5050f3","#a0e8b2","#3c7f49","#f0f050","#5cb85c","darkgrey","lightgrey")
cb_palette <-c("#e14444","darkgrey","navyblue","#E0E0E0")
lb_palette <-c("#3cf065","darkgrey","navyblue","#E0E0E0")
sb_palette <-c("#5050f3","darkgrey","navyblue","#E0E0E0")
ul_palette <-c("#a0e8b2","darkgrey","navyblue","#E0E0E0")
ll_palette <-c("#3c7f49","darkgrey","navyblue","#E0E0E0")
t_palette <-c("#f0f050","darkgrey","navyblue","#E0E0E0")


#2. Plotting the phases in % for individuals between necropolises ----

jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phases_burials.jpeg"), 
     width = 700, height = 400) #save as jpeg


ggplot(
  phy_dataset %>%
    mutate(burial = case_when(
      site == "VFT" ~ factor(burial, levels = c(7, 8, 12, 15, 16,17, 19, 20, "26v", 29, 30, 32)),
      TRUE ~ as.factor(burial)  
    )) %>%
    dplyr::select(site, burial, P1, P2, P3) %>%
    pivot_longer(cols = c(P1, P2, P3), names_to = "P_category", values_to = "count") %>%
    group_by(site, burial, P_category) %>%
    summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(site, burial) %>%
    mutate(percentage = (total_count / sum(total_count)) * 100) %>%
    filter(!is.na(percentage)),  # Data transformation ends here
  aes(x = burial, y = percentage, fill = P_category)
) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site, scales = "free_x") +  
  theme_minimal() +
  labs(title = "Prevalence of P1, P2, and P3 as Percentage of Each Burial",
       x = "Burial",
       y = "Percentage",
       fill = "P Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

#3. Mosaic Plot - Burial ~ Phase_section ----

#dataframe subset for Phases and elements
phases_subset <- phy_dataset %>%
  dplyr::select(site, burial, section, half_A_B, half_C_D, element_1, P1, P2, P3)

#Replace some labels
phases_subset <- phases_subset %>%
  mutate(burial = gsub("\\(US359\\)", "(359)", burial),
         burial = gsub("\\(US360\\)", "(360)", burial))

phases_subset <- phases_subset %>% mutate(across(P1:P3, ~replace_na(., 0)))

#Compute total duplication needed
phases_subset <- phases_subset %>%
  mutate(rep_count = P1 + P2 + P3) %>%  
  uncount(rep_count)  

#Subsetting and organising data  
phases_subset <- phases_subset %>%
  mutate(phase = case_when(
    P1 > 0 ~ "P1",
    P2 > 0 ~ "P2",
    P3 > 0 ~ "P3",
    TRUE ~ NA_character_  
  ),
  phase_section = paste0(phase, 
                         case_when(
                           section == "upper" ~ ".u",
                           section == "central" ~ ".c",
                           section == "lower" ~ ".l",
                           TRUE ~ NA_character_  
                         )))


#Create subsets for VFT and NOGR
v_phases_dist_vft <- subset(phases_subset, site == "VFT")
v_phases_dist_nogr <- subset(phases_subset, site == "NOGR")

#Order the variables in the desired order for each subset
#VFT subset
v_phases_dist_vft$burial <- factor(v_phases_dist_vft$burial, 
                                   levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                              "26v", "29", "30", "32"))

v_phases_dist_vft <- v_phases_dist_vft %>%
  mutate(phase_section = factor(phase_section, levels = c(
    "P1.u", "P1.c", "P1.l", 
    "P2.u", "P2.c", "P2.l", 
    "P3.u", "P3.c", "P3.l"
  )))

#NOGR subset
v_phases_dist_nogr$burial <- factor(v_phases_dist_nogr$burial, 
                                    levels = c( "26n", "33(359)", "33(360)", "36", "56"))

v_phases_dist_nogr <- v_phases_dist_nogr %>%
  mutate(phase_section = factor(phase_section, levels = c(
    "P1.u", "P1.c", "P1.l", 
    "P2.u", "P2.c", "P2.l", 
    "P3.u", "P3.c", "P3.l"
  )))

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_section_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~burial + phase_section, 
            data = v_phases_dist_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))

dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_section_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~burial + phase_section, 
            data = v_phases_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0.2, 0.9),  # Increase first value to move x-axis labels further away
              rot_labels = c(60, 0),  
              gp_labels = gpar(fontsize = 8)))

dev.off()

#4. Plotting the phases in % for elements between necropolises ----

jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phases_element.jpeg"), #save as jpeg
     width = 700, height = 400)

ggplot(phy_dataset %>%
         dplyr::select(site, element_1, P1, P2, P3) %>%
         pivot_longer(cols = c(P1, P2, P3), names_to = "P_category", values_to = "count") %>%
         group_by(site, element_1, P_category) %>%
         summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
         group_by(site, element_1) %>%
         mutate(percentage = (total_count / sum(total_count)) * 100) %>%
         filter(!is.na(percentage)) %>%
         mutate(element_1 = factor(element_1, levels = c("CB", "LB", "SB", "UL", "LL", "L", "T", "C"))), 
       aes(x = element_1, y = percentage, fill = P_category)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site) +
  theme_minimal() +
  labs(title = "Prevalence of P1, P2, and P3 as Percentage of Each Element",
       x = "Element",
       y = "Percentage",
       fill = "P Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

#5. Plotting phases ~ sides in counts ----
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phases_sides.jpeg"), #save as jpeg
     width = 500, height = 400)

ggplot(phy_dataset %>%
         dplyr::select(site, side, P1, P2, P3) %>%
         pivot_longer(cols = c(P1, P2, P3), names_to = "P_category", values_to = "count") %>%
         filter(side != "I") %>%  # Exclude "I"
         group_by(site, side, P_category) %>%
         summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop"),  # Sum counts
       aes(x = side, y = total_count, fill = P_category)) +  # Use total_count instead of percentage
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site) +
  theme_minimal() +
  labs(title = "Prevalence of P1, P2, and P3 as Counts of Sided Elements",
       x = "Side",
       y = "Total Count",
       fill = "P Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

#6. Mosaic Plot - Burial ~ Phase_element ----

#Subsetting and organising data
phases_elem_subset <- phases_subset %>%
  filter(!element_1 %in% c("L", "C")) %>%  
  mutate(phase = case_when(
    P1 > 0 ~ "P1",
    P2 > 0 ~ "P2",
    P3 > 0 ~ "P3",
    TRUE ~ NA_character_  
  ),
  phase_element = paste0(phase, ".", element_1))  

#Create subsets for VFT and NOGR
phases_elem_vft <- subset(phases_elem_subset, site == "VFT")
phases_elem_nogr <- subset(phases_elem_subset, site == "NOGR")

#Order the 'burial' variable in the desired order for each subset
# VFT subset
phases_elem_vft$burial <- factor(phases_elem_vft$burial, 
                                 levels = c("7", "8", "12", "15", "16", "17", "19", "20", 
                                            "26v", "29", "30", "32"))

phases_elem_vft <- phases_elem_vft %>%
  mutate(phase_element = factor(phase_element, levels = c(
    "P1.CB", "P1.LB", "P1.SB","P1.UL","P1.LL","P1.T", 
    "P2.CB", "P2.LB", "P2.SB","P2.UL","P2.LL","P2.T", 
    "P3.CB", "P3.LB", "P3.SB","P3.UL","P3.LL","P3.T" 
  )))

# NOGR subset
phases_elem_nogr$burial <- factor(phases_elem_nogr$burial, 
                                  levels = c( "26n", "33(359)", "33(360)", "36", "56"))

phases_elem_nogr <- phases_elem_nogr %>%
  mutate(phase_element = factor(phase_element, levels = c(
    "P1.CB", "P1.LB", "P1.SB","P1.UL","P1.LL","P1.T", 
    "P2.CB", "P2.LB", "P2.SB","P2.UL","P2.LL","P2.T", 
    "P3.CB", "P3.LB", "P3.SB","P3.UL","P3.LL","P3.T"
  ))) %>% 
  filter(!phase_element %in% c("P1.UL", "P1.LL", "P1.LB", "P2.LB")) %>% #filtering out zero columns
  droplevels()

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_element_VFT.jpeg"), #save as jpeg
     width = 800, height = 600)

vcd::mosaic(~burial + phase_element, 
            data = phases_elem_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(60, 30),  
              gp_labels = gpar(fontsize = 8)
            ))

dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_element_NOGR.jpeg"), #save as jpeg
     width = 800, height = 600)

vcd::mosaic(~burial + phase_element, 
            data = phases_elem_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0.2, 0.9),  
              rot_labels = c(60, 0),  
              gp_labels = gpar(fontsize = 8)
            ))

dev.off()

#Convert data frame to contingency table for Chi-Square extraction
phases_elem_table_vft <- table(phases_elem_vft$burial,phases_elem_vft$phase_element)
phases_elem_table_nogr <- table(phases_elem_nogr$burial, phases_elem_nogr$phase_element)

#Perform chi-square test
cr_phases_elem_nogr <- chisq.test(phases_elem_table_nogr, simulate.p.value = TRUE, B = 10000) 
cr_phases_elem_vft <- chisq.test(phases_elem_table_vft, simulate.p.value = TRUE, B = 10000) 

#Extract standardized Pearson residuals
rdf_phases_elem_vft <- as.data.frame(cr_phases_elem_vft$stdres)  
rdf_phases_elem_nogr <- as.data.frame(cr_phases_elem_nogr$stdres) 

#Heatmap of extracted residuals VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phases_element_VFT_heat.jpeg"), #save as jpeg
     width = 500, height = 600)

ggplot(rdf_phases_elem_vft, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  labs(title = "Vicofertile",
       x = "Burial",
       y = "Phase_element",
       fill = "Residuals") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)  
  )

dev.off()

#Heatmap of extracted residuals NOGR
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phases_element_NOGR_heat.jpeg"), #save as jpeg
     width = 500, height = 600)

ggplot(rdf_phases_elem_nogr, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  labs(title = "Nogarole Rocca",
       x = "Burial",
       y = "Phase_element",
       fill = "Residuals") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5) 
  )

dev.off()

#7. Plotting phases ~ fractures and deformations in counts ----

jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Fractures_phases.jpeg"), #save as jpeg
     width = 600, height = 400)

ggplot(phases_fdsl %>%
         pivot_longer(cols = c("mosaic", "longitudinal", "transversal", "curved", 
                               "bulls_eye", "deformation", "twisting"), 
                      names_to = "trait", values_to = "count"),
       aes(x = trait, y = count, fill = phase)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site) +  # Separate plots for each site
  theme_minimal() +
  labs(title = "Fracture Trait Counts by Phase",
       x = "Fracture Type",
       y = "Count",
       fill = "Phase") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

#8. Mosaic Plot - Fractures_deformations ~ Phases ----

#Subsetting and organising data
phases_fr_def <- phases_fdsl %>%
  select(phase, site, mosaic, longitudinal, transversal, curved, 
         bulls_eye, deformation, twisting) %>%
  pivot_longer(cols = c(mosaic, longitudinal, transversal, curved, 
                        bulls_eye, deformation, twisting), 
               names_to = "fractures_deformations", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(original_count = count) %>%  
  uncount(count) %>%
  rename(count = original_count)

#Create subsets for VFT and NOGR
phases_fr_def_vft <- subset(phases_fr_def, site == "VFT")
phases_fr_def_nogr <- subset(phases_fr_def, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures","Chapter 5", "Thermal_alterations","Phase_fractures_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~fractures_deformations + phase, 
            data = phases_fr_def_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              offset_labels = c(0, 1.5),  
              offset_varnames = c(0, 2),  
              rot_labels = c(0, 0),  
              gp_labels = gpar(fontsize = 8)
            ))

dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures","Chapter 5", "Thermal_alterations","Phase_fractures_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~fractures_deformations + phase, 
            data = phases_fr_def_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0, 1.5),  
              offset_varnames = c(0, 2),  
              rot_labels = c(0, 0),  
              gp_labels = gpar(fontsize = 8)
            ))

dev.off()


#9. Mosaic Plot - Concretion_staining ~ Phases ----

#Subsetting and organising data
phases_stain <- phases_fdsl %>%
  select(site,phase,concretion,red_brown,yellow, green, black) %>%
  pivot_longer(cols = c(concretion,red_brown,yellow, green, black), 
               names_to = "concretion_staining", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(original_count = count) %>%  
  uncount(count) %>%
  rename(count = original_count)

#Create subsets for VFT and NOGR
phases_stain_vft <- subset(phases_stain, site == "VFT")
phases_stain_nogr <- subset(phases_stain, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_staining_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~concretion_staining + phase, 
            data = phases_stain_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              offset_labels = c(0, 1.5),  
              offset_varnames = c(0, 2),  
              rot_labels = c(0, 0),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Phase_staining_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~concretion_staining + phase, 
            data = phases_stain_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0, 1.5),  
              offset_varnames = c(0, 2),  
              rot_labels = c(90, 0),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#10. Mosaic Plot - Fractures_deformations ~ element ----

#Subsetting and organising data
data_fra_def_elem <- phy_dataset %>%
  select(site, section, element_1, mosaic, longitudinal, transversal, bulls_eye, curved, deformation, twisting) %>%
  filter(!element_1 %in% c("L", "C")) %>%  
  pivot_longer(cols = c(mosaic, longitudinal, transversal, bulls_eye, curved, deformation, twisting), 
               names_to = "fractures_deformations", values_to = "count") %>%
  filter(count > 0) %>%  
  uncount(weights = count, .remove = FALSE) %>%
  mutate(
    section = factor(section, levels = c("upper", "central", "lower")), 
    fractures_deformations = factor(fractures_deformations, 
                                    levels = c("mosaic", "longitudinal", "transversal", "bulls_eye", "curved", "deformation", "twisting")),  # Reorder fracture types
    element_1 = factor(element_1, levels = c("CB", "LB", "SB", "UL", "LL", "T"))  
  ) 

#Create subsets for VFT and NOGR
data_fra_def_elem_vft <- subset(data_fra_def_elem, site == "VFT")
data_fra_def_elem_nogr <- subset(data_fra_def_elem, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","fractures_element_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~fractures_deformations + element_1, 
            data = data_fra_def_elem_vft, 
            shade = TRUE,
            main = "Vicofertile",
            set_varnames = c(element_1 = "element"), 
            labeling_args = list(
              rot_labels = c(0, 0), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures","Chapter 5","Thermal_alterations", "fractures_element_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~fractures_deformations + element_1, 
            data = data_fra_def_elem_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            set_varnames = c(element_1 = "element"),  
            labeling_args = list(
              rot_labels = c(0, 0), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#11. Plotting staining in counts between necropolises ----
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","Staining_phases.jpeg"), #save as jpeg
     width = 700, height = 400)

ggplot(phases_fdsl %>%
         pivot_longer(cols = c("concretion", "red_brown", "yellow", "green", 
                               "black"), 
                      names_to = "trait", values_to = "count"),
       aes(x = trait, y = count, fill = phase)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site) +  # Separate plots for each site
  theme_minimal() +
  labs(title = "Concretion and staining Counts by Phase",
       x = "Staining",
       y = "Count",
       fill = "Phase") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

 
#12. Mosaic Plot - Staining ~ section ----

#Subsetting and organising data
data_staining_dist <- phy_dataset %>%
  rename(black = manganese) %>%  
  select(site, section, concretion, red_brown, yellow, green, black) %>%
  pivot_longer(cols = c(concretion, red_brown, yellow, green, black), 
               names_to = "staining", values_to = "count") %>%
  filter(count > 0) %>%
  uncount(weights = count, .remove = FALSE)%>%
  mutate(section = factor(section, levels = c("upper", "central", "lower"))) 

#Create subsets for VFT and NOGR
data_staining_dist_vft <- subset(data_staining_dist, site == "VFT")
data_staining_dist_nogr <- subset(data_staining_dist, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","staining_section_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~staining + section, 
            data = data_staining_dist_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(0, 0), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for VFT
jpeg(file = here("figures","Chapter 5", "Thermal_alterations","staining_section_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~staining + section, 
            data = data_staining_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              rot_labels = c(0, 0), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#13. Plot green staining prevalence by element ----

#VFT
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","green_staining_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

(ggplot(data = phy_dataset %>%
          filter(site=="VFT") %>%
          filter(burial %in% c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32")) %>%
          filter(!is.na(element_1)) %>%
          mutate(burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32"))),
        aes(x = burial, y = green, fill = element_1)) +
   geom_col() +
   ggtitle("Vicofertile") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5, vjust = 2)) +
   labs(x = "Burial", y = "Nr stained elements") +
   scale_fill_manual(name = NULL, 
                     values = physical_palette,
                     breaks = c("CB", "LB", "SB", "UL", "LL", "T", "L", "C", "I"),
                     labels = c("(CB) Cranial Bones", "(LB) Long Bones", "(SB) Spongy Bone", 
                                "(UL) Upper Limbs", "(LL) Lower Limbs", "(T) Trunk", "(L) Limbs", 
                                "(C) Cortical bone", "(I) Indeterminate")) +
   scale_y_continuous(breaks = seq(0, 20, by = 2)) +
   coord_cartesian(ylim = c(0, 20)))
dev.off()

#NOGR
jpeg(file = here("figures", "Chapter 5","Thermal_alterations","green_staining_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

  (ggplot(data = phy_dataset %>%
            filter(site=="NOGR") %>%
            filter(burial %in% c("26n", "33(US359)", "33(US360)", "36", "56")) %>%
            filter(!is.na(element_1)) %>%
            mutate(burial = factor(burial, levels = c("26n", "33(US359)", "33(US360)", "36", "56"))),
          aes(x = burial, y = green, fill = element_1)) +
     geom_col() +
     ggtitle("Nogarole Rocca") +
     theme_minimal() +
     theme(plot.title = element_text(hjust = 0.5, vjust = 2)) +
     labs(x = "Burial", y = "Nr stained elements") +
     scale_fill_manual(name = NULL, 
                       values = physical_palette,
                       breaks = c("CB", "LB", "SB", "UL", "LL", "T", "L", "C", "I"),
                       labels = c("(CB) Cranial Bones", "(LB) Long Bones", "(SB) Spongy Bone", 
                                  "(UL) Upper Limbs", "(LL) Lower Limbs", "(T) Trunk", "(L) Limbs", 
                                  "(C) Cortical bone", "(I) Indeterminate")) +
     scale_y_continuous(breaks = seq(0, 20, by = 2)) +
     coord_cartesian(ylim = c(0, 20))+
     theme(legend.position = "none"))
dev.off()
