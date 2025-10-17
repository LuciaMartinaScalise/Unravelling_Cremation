#1. Load libraries and data ----
library(vcd)
library(tidyverse)
library(ggplot2)
library(here)

phy_dataset <- read.csv(
  here("data", "Chapter 5","Database_Vicofertile_Nogarole.csv"),
  na.strings = c("", "NA", "N/A")
) #general data set 

phases_fdsl <- read.csv(
  here("data","Chapter 5", "Phases_fract_stai_dim.csv"),
  na.strings = c("", "NA", "N/A")
) #subset specific for phases, fractures, staining, and dimensional categories

frag_data_burial<- read.csv(
  here("data","Chapter 5", "Database_Fragmentation.csv"),
  na.strings = c("", "NA", "N/A")
) #subset for fragmentation index calculation


#2. Plotting phases ~ maximum_length in counts ----
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight", "Sizes_phases.jpeg"), #save as jpeg
     width = 700, height = 400)

ggplot(phases_fdsl %>%
         pivot_longer(cols = c("cm_0_1",	"cm_1_2",	"cm_2_3",	"cm_3_4",	"cm_4_5",	"over_cm_5"), 
                      names_to = "trait", values_to = "count"),
       aes(x = trait, y = count, fill = phase)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ site) +  # Separate plots for each site
  theme_minimal() +
  labs(title = "Number of elements in each lenght class by Phase",
       x = "Staining",
       y = "Count",
       fill = "Phase") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("P1" = "black", "P2" = "grey", "P3" = "white"))

dev.off()

#3. Mosaic Plot - Maximum_lenght ~ phase ----

#Subsetting and organising data
phases_lenght <- phases_fdsl %>%
  select(site,phase,cm_0_1,	cm_1_2,	cm_2_3,	cm_3_4,	cm_4_5,	over_cm_5) %>%
  pivot_longer(cols = c(cm_0_1,	cm_1_2,	cm_2_3,	cm_3_4,	cm_4_5,	over_cm_5), 
               names_to = "maximum_lenght", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(original_count = count) %>%  # Store original count before expanding
  uncount(count) %>%
  rename(count = original_count)

#Create subsets for VFT and NOGR
phases_lenght_vft <- subset(phases_lenght, site == "VFT")
phases_lenght_nogr <- subset(phases_lenght, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","size_phase_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~maximum_lenght + phase, 
            data = phases_lenght_vft, 
            shade = TRUE,
            main ="Vicofertile",
            labeling_args = list(
              offset_labels = c(0, 1.5),  # Moves category labels
              offset_varnames = c(0, 2),  # Moves the variable name (y-axis label) further away
              rot_labels = c(0, 0),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","size_phase_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~maximum_lenght + phase, 
            data = phases_lenght_nogr , 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              offset_labels = c(0, 1.5),  # Moves category labels
              offset_varnames = c(0, 2),  # Moves the variable name (y-axis label) further away
              rot_labels = c(0, 0),  
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#4. Mosaic Plot - Maximum_lenght ~ section ----

#Subsetting and organising data
data_lenght_dist <- phy_dataset %>%
  select(site, section, cm_0_1,	cm_1_2,	cm_2_3,	cm_3_4,	cm_4_5,	over_cm_5) %>%
  pivot_longer(cols = c(cm_0_1,	cm_1_2,	cm_2_3,	cm_3_4,	cm_4_5,	over_cm_5), 
               names_to = "max_lenght", values_to = "count") %>%
  filter(count > 0) %>%
  uncount(weights = count, .remove = FALSE)%>%
  mutate(section = factor(section, levels = c("upper", "central", "lower")))  # Reorder section levels

#Create subsets for VFT and NOGR
data_lenght_dist_vft <- subset(data_lenght_dist, site == "VFT")
data_lenght_dist_nogr <- subset(data_lenght_dist, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight", "size_section_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~max_lenght + section, 
            data = data_lenght_dist_vft, 
            shade = TRUE,
            main = "Vicofertile",
            labeling_args = list(
              rot_labels = c(0, 30), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","size_section_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~max_lenght + section, 
            data = data_lenght_dist_nogr, 
            shade = TRUE,
            main = "Nogarole Rocca",
            labeling_args = list(
              rot_labels = c(0, 30), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#5. Mosaic Plot - Maximum_length ~ element ----

#Subsetting and organising data
data_lenght_elem <- phy_dataset %>%
  select(site, element_1, cm_0_1, cm_1_2, cm_2_3, cm_3_4, cm_4_5, over_cm_5) %>%
  pivot_longer(cols = c(cm_0_1, cm_1_2, cm_2_3, cm_3_4, cm_4_5, over_cm_5), 
               names_to = "max_lenght", values_to = "count") %>%
  filter(count > 0) %>%
  uncount(weights = count, .remove = FALSE) %>%
  filter(!element_1 %in% c("L", "C")) %>%  # Exclude L and C
  mutate(element_1 = factor(element_1, levels = c("CB", "LB", "SB", "UL", "LL", "T")))  # Reorder levels

#Create subsets for VFT and NOGR
data_lenght_elemt_vft <- subset(data_lenght_elem, site == "VFT")
data_lenght_elem_nogr <- subset(data_lenght_elem, site == "NOGR")

#Create the mosaic plot for VFT
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","size_element_VFT.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~max_lenght + element_1, 
            data = data_lenght_elemt_vft, 
            shade = TRUE,
            main = "Vicofertile",
            set_varnames = c(max_lenght = "max_lenght", element_1 = "element"),
            labeling_args = list(
              rot_labels = c(30, 30), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()

#Create the mosaic plot for NOGR
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","size_element_NOGR.jpeg"), #save as jpeg
     width = 700, height = 400)

vcd::mosaic(~max_lenght + element_1, 
            data = data_lenght_elem_nogr, 
            shade = TRUE,
            main = "Nogaole Rocca",
            set_varnames = c(max_lenght = "max_lenght", element_1 = "element"),
            labeling_args = list(
              rot_labels = c(30, 30), 
              offset_labels = c(0, 1.5),
              offset_varnames = c(0, 2),
              gp_labels = gpar(fontsize = 8)
            ))
dev.off()


#6. Fragmentation Index ----
##6.1 FI per Burial ----
#Calculate MNE per burial (number of unique identified bones)
mne_burial <- frag_data_burial %>%
  filter(!is.na(bone)) %>%
  distinct(burial, bone, side) %>%
  group_by(burial, bone) %>%
  summarise(
    has_L = any(side == "L"),
    has_R = any(side == "R"),
    has_I = any(side == "I"),
    .groups = "drop"
  ) %>%
  mutate(
    bone_count = case_when(
      has_L & has_R ~ 2,
      has_L | has_R ~ 1,
      has_I ~ 1,
      TRUE ~ 0
    )
  ) %>%
  group_by(burial) %>%
  summarise(MNE = sum(bone_count), .groups = "drop") %>%
  mutate(
    MNE = if_else(burial == "19", MNE * 2, MNE)  # Adjust MNE only for burial 19
  )

#Calculate NISP per burial
nisp_burial <- phy_dataset %>%
  group_by(burial) %>%
  summarise(NISP = sum(total_element))

#Join the two summaries and calculate the Fragmentation Index and rescale it (min-max normalization)
frag_index_burial <- left_join(nisp_burial, mne_burial, by = "burial") %>%
  mutate(Fragmentation_Index = NISP / MNE)

#Plotting the FI per burial
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","FI_burials.jpeg"), #save as jpeg
     width = 700, height = 400)

ggplot(
  frag_index_burial %>%
    mutate(
      site = case_when(
        burial %in% c("26n", "33(US359)", "33(US360)", "36", "56") ~ "Nogarole Rocca",
        TRUE ~ "Vicofertile"
      ),
      age_group = case_when(
        burial %in% c("33(US359)", "20", "16", "26v") ~ "Non-adult",
        burial == "19" ~ "Adult + Non-Adult",
        TRUE ~ "Adult"
      ),
      burial = factor(
        burial,
        levels = c(
          # Vicofertile in custom order
          "7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32",
          # Then Nogarole Rocca
          "26n", "33(US359)", "33(US360)", "36", "56"
        )
      )
    ),
  aes(x = burial, y = Fragmentation_Index, fill = age_group)
) +
  geom_col(position = "stack") +
  geom_hline(
    data = frag_index_burial %>%
      mutate(
        site = case_when(
          burial %in% c("26n", "33(US359)", "33(US360)", "36", "56") ~ "Nogarole Rocca",
          TRUE ~ "Vicofertile"
        )
      ) %>%
      group_by(site) %>%
      summarise(mean_FI = mean(Fragmentation_Index)),
    aes(yintercept = mean_FI),
    linetype = "dashed", color = "black", size = 1
  ) +
  facet_wrap(~site, scales = "free_x") +
  scale_fill_manual(32,
                    values = c(
                      "Adult" = "steelblue",
                      "Non-adult" = "salmon",
                      "Adult + Non-Adult" = "navyblue"
                    )
  ) +
  labs(
    title = "Fragmentation Index per Burial by Age Group and Site",
    x = "Burial",
    y = "Fragmentation Index",
    fill = "Age Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()


##6.2 FI per Burial and Element -----
#Calculate MNE
mne_burial_elem <- phy_dataset %>%
  filter(element_1 %in% c("CB","UL","LL","T")) %>%
  group_by(burial, element_1) %>%
  summarise(MNE = sum(total_element, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    MNE = if_else(burial == "19", MNE * 2, MNE)
  )

#Compute Fragmentation Index
frag_index_burial_elem <- mne_burial_elem %>%
  left_join(nisp_burial, by = "burial") %>%
  mutate(FI_element = NISP / MNE)

#Plotting the FI per burial and element by site

#NOGR 
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","FI_element_NOGR.jpeg"), #save as jpeg
     width = 800, height = 500)

ggplot(
  frag_index_burial_elem %>%
    filter(burial %in% c("26n", "33(US359)", "33(US360)", "36", "56")) %>%
    mutate(
      burial = factor(burial, levels = c("26n", "33(US359)", "33(US360)", "36", "56")),
      element_1 = factor(element_1, levels = c("CB", "UL", "LL", "T")),
      pattern = case_when(
        burial %in% c("33(US359)") ~ "stripe",
        TRUE ~ "none"
      )
    ),
  aes(x = burial, y = FI_element, fill = element_1, pattern = pattern)
) +
  ggpattern::geom_col_pattern(
    pattern_fill = "black",
    pattern_density = 0.1,
    pattern_spacing = 0.03,
    pattern_angle = 45,
    color = "black"
  ) +
  geom_hline(
    data = frag_index_burial_elem %>%
      filter(burial %in% c("26n", "33(US359)", "33(US360)", "36", "56")) %>%
      group_by(element_1) %>%
      summarise(mean_fi = mean(FI_element, na.rm = TRUE)) %>%
      mutate(element_1 = factor(element_1, levels = c("CB", "UL", "LL", "T"))),
    aes(yintercept = mean_fi),
    linetype = "dashed", color = "black", size = 1
  ) +
  facet_wrap(~element_1, scales = "free_x") +
  scale_fill_manual(values = c(
    "CB" = "#e14444",
    "UL" = "#a0e8b2",
    "LL" = "#3c7f49",
    "T"  = "#f0f050"
  )) +
  labs(
    title = "Fragmentation Index by Element – Nogarole Rocca Burials",
    x = "Burial",
    y = "Fragmentation Index"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
dev.off()

#VFT
jpeg(file = here("figures", "Chapter 5","Fragmentation_weight","FI_element_VFT.jpeg"), #save as jpeg
     width = 800, height = 500)

ggplot(
  frag_index_burial_elem %>%
    filter(burial %in% c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32")) %>%
    mutate(
      burial = factor(burial, levels = c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32")),
      element_1 = factor(element_1, levels = c("CB", "UL", "LL", "T")),
      pattern = case_when(
        burial %in% c("20", "16", "26v") ~ "stripe",
        burial == "19" ~ "circle",
        TRUE ~ "none"
      )
    ),
  aes(x = burial, y = FI_element, fill = element_1, pattern = pattern)
) +
  ggpattern::geom_col_pattern(
    pattern_fill = "black",
    pattern_density = 0.1,
    pattern_spacing = 0.03,
    pattern_angle = 45,
    color = "black"
  ) +
  geom_hline(
    data = frag_index_burial_elem %>%
      filter(burial %in% c("7", "8", "12", "15", "16", "17", "19", "20", "26v", "29", "30", "32")) %>%
      group_by(element_1) %>%
      summarise(mean_fi = mean(FI_element, na.rm = TRUE)) %>%
      mutate(element_1 = factor(element_1, levels = c("CB", "UL", "LL", "T"))),
    aes(yintercept = mean_fi),
    linetype = "dashed", color = "black", size = 1
  ) +
  facet_wrap(~element_1, scales = "free_x") +
  scale_fill_manual(values = c(
    "CB" = "#e14444",
    "UL" = "#a0e8b2",
    "LL" = "#3c7f49",
    "T"  = "#f0f050"
  )) +
  labs(
    title = "Fragmentation Index by Element – Vicofertile Burials",
    x = "Burial",
    y = "Fragmentation Index"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
dev.off()
