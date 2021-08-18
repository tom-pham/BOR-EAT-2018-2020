##%######################################################%##
#                                                          #
####      Survival Analysis for BOR-EAT 2018-2020       ####
#                                                          #
##%######################################################%##
# Author: Tom Pham
# Script and data info: This script performs survival analysis using CJS model
# and outputs detection probabilities using RMark
# Data source: Data extracted from NOAA ERDDAP data server (oceanview)
# Data details: StudyIDs examined: BC-Spring-2018, CNFH_FMR 2019/2020,
# ColemanLateFall 2018/2019, DeerCk_SH_Wild_2018, DeerCk_Wild_2018
# FR_Spring_2019, MillCk_Wild_2018, Mok_Fall_2018, Nimbus_Fall_2018,
# RBDD_2018, RBDD_WR_2018, SB_Spring 2018/2019, Winter_H 2018/2019/2020


source('helper.R')

#### ColemanLateFall 2018/2019  -----------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "ColemanLateFall_2018_19"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE)  

studyIDs <- c("ColemanLateFall_2018", "ColemanLateFall_2019")

# name for paths and files
name <- "ColemanLateFall_2018_19"

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/ColemanLateFall_2018_19/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# B/c two studyID find only common sites between them
common_sites <- all_detections %>% 
  bind_rows() %>% 
  select(StudyID, GEN, GenRKM) %>% 
  distinct() %>% 
  select(GEN, GenRKM) %>% 
  group_by(GEN) %>% 
  filter(n()>1) %>% 
  distinct()

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>%
  filter(
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
    # Choose only common sites OR release sites
    GEN %in% common_sites$GEN 
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

##### Detection Probability ----------------------------------------------------
outputs <- get_mark_model(all.inp, multiple = T, standardized = F)

p <- format_p(outputs, multiple = T)

# Identify receivers that had poor detection efficiencies
# bad_rec <- p %>% 
#   filter(estimate < 0.7) %>% 
#   filter(!(reach_end %in% c("GoldenGateE", "GoldenGateW", "BeniciaE", 
#                             "BeniciaW")))

# Compare estimates for each studyID
p_comparison <- compare_detection_probability_multi_studyid(p)

# Manually select sites to use and re-run
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BattleCk_CNFH_Rel", "BattleCk3", "Abv_Altube2", "Mill_Ck_Conf", 
               "Abv_WoodsonBr", "Blw_IrvineFinch", "BlwOrd", "ButteBr", 
               "Colusa AC2", "Colusa BC4", "Blw_Knights_GS3", 
               "Blw_Elkhorn_GS1", "TowerBridge","Hood",  "ChippsE", "ChippsW", "BeniciaE", 
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

outputs <- get_mark_model(all.inp, multiple = T, standardized = F)

p <- format_p(outputs, multiple = T)
write_csv(p_report, paste0("./Outputs/", name, "/", name, 
                           " Detection Probability.csv"))

p_report <- format_p_report(p, multiple = T)
write_csv(p_report, paste0("./Outputs/", name, "/", name, 
                           " Detection Probability Format.csv"))

plot_p(p, multiple = T)
ggsave(paste0("./Outputs/", name, "/", name, " Detection Probability.png"),
       width = 14, height = 10)

##### Reach survival per 10km --------------------------------------------------

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# write_csv(reach.meta.aggregate, "./Outputs/ColemanLateFall_2018_19/ColemanLateFall_2018_19_Sites.csv")

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

# Manually change regions
phi <- phi %>% 
  mutate(
    Region = case_when(
      reach_num == 1 ~ "Battle Ck",
      (reach_num >= 2 & reach_num <= 7) ~ "Upper Sac R",
      (reach_num >= 8 & reach_num <= 13) ~"Lower Sac R",
      reach_num == 14 ~ "Delta",
      reach_num >= 15 ~ "Bay",
      TRUE ~ Region
    )
  )

region_breaks <- c(1.5, 7.5, 13.5, 14.5)

plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Reach", multiple = T, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"), 
        width = 14, height = 10)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, reach_num, 'Count start' = count_start, "Count end" = count_end) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival rate per 10km (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Reach Survival per10km.csv"))



##### Cumulative Survival By Year-----------------------------------------------

cum_survival_2018 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "ColemanLateFall_2018")), T)
cum_survival_2019 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "ColemanLateFall_2019")), T)
cum_survival_all <- bind_rows(cum_survival_2018, cum_survival_2019)

cleanup(ask = FALSE)

cum_survival_all <- cum_survival_all %>% 
  add_column(
    GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], length(studyIDs)),
    GenRKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], length(studyIDs)),
    reach_num = rep(seq(0, (nrow(reach.meta.aggregate)-1), 1), length(studyIDs))
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>% 
      distinct()
  ) %>% 
  mutate_at(
    vars("cum.phi", "cum.phi.se", "LCI", 'UCI'), round, digits= 2
  ) %>% 
  mutate(
    Region = case_when(
      GEN == "Blw_Paynes_Ck" ~ "Upper Sac R",
      GEN == "ButteBr" ~ "Upper Sac R",
      GEN == "Chipps" ~ "Delta",
      GEN %in% c("Benicia", "GoldenGateE", "GoldenGateW") ~ "Bay",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("Battle Ck", "Upper Sac R", "Lower Sac R", "Delta", "Bay")),
    'Survival estimate (SE)' = paste0(cum.phi," (",  cum.phi.se, ")")
  ) %>% 
  filter(
    GEN != "GoldenGateW"
  )

region_breaks <- c(2.5, 7.5, 13.5, 14.5)

plot_cum_surv(cum_survival_all, add_breaks = T, multiple = T)


ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"), 
        width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(StudyID, 'Reach #' = reach_num, GEN, GenRKM, Region, 'Survival estimate (SE)', 
         LCI, UCI)

write_csv(cum_survival_all, paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv"))


##### Region Survival per 10km -------------------------------------------------

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BattleCk_CNFH_Rel", "ButteBr", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

# Manually change regions
phi$Region <- rep(c("Upper Sacramento", "Lower Sacramento", "Delta", "Bay", "Bay"), 2)
region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Region", multiple = T, text_size = 16)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km (+year).png"), 
        width = 10, height = 6)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival rate per 10km (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, 
          paste0("./Outputs/", name, "/", name, " Region Survival per10km (+year).csv"))


##### Region Survival ----------------------------------------------------------

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BattleCk_CNFH_Rel", "ButteBr", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

outputs <- get_mark_model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

# Manually change regions
phi$Region <- rep(c("Upper Sacramento", "Lower Sacramento", "Delta", "Bay", "Bay"), 2)
region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival", 
         xlabel = "Region", multiple = T, text_size = 16)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival (+year).png"), 
        width = 10, height = 6)


phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival rate (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, 
          paste0("./Outputs/", name, "/", name, " Region Survival (+year).csv"))



#### RBDD_WR_2018--------------------------------------------------------------

##### Create INP  --------------------------------------------------------------
name <- "RBDD_WR_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE)  

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("Benicia"),
                                         c("SacTrawl")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("BeniciaE", "BeniciaW"),
                                         c("SacTrawl1", "SacTrawl2")))

studyIDs <- c("RBDD_WR_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

names <- lapply(studyIDs, function(x) paste0("./outputs/RBDD_WR_2018/",
                                             x, ".csv"))

# Save detections to CSV 
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
    GEN != "I80-50_Br" # Only 1 fish, looks like false detect
  )

# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

# Select specific sites and re-run
reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("RBDD_Rel", "Blw_Salt", "Mill_Ck_Conf",
               "Abv_WoodsonBr", "GCID_abv", "Blw_IrvineFinch", "BlwOrd")
  ) 

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

outputs <- get_mark_model(all.inp, F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F) %>% 
  # Remove last reach
  slice(1:(n()-1))

# Plot detection probability
plot_p(p, multiple = F)

ggsave(paste0(path, name, " detection probability.png"),
       width = 12, height = 9, dpi = 500)

write_csv(p, paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv"))

##### Reach survival per 10km --------------------------------------------------

# Select specific sites
reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("RBDD_Rel", "Blw_Salt", "Mill_Ck_Conf",
               "Abv_WoodsonBr", "GCID_abv", "Blw_IrvineFinch", "BlwOrd")
  ) 


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1)) # Remove the last estimate which is not estimable

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 16)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

# Save the combined sites
write_csv(reach.meta.aggregate, 
          paste0(path, name, "_sites.csv"))

##### Cumulative Survival By Year-----------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), add_release = T)
cum_survival_2018

cum_survival_all <- format_cum_surv(cum_survival_2018) %>% 
  slice(1:(n()-1)) # Remove the last estimate which is not estimable

# Plot
plot_cum_surv(cum_survival_all, add_breaks = F, multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 7, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))


#### Winter_H 2018/2019  ------------------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "Winter_H_2018_19"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("Winter_H_2018", "Winter_H_2019")

# # Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

names <- lapply(studyIDs, function(x) paste0("./outputs/Winter_H_2018_19/",
                                             x, ".csv"))
# Save detections to CSV 
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# B/c two studyID find only common sites between them
common_sites <- all_detections %>% 
  bind_rows() %>% 
  select(StudyID, GEN, GenRKM) %>% 
  distinct() %>% 
  select(GEN, GenRKM) %>% 
  group_by(GEN) %>% 
  filter(n()>1) %>% 
  distinct()

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>%
  filter(
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
    # Choose only common sites OR release sites
    GEN %in% common_sites$GEN | GEN %in% c("Caldwell_Park_Rel", "Bonnyview_Rel")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, multiple = T, standardized = F)

p <- format_p(outputs, multiple = T)

# Identify receivers that had poor detection efficiencies
# bad_rec <- p %>%
#   filter(estimate < 0.7) %>%
#   filter(!(reach_end %in% c("GoldenGateE", "GoldenGateW", "BeniciaE",
#                             "BeniciaW")))

# Compare estimates for each studyID
p_comparison <- compare_detection_probability_multi_studyid(p)

# Manually select sites to use and re-run
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("Caldwell_Park_Rel", "Bonnyview_Rel", "Blw_ClearCr", "BlwCowCr", 
               "Blw_Paynes_Ck", "Blw_Salt", "GCID_abv", "Colusa AC3", "Colusa AC2",
               "AbvTisdale","BlwTisdale",  "BlwChinaBend", "Knights_RST", 
               "ChippsE", "ChippsW", "BeniciaE", "BeniciaW", "GoldenGateE", 
               "GoldenGateW"
    )
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

outputs <- get_mark_model(all.inp, multiple = T, standardized = F)
p <- format_p(outputs, multiple = T) %>% 
  filter(!(reach_num %in% c(1, 2)))
write_csv(p, paste0(path, name, 
                    " Detection Probability.csv"))

p_report <- format_p_report(p, multiple = T)
write_csv(p_report, paste0(path, name, 
                    " Detection Probability Format.csv"))

plot_p(p, multiple = T)
ggsave(paste0(path, name, " Detection Probability.png"),
       width = 12, height = 9, dpi = 500)


##### Reach survival per 10km --------------------------------------------------

# Must run survival analyses separately because of two different release locations
# with no common site at second release loc

###### Winter_H_2018 -----------------------------------------------------------

# Make a copy of reach.meta
reach.meta2 <- reach.meta
reach.meta <- reach.meta2 %>% filter(GEN != "Caldwell_Park_Rel")
all_detections <- vroom(names[1])

all_aggregated <- lapply(list(all_detections), aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp_2018 <- pmap(list(list(all_detections),all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp_2018, standardized = T, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2018"
phi_18 <- format_phi(outputs, multiple = F)

# Save the combined sites
write_csv(reach.meta.aggregate %>% 
            add_row(reach.meta2[1,]) %>% 
            arrange(desc(GenRKM)),
          paste0(path, name, "_sites.csv"))

###### Winter_H_2019 -----------------------------------------------------------
reach.meta <- reach.meta2
reach.meta <- reach.meta2 %>% filter(GEN != "Bonnyview_Rel")
all_detections <- vroom(names[2])

all_aggregated <- lapply(list(all_detections), aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp_2019 <- pmap(list(list(all_detections),all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp_2019, standardized = T, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2019"
phi_19 <- format_phi(outputs, multiple = F)

###### Combine 2018/2019  -----------------------------------------------------
phi <- rbind(phi_18, phi_19)

phi <- phi %>% 
  mutate(
    Reach = ifelse(reach_num == 1, "Bonnyview_Rel / \nCaldwell_Park_Rel to \n Blw_ClearCr", Reach),
    reach_start = ifelse(reach_num == 1, "Bonnyview_Rel/Caldwell_Park_Rel", reach_start),
    rkm_start = ifelse(reach_num == 1, "540.25/551.28", rkm_start),
    RKM = ifelse(reach_num == 1, "540.25/551.28 to 535.62", RKM),
  )

# Manually change regions
phi <- phi %>% 
  mutate(
    Region = case_when(
      reach_num <= 5 ~ "Upper Sac R",
      (reach_num >= 6 & reach_num <= 11) ~ "Lower Sac R",
      reach_num == 12 ~ "Delta",
      reach_num >= 13 ~ "Bay",
      TRUE ~ Region
    )
  )

write_csv(phi, paste0(path, name, 
                      " Reach Survival per10km.csv"))

region_breaks <- c(5.5, 11.5, 12.5)

plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Reach", multiple = T)

ggsave(paste0(path, name, " Reach Survival per10km.png"),
       width = 12, height = 9, dpi = 500)

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                           " Reach Survival per10km Format.csv"))

##### Cumulative Survival -----------------------------------------------------

cum_survival_2018 <- get_cum_survival(all.inp_2018, add_release = T)
cum_survival_2019 <- get_cum_survival(all.inp_2019, add_release = T)

cleanup(ask = F)

### Now combine cumulative survival
studyIDs <- c("Winter_H_2018", "Winter_H_2019")
cum_survival_all <- bind_rows(cum_survival_2018, cum_survival_2019)

cum_survival_all <- cum_survival_all %>% 
  add_column(
    GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], length(studyIDs)),
    GenRKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], length(studyIDs)),
    reach_num = rep(seq(1, (nrow(reach.meta.aggregate)), 1), length(studyIDs))
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>% 
      distinct()
  ) %>% 
  mutate_at(
    vars("cum.phi", "cum.phi.se", "LCI", 'UCI'), round, digits= 2
  ) %>% 
  mutate(
    Region = case_when(
      reach_num <= 6 ~ "Upper Sac R",
      (reach_num >= 6 & reach_num <= 12) ~ "Lower Sac R",
      reach_num == 13 ~ "Delta",
      reach_num >= 14 ~ "Bay",
      TRUE ~ Region
    ),
    GEN = ifelse(reach_num == 1, "Bonnyview_Rel /\nCaldwell_Park_Rel", GEN),
    GenRKM = ifelse(reach_num == 1, "540.25/551.28", GenRKM)
  ) %>% 
  filter(
    GEN != 'GoldenGateW'
  )

region_breaks <- c(6.5, 12.5, 13.5)

plot_cum_surv(cum_survival_all, add_breaks = T, multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 8, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(StudyID, 'Reach #' = reach_num, GEN, GenRKM, Region, 'Survival estimate (SE)', 
         LCI, UCI)

write_csv(cum_survival_all, "./Outputs/Winter_H_2018_2019/Winter_H_2018_2019 Cumulative Survival.csv")


##### Regional survival per 10km ----------------------------------------------

###### Winter_H_2018 -----------------------------------------------------------

reach.meta <- reach.meta2 %>% 
  filter(GEN %in% c("Bonnyview_Rel", "Colusa AC3", "Knights_RST", "ChippsE",
                    "ChippsW", "GoldenGateE", "GoldenGateW"))

all_aggregated_18 <- lapply(list(all_detections[[1]]), aggregate_GEN)
all_EH_18 <- lapply(all_aggregated_18, make_EH)
all.inp_18 <- pmap(list(list(all_detections[[1]]),all_EH_18), create_inp) %>% 
  bind_rows()

# Reach survival for 2018
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp_18, standardized = T, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2018"
phi_18 <- format_phi(outputs, multiple = F)


###### Winter_H_2019 -----------------------------------------------------------

reach.meta <- reach.meta2 %>% 
  filter(GEN %in% c("Caldwell_Park_Rel", "Colusa AC3", "Knights_RST", "ChippsE",
                    "ChippsW", "GoldenGateE", "GoldenGateW"))
all_aggregated_19 <- lapply(list(all_detections[[2]]), aggregate_GEN)
all_EH_19 <- lapply(all_aggregated_19, make_EH)
all.inp_19 <- pmap(list(list(all_detections[[2]]),all_EH_19), create_inp) %>% 
  bind_rows()

# Reach survival for 2019
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp_19, standardized = T, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2019"
phi_19 <- format_phi(outputs, multiple = F)

phi_combined <- bind_rows(phi_18, phi_19) %>% 
  mutate(
    reach_start = gsub("Bonnyview_Rel", "Release", reach_start),
    reach_start = gsub("Caldwell_Park_Rel", "Release", reach_start),
    Reach = gsub("Bonnyview_Rel", "Release", Reach),
    Reach = gsub("Caldwell_Park_Rel", "Release", Reach),
    RKM = gsub("540.25", "551.28/540.25", RKM)
  )

phi_combined <- bind_rows(phi_18, phi_19) %>% 
  mutate(
    reach_start = gsub("Bonnyview_Rel", "Release", reach_start),
    reach_start = gsub("Caldwell_Park_Rel", "Release", reach_start),
    Reach = gsub("Bonnyview_Rel", "Release", Reach),
    Reach = gsub("Caldwell_Park_Rel", "Release", Reach),
    RKM = gsub("540.25", "551.28/540.25", RKM)
  )


# Manually change regions
phi_combined <- phi_combined %>% 
  mutate(
    Region = case_when(
      reach_num == 3 ~ "Delta",
      Region %in% c("West Delta", "SF Bay") ~ "Bay",
      Region == "Upper Sac R" ~ "Upper Sacramento",
      Region == "Lower Sac R"~ "Lower Sacramento",
      TRUE ~ Region
    )
  )

region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi_combined, type = "Region", add_breaks = T, 
         ylabel = "Survival per 10km", xlabel = "Region", multiple = T,
         text_size = 16)

ggsave("./Outputs/Winter_H_2018_19/Winter_H_2018_19 Region Survival per10km.png", 
        width = 10, height = 6)


phi_table <- phi_combined %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename(
    'Survival rate per 10km (SE)' = Estimate,
    'Reach #' = reach_num
  ) %>% 
  select(-SE)

write_csv(phi_table, "./Outputs/Winter_H_2018_19/Winter_H_2018_19 Region Survival per10km.csv")


##### Regional Survival -------------------------------------------------------

# 2018
outputs <- get_mark_model(all.inp_18, standardized = F, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2018"
phi_18 <- format_phi(outputs, multiple = F)

# 2019
outputs <- get_mark_model(all.inp_19, standardized = F, multiple = F)
cleanup(ask = F)

studyIDs <- "Winter_H_2019"
phi_19 <- format_phi(outputs, multiple = F)

phi_combined <- bind_rows(phi_18, phi_19) %>% 
  mutate(
    reach_start = gsub("Bonnyview_Rel", "Release", reach_start),
    reach_start = gsub("Caldwell_Park_Rel", "Release", reach_start),
    Reach = gsub("Bonnyview_Rel", "Release", Reach),
    Reach = gsub("Caldwell_Park_Rel", "Release", Reach),
    RKM = gsub("540.25", "551.28/540.25", RKM)
  )

phi_combined <- bind_rows(phi_18, phi_19) %>% 
  mutate(
    reach_start = gsub("Bonnyview_Rel", "Release", reach_start),
    reach_start = gsub("Caldwell_Park_Rel", "Release", reach_start),
    Reach = gsub("Bonnyview_Rel", "Release", Reach),
    Reach = gsub("Caldwell_Park_Rel", "Release", Reach),
    RKM = gsub("540.25", "551.28/540.25", RKM)
  )


# Manually change regions
phi_combined <- phi_combined %>% 
  mutate(
    Region = case_when(
      reach_num == 3 ~ "Delta",
      Region %in% c("West Delta", "SF Bay") ~ "Bay",
      Region == "Upper Sac R" ~ "Upper Sacramento",
      Region == "Lower Sac R"~ "Lower Sacramento",
      TRUE ~ Region
    )
  )

region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi_combined, type = "Region", add_breaks = T, ylabel = "Survival", 
         xlabel = "Region", multiple = T, text_size = 16)

ggsave("./Outputs/Winter_H_2018_19/Winter_H_2018_19 Region Survival.png", width = 10, height = 6, dpi = 500)


phi_table <- phi_combined %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename(
    'Survival rate per 10km (SE)' = Estimate,
    'Reach #' = reach_num
  ) %>% 
  select(-SE)

write_csv(phi_table, "./Outputs/Winter_H_2018_2019/Winter_H_2018_2019 Region Survival.csv")


#### Mok_Fall_2018 ------------------------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "Mok_Fall_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("Mok_Fall_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/Mok_Fall_2018/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))



reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("Sherman_Island_Rel", "ChippsE", "ChippsW", "BeniciaE", "BeniciaW",
               "GoldenGateE", "GoldenGateW")
  )

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("BeniciaE", "BeniciaW")))

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

p_format <- format_tbl_report(p, multiple = F)

# Plot detection probability
plot_p(p, multiple = F, text_size = 20)

ggsave(paste0(path, name, " Detection Probability.png"),
       width = 10, height = 8, dpi = 500)


# dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p_format, paste0(paste0("./Outputs/", name, "/", name,
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

region_breaks <- c(1.5, 2.5)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))


##### Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format_cum_surv(cum_survival_all)

# Plot
plot_cum_surv(cum_survival_all, add_breaks = F, multiple = F)



ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 8, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))


#### SB_Spring_2018_19 -------------------------------------------------------

##### Create INP  -----------------------------------
name <- "SB_Spring_2018_2019"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("SB_Spring_2018", "SB_Spring_2019")

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

names <- lapply(studyIDs, function(x) paste0(path,
                                             x, ".csv"))

# Save detections to CSV 
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    # GEN %in% common_sites$GEN,
    # !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
    #   GEN %in% c("ChippsE", "ChippsW")
    GEN %in% c("SutterBypass_Weir2_RST_Rel", "Butte1", "Butte2", "Butte3",
               "Butte6", "TowerBridge", "Hood", "ChippsE", "ChippsW", 
               "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

leaflet(data = reach.meta) %>% addTiles() %>%
  addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
             label = ~as.character(GEN),
             labelOptions = labelOptions(noHide = T, textOnly = TRUE))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, F, multiple = T)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = T)

p_comparison <- compare_detection_probability_multi_studyid(p)

# Plot detection probability
plot_p(p, multiple = TRUE)

# Identify and remove sites with low p <0.7
reach.meta <- reach.meta %>% 
  filter(
    !(GEN %in% c("Abv_FremontWeir", "Blw_FremontWeir"))
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Then rerun 
outputs <- get_mark_model(all.inp, F, multiple = T)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = T)

p <- p %>% 
  mutate(
    reach_start = ifelse(reach_num == 1, 'Release',
                         reach_start),
    Reach = ifelse(reach_num == 1, 
                   str_replace(Reach, 'SutterBypass_Weir2_RST_Rel', 'Release'),
                   Reach)
  )

# Plot detection probability
plot_p(p, multiple = T)

ggsave(paste0(path, name, " Detection Probability.png"),
       width = 12, height = 9, dpi = 500)

p <- format_tbl_report(p, multiple = T)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("SutterBypass_Weir2_RST_Rel", "Butte1", "Butte2", "Butte3",
               "Butte6", "TowerBridge", "Hood", "ChippsE", "ChippsW", 
               "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
  ) %>% 
  mutate(GEN = ifelse(GEN == "SutterBypass Weir2 RST", "SutterBypass \nWeir2 RST", GEN))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

region_breaks <- c(4.5, 6.5, 7.5)

# Rename first site to 'Release' bc it is too long
phi <- phi %>% 
  mutate(
    reach_start = ifelse(reach_num == 1, 'Release',
                         reach_start),
    Reach = ifelse(reach_num == 1, 
                   str_replace(Reach, 'SutterBypass_Weir2_RST_Rel', 'Release'),
                   Reach)
  )

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = T, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 12, height = 9, dpi = 500)

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))



##### Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "SB_Spring_2018")), T)
cum_survival_2018

cum_survival_2019 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "SB_Spring_2019")), T)
cum_survival_2019

cum_survival_all <- bind_rows(cum_survival_2018, cum_survival_2019) %>% 
  mutate_at(vars("cum.phi", "cum.phi.se", "LCI", "UCI"), round, digits = 2)

cum_survival_all <- format_cum_surv(cum_survival_all)

# Rename first site to 'Release' bc it is too long
cum_survival_all <- cum_survival_all %>% 
  mutate(
    GEN = ifelse(reach_num == 0, 'Release',
                         GEN)
  )

# Plot
region_breaks <- c(5.5, 7.5, 8.5)
plot_cum_surv(cum_survival_all, add_breaks = T, multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 7, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    StudyID, 'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )


cum_survival_all <- pivot_cum_surv(cum_survival_all)
write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))
##### Region survival  ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("SutterBypass_Weir2_RST_Rel", "Butte6", "Hood", "ChippsE", 
               "ChippsW", "GoldenGateE", "GoldenGateW")
  ) 

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()


###### Survival per 10km -------------------------------------------------------
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get_mark_model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_per10km$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot_phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival per 10km", multiple = T)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 8, height = 6, dpi = 500)

phi_per10km_table <- make_phi_table(phi_per10km, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                           " Region Survival per 10km.csv")))


###### Survival per region -----------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_region$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

plot_phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival",
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 8, height = 6, dpi = 500)


phi_region_table <- make_phi_table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))

#### Nimbus_Fall_2018 --------------------------------------------------------

##### Create INP -------------------------------------------------------------
name <- "Nimbus_Fall_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("Nimbus_Fall_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/ColemanLateFall_2018_19/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Remove Mare Island release fish
all_detections <- lapply(all_detections, function(x) {
  x %>% 
    bind_rows() %>% 
    left_join(
      TaggedFish %>% 
        select(fish_id, release_location),
      by = c("FishID" = "fish_id")
    ) %>% 
    filter(release_location != "Mare_Island_Rel") 
})

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("AR_Sunrise_Ramp_Rel", "TowerBridge", "SacTrawl1",
               "SacTrawl2", "Hood", "ChippsE", "ChippsW", "BeniciaE", 
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)



# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))
# 

### Manually create EH with slight variation to ignore a release group ###

min_detects <- all_aggregated %>% 
  bind_rows() %>% 
  filter(GEN %in% reach.meta.aggregate$GEN) %>% 
  group_by(FishID, GEN, GenRKM) %>% 
  summarise(
    min_time = min(time)
  ) %>% 
  arrange(
    FishID, min_time
  )

# Get list of all tagged fish for the studyID
fish <- TaggedFish %>% 
  filter(
    study_id == studyIDs,
    release_location != 'Mare_Island_Rel'
  ) %>% 
  arrange(fish_id) %>% 
  pull(fish_id)

# Create matrix of all combinations of fish and GEN
EH <- expand.grid(
  fish,
  reach.meta.aggregate$GEN, stringsAsFactors = FALSE 
)

names(EH) <- c('FishID', 'GEN')  

# Add col detect to min_detects, these fish get a 1
min_detects$detect <- 1

# Join in detections to the matrix, fish detected a GEN will be given a 1
# otherwise it will be given a 0
EH <- EH %>% 
  left_join(
    min_detects %>% 
      select(
        FishID, GEN, detect
      ), by = c("FishID", "GEN")
  ) %>% 
  # Replace NA with 0 https://stackoverflow.com/questions/28992362/dplyr-join-define-na-values
  mutate_if(
    is.numeric, coalesce, 0
  )

# Reshape the df wide, so that columns are GEN locations, rows are fish, 
# values are 1 or 0 for presence/absence
EH <- reshape(EH, idvar = 'FishID', timevar = 'GEN', direction = 'wide')
colnames(EH) <- gsub('detect.', '', colnames(EH))
# Manually make the release column a 1 because all fish were released there
# sometimes detections df does not reflect that accurately
EH[2] <- 1
all_EH <- list(EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

region_breaks <- c(1.5, 3.5, 4.5)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))


##### Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018) %>% 
  mutate_at(
    vars("cum.phi", "cum.phi.se", "LCI", 'UCI'), round, digits= 2
  ) 

cum_survival_all <- format_cum_surv(cum_survival_all)

# Plot
region_breaks <- c(1.5, 4.5, 5.5)
plot_cum_surv(cum_survival_all, add_breaks = T, multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 8, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

##### Region survival  ----------------------------------------------------------------

region_sites <-  c("AR_Sunrise_Ramp_Rel", "TowerBridge", "Hood", "Chipps", 
                   "GoldenGateE", "GoldenGateW")

reach.meta.aggregate <- reach.meta.aggregate %>% 
  filter(GEN %in% region_sites)
 
all_EH <- all_EH[[1]] %>% 
  select(FishID, AR_Sunrise_Ramp_Rel, TowerBridge, Hood, Chipps, GoldenGateE,
         GoldenGateW)

all_EH <- list(all_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()


###### Survival per 10km -------------------------------------------------------
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = F) %>% 
  filter(reach_num != 5) %>% 
  mutate(Region = c("American River", "Lower Sacramento", "Delta", "Bay"))

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot_phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = F)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

# phi_per10km_table <- make_phi_table(phi_per10km, standardized = T)
# write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
#                                            " Region Survival per 10km.csv")))

###### Survival per region -----------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = F) %>% 
  filter(reach_num != 5) %>% 
  mutate(Region = c("American River", "Lower Sacramento", "Delta", "Bay"))

plot_phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 10, height = 6, dpi = 500)


phi_region_table <- make_phi_table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))


##### Mare Island Release Reach Survival per 10km -----------------------------

# Reread in the original detections
all_detections <- lapply(names, vroom)

# Use Mare Island release fish
all_detections <- lapply(all_detections, function(x) {
  x %>% 
    bind_rows() %>% 
    left_join(
      TaggedFish %>% 
        select(fish_id, release_location),
      by = c("FishID" = "fish_id")
    ) %>% 
    filter(release_location == "Mare_Island_Rel") 
})

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("Mare_Island_Rel", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

min_detects <- all_aggregated %>% 
  bind_rows() %>% 
  filter(GEN %in% reach.meta.aggregate$GEN) %>% 
  group_by(FishID, GEN, GenRKM) %>% 
  summarise(
    min_time = min(time)
  ) %>% 
  arrange(
    FishID, min_time
  )

# Get list of all tagged fish for the studyID
fish <- TaggedFish %>% 
  filter(
    study_id == studyIDs,
    release_location == 'Mare_Island_Rel'
  ) %>% 
  arrange(fish_id) %>% 
  pull(fish_id)

# Create matrix of all combinations of fish and GEN
EH <- expand.grid(
  fish,
  reach.meta.aggregate$GEN, stringsAsFactors = FALSE 
)

names(EH) <- c('FishID', 'GEN')  

# Add col detect to min_detects, these fish get a 1
min_detects$detect <- 1

# Join in detections to the matrix, fish detected a GEN will be given a 1
# otherwise it will be given a 0
EH <- EH %>% 
  left_join(
    min_detects %>% 
      select(
        FishID, GEN, detect
      ), by = c("FishID", "GEN")
  ) %>% 
  # Replace NA with 0 https://stackoverflow.com/questions/28992362/dplyr-join-define-na-values
  mutate_if(
    is.numeric, coalesce, 0
  )

# Reshape the df wide, so that columns are GEN locations, rows are fish, 
# values are 1 or 0 for presence/absence
EH <- reshape(EH, idvar = 'FishID', timevar = 'GEN', direction = 'wide')
colnames(EH) <- gsub('detect.', '', colnames(EH))
# Manually make the release column a 1 because all fish were released there
# sometimes detections df does not reflect that accurately
EH[2] <- 1
all_EH <- list(EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# phi_table <- make_phi_table(phi, standardized = F)

phi_table <- format_tbl_report(phi, multiple = F)

write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Mare Island Reach Survival.csv")))

#### BC-Spring-2018 ----------------------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "BC-Spring-2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("BC-Spring-2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/BC-Spring-2018/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Remove duplicated release site
reach.meta <- reach.meta %>% 
  filter(
    GEN != "UpperButte_RST",
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)
p <- format_p(outputs, multiple = F)

plot_p(p, multiple = F)

ggsave(paste0(path, name, " Detection Probability.png"),
       width = 10, height = 8, dpi = 500)

write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)

phi <- format_phi(outputs, multiple = F)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 8, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))

##### Cumulative Survival ----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- cum_survival_2018

cum_survival_all <- format_cum_surv(cum_survival_all)

# Plot
plot_cum_surv(cum_survival_all, add_breaks = F, multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 8, height = 8, dpi = 500)

write_csv(cum_survival_all, "./Outputs/BC-Spring-2018/BC-Spring-2018 Cumulative Survival.csv")

#### DeerCk_Wild_CHK_2018 -----------------------------------------------------


##### Create INP ------------------------------------------------------------------
name <- "DeerCk_Wild_CHK_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("DeerCk_Wild_CHK_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/DeerCk_Wild_CHK_2018/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>% 
  filter(
    GEN != "DeerCk_RST_Rel",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("DeerCk_RST","DeerCk3","Abv_WoodsonBr",
               "GCID_abv", "Blw_IrvineFinch", "BlwOrd", "ButteBr", "Colusa AC2")
  )


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Then rerun 
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

p_format <- format_tbl_report(p, multiple = F)

write_csv(p_format, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------

KM <- reach.meta.aggregate$GenRKM
# reach_length <- abs(diff(KM))/10

# Have to manually enter reach_length for some reason using abs(diff(KM))
# gives issue when running mark model
reach_length <- c(0.886, 0.801, 1.904, 1.145, 3.294, 1.762, 2.549)
outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1))

region_breaks <- c(2.5)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

phi_report <- format_tbl_report(phi, multiple = F)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))

##### Cumulative Survival ----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- cum_survival_2018

cum_survival_all <- format_cum_surv(cum_survival_all)

# Plot
region_breaks <- 2.5
plot_cum_surv(cum_survival_all%>% 
                slice(1:(n()-1)), add_breaks = T, multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 7, dpi = 500)

write_csv(cum_survival_all, paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv"))


##### Region Survival ---------------------------------------------------------

###### Survival per 10km-------------------------------

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("DeerCk_RST", "Abv_WoodsonBr", "ButteBr", "Colusa AC3")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1))

region_breaks <- c(1.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"), 
        width = 8, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival rate per 10km (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Region Survival per10km.csv"))

###### Region Survival --------------------------------------------------------

outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1)) # remove last row

# Manually change regions
region_breaks <- c(1.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival", 
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"), 
        width = 8, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Region Survival.csv"))


#### MillCk_Wild_CHK_2018  -----------------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "MillCk_Wild_CHK_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("MillCk_Wild_CHK_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/BC-Spring-2018/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Remove duplicated release site
reach.meta <- reach.meta %>% 
  filter(
    GEN != "MillCk_RST",
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)
p <- format_p(outputs, multiple = F)

plot_p(p, multiple = F)

ggsave(paste0(path, name, " Detection Probability.png"),
       width = 10, height = 8, dpi = 500)

write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)

phi <- format_phi(outputs, multiple = F)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F, text_size = 20)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 8, height = 8, dpi = 500)

phi_table <- make_phi_table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km Format.csv"))

##### Cumulative Survival ----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- cum_survival_2018

cum_survival_all <- format_cum_surv(cum_survival_all)

# Plot
plot_cum_surv(cum_survival_all, add_breaks = F, multiple = F)

cum_survival_all <- format_cum_surv(cum_survival_all)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 8, height = 8, dpi = 500)

write_csv(cum_survival_all, paste0(path, name, 
                             " Cumulative Survival.csv"))

#### RBDD_2018 ----------------------------------------------------------------

##### Create INP --------------------------------------------------------------
name <- "RBDD_2018"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("RBDD_2018")

# Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/RBDD_2018/",
                                             x, ".csv"))
# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>%
  filter(
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

# Select sites
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("RBDD_Rel", "Blw_Salt",  "Mill_Ck_Conf", "Abv_WoodsonBr", "GCID_abv", 
               "Blw_IrvineFinch", "BlwOrd", "ButteBr", "Colusa AC2", "Colusa BC2",
               "AbvTisdale", "Knights_RST", "Blw_FRConf", "I80-50_Br", "Freeport",
               "Hood")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

write_csv(reach.meta.aggregate, 
          paste0(path, name, "_sites.csv"))

outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)
p <- format_p(outputs, multiple = F)

# Plot detection probability
plot_p(p, multiple = F)

# dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
# write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
#                            " detection probability.csv")))

##### Reach survival per 10km ----------------------------------------------------------------

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) 

region_breaks <- c(7.5)

# Plot survival probability
plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = F)

# Create folder for studyid, save plot
ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 12, height = 8, dpi = 500)

phi_report <- format_tbl_report(phi, multiple = F) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
  )

write_csv(phi_report, paste0(path, name, 
                           " Reach Survival per10km.csv"))
##### Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format_cum_surv(cum_survival_all) %>% 
  slice(1:(n()-1))

# Plot
region_breaks <- 8.5
plot_cum_surv(cum_survival_all, add_breaks = T, multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"), 
       width = 10, height = 7, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

##### Region survival  ----------------------------------------------------------------

###### Survival per 10km -------------------------------------------------------

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("RBDD_Rel", "ButteBr", "Freeport", "Hood")
  ) 

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1))

region_breaks <- c(1.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per 10km.png"), 
        width = 8, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival per 10km (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Region Survival per 10km.csv"))


###### Region Survival  -------------------------------------------------------
outputs <- get_mark_model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F) %>% 
  slice(1:(n()-1))

region_breaks <- c(1.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival", 
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"), 
        width = 8, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Region Survival.csv"))

#### CNFH_FMR 2019_2020 ------------------------------------------------------

##### Create INP -----------------------------------
name <- "CNFH_FMR_2019_20"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("CNFH_FMR_2019", "CNFH_FMR_2020")

# # Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# # Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/CNFH_FMR_2019_20/",
                                             x, ".csv"))

# pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# B/c two studyID find only common sites between them
common_sites <- all_detections %>% 
  bind_rows() %>% 
  select(StudyID, GEN, GenRKM) %>% 
  distinct() %>% 
  select(GEN, GenRKM) %>% 
  group_by(GEN) %>% 
  filter(n()>1) %>% 
  distinct()

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>%
  filter(
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
    # Choose only common sites OR release sites
    GEN %in% common_sites$GEN | GEN %in% c("Caldwell_Park_Rel", "Bonnyview_Rel")
  )


all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()


##### Detection probability ---------------------------------------------------
outputs <- get_mark_model(all.inp, multiple = T, standardized = F)

p <- format_p(outputs, multiple = T)


# Identify receivers that had poor detection efficiencies
bad_rec <- p %>%
  filter(estimate < 0.7) %>%
  filter(!(reach_end %in% c("GoldenGateE", "GoldenGateW", "BeniciaE",
                            "BeniciaW")))

# Compare estimates for each studyID
df <- p %>% 
  select(-c("se", "reach_num", "Reach", "RKM", "Region", "lcl", "ucl")) %>% 
  pivot_wider(names_from = StudyID, values_from = c("estimate", "count_start", "count_end"),
              names_glue = "{StudyID} {.value}")

prefixes <- unique(p$StudyID)

names_to_order <- map(prefixes, ~ names(df)[grep(paste0(.x, " "), names(df))]) %>% unlist
names_id <- setdiff(names(df), names_to_order)

df <- df %>%
  select(names_id, names_to_order)

# Manually select sites to use
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("RBDD_Rel", "Blw_Salt", "Mill_Ck_Conf", "Abv_WoodsonBr",
               "GCID_abv", "GCID_blw", "Blw_IrvineFinch", "BlwOrd", "ButteBr", 
               "Colusa AC3", "AbvTisdale", "BlwChinaBend", 
               "Knights_RST", "Blw_FRConf", "TowerBridge",
               "TowerBr", "BeniciaE", 
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )


all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

write_csv(reach.meta.aggregate, 
          paste0(path, name, "_sites.csv"))


outputs <- get_mark_model(all.inp, multiple = T, standardized = F)

p <- format_p(outputs, multiple = T)

plot_p(p, multiple = T)

ggsave(paste0(path, name, " Detection Probability.png"),
       width = 12, height = 8, dpi = 500)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

p_report <- format_tbl_report(p, multiple = TRUE)
write_csv(p_report, paste0(path, name, 
                             " Detection Probability Format.csv"))



##### Reach survival per 10km ----------------------------------------------------------------

write_csv(reach.meta.aggregate, "./Outputs/CNFH_FMR_2019_20/CNFH_FMR_2019_20_Sites.csv")

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

# Manually change regions
phi <- phi %>% 
  mutate(
    Region = case_when(
      rkm_end >= 344.100 ~ "Upper Sac R",
      rkm_end < 344.100 & rkm_end >= 138.220 ~ "Lower Sac R",
      rkm_end < 138.220 & rkm_end >= 52.140 ~ "Delta",
      rkm_end < 52.140 ~ "Bay",
      TRUE ~ Region
    )
  )

region_breaks <- c(8.5, 14.5, 15.5)

plot_phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Reach", multiple = T)

ggsave("./Outputs/CNFH_FMR_2019_20/CNFH_FMR_2019_20 Reach Survival per10km (+year).png", 
        width = 12, height = 8, dpi = 500)

phi_report <- format_tbl_report(phi, multiple = TRUE)
write_csv(phi_report, paste0(path, name, 
                           " Reach Survival per10km (+year).csv"))

##### Cumulative Survival By Year----------------------------------------------------------------
# reach.meta <- get_receiver_GEN(all_detections)

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()


cum_survival_2019 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "CNFH_FMR_2019")), T)
cum_survival_2020 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "CNFH_FMR_2020")), T)

cum_survival_all <- bind_rows(cum_survival_2019, cum_survival_2020)

cleanup(ask = FALSE)

cum_survival_all <- cum_survival_all %>% 
  add_column(
    GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], length(studyIDs)),
    GenRKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], length(studyIDs)),
    reach_num = rep(seq(1, (nrow(reach.meta.aggregate)), 1), length(studyIDs))
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>% 
      distinct()
  ) %>% 
  mutate_at(
    vars("cum.phi", "cum.phi.se", "LCI", 'UCI'), round, digits= 2
  ) %>% 
  mutate(
    Region = case_when(
      GenRKM >= 344.100 ~ "Upper Sac R",
      GenRKM < 344.100 & GenRKM >= 138.220 ~ "Lower Sac R",
      GenRKM < 138.220 & GenRKM >= 52.140 ~ "Delta",
      GenRKM < 52.140 ~ "Bay",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("Upper Sac R", "Lower Sac R", "Delta", "Bay")),
    'Survival estimate (SE)' = paste0(cum.phi," (",  cum.phi.se, ")")
  ) %>% 
  filter(
    GEN != "GoldenGateW"
  )

region_breaks <- c(9.5, 15.5, 16.5)

plot_cum_surv(cum_survival_all, add_breaks = T, multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"), 
        width = 10, height = 7, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(StudyID, 'Reach #' = reach_num, GEN, GenRKM, Region, 'Survival estimate (SE)', 
         LCI, UCI)

write_csv(cum_survival_all, paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv"))

##### Region Survival per 10km-------------------------------

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("RBDD_Rel", "ButteBr", "TowerBridge", "BeniciaE", "BeniciaW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)


# Manually change regions
phi$Region <- rep(c("Upper Sacramento", "Lower Sacramento", "Delta", "Bay", "Bay"), 2)
region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival per 10km", 
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"), 
        width = 10, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival rate per 10km (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                                   " Region Survival per10km.csv"))

##### Region Survival --------------------------------------------------------

outputs <- get_mark_model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)


# Manually change regions
phi$Region <- rep(c("Upper Sacramento", "Lower Sacramento", "Delta", "Bay", "Bay"), 2)
region_breaks <- c(1.5, 2.5, 3.5)

plot_phi(phi, type = "Region", add_breaks = T, ylabel = "Survival", 
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"), 
        width = 10, height = 6, dpi = 500)

phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, 
         UCI = ucl, 'Reach #' = reach_num) %>% 
  mutate(
    Reach = str_remove_all(Reach, "\n"),
    Estimate = round(Estimate, 2),
    SE = round(SE, 2),
    LCI = round(LCI, 2),
    UCI = round(UCI, 2),
    Estimate = paste0(Estimate, " (", as.character(SE), ")")
  ) %>% 
  rename('Survival (SE)' = Estimate) %>% 
  select(-SE)

write_csv(phi_table, paste0("./Outputs/", name, "/", name, 
                            " Region Survival.csv"))

#### FR_Spring_2019 ------------------------------------------------------

##### Create INP -----------------------------------
name <- "FR_Spring_2019"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 
studyIDs <- c("FR_Spring_2019")

# # Retreive ERDDAP data if first time
# all_detections <- lapply(studyIDs, get_detections)

# # Save detections to CSV 
names <- lapply(studyIDs, function(x) paste0("./outputs/CNFH_FMR_2019_20/",
                                             x, ".csv"))

pmap(list(all_detections, names), write_csv)

# Load CSV's if not first time
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>%
  filter(
    GEN != "Colusa BC2", # Above release 
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
  )


###### Boyds -------------------------------------------------------------------
# Split up INP by release groups Boyds and Gridley

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BoydsPump", "BC_Beach",
               "Blw_FRConf", "I80-50_Br", "Hood", "ChippsE", "ChippsW", 
               "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

reach.meta.boyds <- reach.meta.aggregate

boyds <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Boyds_Rel")


###### Gridley -----------------------------------------------------------------
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("FR_Gridley_Rel", "SunsetPumps", "BoydsPump", "BC_Beach",
               "Blw_FRConf", "I80-50_Br", "Hood", "ChippsE", "ChippsW", 
               "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

reach.meta.gridley <- reach.meta.aggregate

gridley <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Gridley_Rel")

##### Detection probability ---------------------------------------------------
# outputs <- get_mark_model(all.inp, multiple = F, standardized = F)
# 
# p <- format_p(outputs, multiple = F)
# 
# 
# # Identify receivers that had poor detection efficiencies
# bad_rec <- p %>%
#   filter(estimate < 0.7) %>%
#   filter(!(reach_end %in% c("GoldenGateE", "GoldenGateW", "BeniciaE",
#                             "BeniciaW")))
# 
# # Compare estimates for each studyID
# df <- p %>% 
#   select(-c("se", "reach_num", "Reach", "RKM", "Region", "lcl", "ucl")) %>% 
#   pivot_wider(names_from = StudyID, values_from = c("estimate", "count_start", "count_end"),
#               names_glue = "{StudyID} {.value}")
# 
# prefixes <- unique(p$StudyID)
# 
# names_to_order <- map(prefixes, ~ names(df)[grep(paste0(.x, " "), names(df))]) %>% unlist
# names_id <- setdiff(names(df), names_to_order)
# 
# df <- df %>%
#   select(names_id, names_to_order)
# 
# # Manually select sites to use
# reach.meta <- reach.meta %>% 
#   filter(
#     GEN %in% c("RBDD_Rel", "Blw_Salt", "Mill_Ck_Conf", "Abv_WoodsonBr",
#                "GCID_abv", "GCID_blw", "Blw_IrvineFinch", "BlwOrd", "ButteBr", 
#                "Colusa AC3", "AbvTisdale", "BlwChinaBend", 
#                "Knights_RST", "Blw_FRConf", "TowerBridge",
#                "TowerBr", "BeniciaE", 
#                "BeniciaW", "GoldenGateE", "GoldenGateW")
#   )
# 
# 
# all_aggregated <- lapply(all_detections, aggregate_GEN)
# 
# all_EH <- lapply(all_aggregated, make_EH)
# 
# all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
#   bind_rows()
# 
# write_csv(reach.meta.aggregate, 
#           paste0(path, name, "_sites.csv"))
# 
# 
# outputs <- get_mark_model(all.inp, multiple = T, standardized = F)
# 
# p <- format_p(outputs, multiple = T)
# 
# plot_p(p, multiple = T)
# 
# ggsave(paste0(path, name, " Detection Probability.png"),
#        width = 12, height = 8, dpi = 500)
# 
# dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
# write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
#                            " detection probability.csv")))
# 
# p_report <- format_tbl_report(p, multiple = TRUE)
# write_csv(p_report, paste0(path, name, 
#                            " Detection Probability Format.csv"))
# 
# 

##### Reach survival per 10km ----------------------------------------------------------------

# write_csv(reach.meta.aggregate, "./Outputs/CNFH_FMR_2019_20/CNFH_FMR_2019_20_Sites.csv")

###### Boyds -------------------------------------------------------------------

reach.meta.aggregate <- reach.meta.boyds %>% 
  filter(
    GEN %in% c("BoydsPump", "BC_Beach",
               "Blw_FRConf", "I80-50_Br", "Hood", "Chipps", 
               "Benicia", "GoldenGateE", "GoldenGateW")
  )

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(boyds, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Manually change regions
phi_boyds <- phi %>% 
  mutate(
    Region = case_when(
      rkm_end < 138.220 & rkm_end >= 52.140 ~ "Delta",
      rkm_end < 52.140 ~ "Bay",
      TRUE ~ Region
    )
  )


###### Gridley -----------------------------------------------------------------

reach.meta.aggregate <- reach.meta.gridley

# Get river KM
KM <- reach.meta.gridley$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(gridley, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Manually change regions
phi_gridley <- phi %>% 
  mutate(
    Region = case_when(
      rkm_end < 138.220 & rkm_end >= 52.140 ~ "Delta",
      rkm_end < 52.140 ~ "Bay",
      TRUE ~ Region
    )
  )


###### Combine Survival Estimates ----------------------------------------------
phi_boyds$Release <- "Boyds"
phi_gridley$Release <- "Gridley"

phi_combined <- rbind(phi_boyds, phi_gridley)

# Create the levels order 
lvls <- phi_combined %>% 
  select(Reach, rkm_start) %>% 
  distinct() %>% 
  arrange(desc(rkm_start)) %>% 
  pull(Reach)

p <- phi_combined
p <- p %>% 
  # Filter out the GGE to GGW estimate, not really useful
  filter(reach_end != "GoldenGateW") %>% 
  mutate(
    Reach = factor(Reach, levels = lvls)
  )

# p["Reach"] <- factor(rep(lvls, length(studyIDs)), levels = lvls)

angle <- 90
hjust <- 1
vjust <- 0.5

add_breaks = TRUE
ylabel = "Survival per 10km"
xlabel = "Reach"
text_size = 20
padding = 20
region_breaks <- c(4.5, 6.5, 7.5)

sites <- rbind(reach.meta.boyds, reach.meta.gridley)
sites <- sites %>%
  distinct() %>% 
  arrange(desc(GenRKM))

write_csv(sites, paste0("./Outputs/", name, "/", name, 
                        " Sites.csv"))


# leaflet(data = x) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))
# 

ggplot(data = p, mapping = aes(x = Reach, y = estimate, group = Release)) +
  geom_point(aes(color = Release), position = position_dodge(.5), size = 3) +
  geom_errorbar(mapping = aes(x = Reach, ymin = lcl, ymax = ucl, 
                              color = Release),  width = .1,
                position = position_dodge(.5)) +
  # Conditionally add breaks 
  {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
  ylab(ylabel) +
  xlab(xlabel) +
  ## WILL NEED TO FIX THIS, ONLY WORKS FOR 2 STUDYIDS
  scale_color_manual(values=c("#007EFF", "#FF8100")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
    plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
    legend.position = "top",
    text = element_text(size=text_size)
  )

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km (+year).png"), 
        width = 10, height = 7, dpi = 500)


# Format estimates for report table
prefixes <- c("Gridley", "Boyds")

phi_report <- phi_combined %>% 
  select(-reach_num) %>%
  mutate(
    estimate_se = paste0(round(estimate, digits = 2), " (", 
                         round(se, digits = 2), ")"),
    lcl = round(lcl, digits = 2),
    ucl = round(ucl, digits = 2),
    Reach = str_replace(Reach, "([\n])", "")
  ) %>% 
  select(-c("estimate", "se")) %>% 
  pivot_wider(names_from = Release, values_from = c("estimate_se",  
                                                    "lcl", "ucl"),
              names_glue = "{Release} {.value}") %>% 
  arrange(desc(rkm_start)) %>% 
  select(-c("reach_start", "reach_end", "rkm_start", 'rkm_end', 'count_start',
            'count_end'))

names_to_order <- map(prefixes, ~ names(phi_report)[grep(paste0(.x, " "), names(phi_report))]) %>% unlist
names_id <- setdiff(names(phi_report), names_to_order)

phi_report <- phi_report %>%
  select(names_id, names_to_order)

write_csv(phi_report, paste0(path, name, 
                             " Reach Survival per10km (+year).csv"))

##### Cumulative Survival By Year-----------------------------------------------
cum_survival_boyds <- get_cum_survival(boyds, T)

cum_survival_gridley <- get_cum_survival(gridley, T)

cum_survival_all <- bind_rows(cum_survival_boyds, cum_survival_gridley) %>% 
  select(-StudyID)


cleanup(ask = FALSE)

cum_survival_all <- cum_survival_all %>% 
  mutate(Release = c(rep("Boyds", 9), rep("Gridley", 11))) %>% 
  add_column(
    GEN = c(reach.meta.boyds$GEN, reach.meta.gridley$GEN),
    GenRKM = c(reach.meta.boyds$GenRKM, reach.meta.gridley$GenRKM),
    Region = c(reach.meta.boyds$Region, reach.meta.gridley$Region)
  ) %>% 
  mutate_at(
    vars("cum.phi", "cum.phi.se", "LCI", 'UCI'), round, digits= 2
  ) %>% 
  mutate(
    Region = case_when(
      GenRKM >= 344.100 ~ "Upper Sac R",
      GenRKM < 344.100 & GenRKM >= 138.220 ~ "Lower Sac R",
      GenRKM < 138.220 & GenRKM >= 52.140 ~ "Delta",
      GenRKM < 52.140 ~ "Bay",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("Upper Sac R", "Lower Sac R", "Delta", "Bay")),
    'Survival estimate (SE)' = paste0(cum.phi," (",  cum.phi.se, ")")
  ) %>% 
  filter(
    GEN != "GoldenGateW"
  )

all_sites <- bind_rows(reach.meta.boyds, reach.meta.gridley) %>% 
  distinct() %>% 
  arrange(desc(GenRKM))

region_breaks <- c(4.5, 6.5, 7.5)

add_breaks = TRUE
padding = 5.5
text_size = 20

lvls <- all_sites %>% 
  pull(GEN)

cum_survival_all %>% 
  mutate(
    GEN = factor(GEN, levels = lvls)
  ) %>%
  ggplot(mapping = aes(x = GEN, y = cum.phi, group = Release)) +
  geom_point(size = 2, aes(color = Release)) +
  geom_errorbar(mapping = aes(x= GEN, ymin = LCI, ymax = UCI, 
                              color = Release),  width = .1) +
  geom_line(size = 0.7, aes(color = Release)) +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  ylab("Cumulative survival") +
  xlab("Receiver Location") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_color_manual(values=c("#007EFF", "#FF8100")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(5.5, 5.5, 5.5, 20, "pt"),
    legend.position = "top",
    text = element_text(size=20)
  ) 

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"), 
        width = 10, height = 7, dpi = 500)

# cum_survival_all <- cum_survival_all %>% 
#   select(StudyID, 'Reach #' = reach_num, GEN, GenRKM, Region, 'Survival estimate (SE)', 
#          LCI, UCI)

# Format estimates for report table
prefixes <- c("Gridley", "Boyds")

phi_report <- cum_survival_all %>% 
  select(-c("cum.phi", "cum.phi.se")) %>% 
  pivot_wider(names_from = Release, values_from = c("Survival estimate (SE)",  
                                                    "LCI", "UCI"),
              names_glue = "{Release} {.value}") %>% 
  arrange(desc(GenRKM)) 

names_to_order <- map(prefixes, ~ names(phi_report)[grep(paste0(.x, " "), names(phi_report))]) %>% unlist
names_id <- setdiff(names(phi_report), names_to_order)

phi_report <- phi_report %>%
  select(names_id, names_to_order)

write_csv(phi_report, paste0(path, name, 
                             " Cumulative Survival.csv"))

##### Region Survival per 10km -------------------------------------------------

###### Gridley -----------------------------------------------------------------

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("FR_Gridley_Rel", "Blw_FRConf", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# reach.meta.gridley <- reach.meta.aggregate

gridley <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Gridley_Rel")

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(gridley, standardized = T, multiple = F)
cleanup(ask = F)

phi_gridley <- format_phi(outputs, multiple = F) %>% 
  filter(reach_end != "GoldenGateW")

# Manually change regions
phi_gridley$Region <- c("Feather River", "Lower Sacramento", "Delta", "Bay")
# region_breaks <- c(1.5, 2.5, 3.5)


###### Boyds -------------------------------------------------------------------
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BoydsPump", "Blw_FRConf", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

boyds <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Boyds_Rel")

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

outputs <- get_mark_model(boyds, standardized = T, multiple = F)
cleanup(ask = F)

phi_boyds <- format_phi(outputs, multiple = F) %>% 
  filter(reach_end != "GoldenGateW")

# Manually change regions
phi_boyds$Region <- c("Feather River", "Lower Sacramento", "Delta", "Bay")
region_breaks <- c(1.5, 2.5, 3.5)


###### Combined Survival Estimates ---------------------------------------------

combined_region_survival <- bind_rows(phi_gridley, phi_boyds) %>% 
  mutate(
    Release = c(rep("Gridley", 4), rep("Boyds", 4))
  )

# Create the levels order 
lvls <- combined_region_survival %>% 
  select(Region, rkm_start) %>% 
  arrange(desc(rkm_start)) %>% 
  pull(Region) %>% 
  unique()

combined_region_survival <- combined_region_survival %>% 
  mutate(
    Region = factor(Region, levels = lvls)
  )

angle <- 0
hjust <- 0.5
vjust <- 0

add_breaks = TRUE
ylabel = "Survival per 10km"
xlabel = "Region"
text_size = 20
padding = 5.5

ggplot(data = combined_region_survival, mapping = aes(x = Region, y = estimate, 
                                                      group = Release)) +
  geom_point(aes(color = Release), position = position_dodge(.5), size = 3) +
  geom_errorbar(mapping = aes(x = Region, ymin = lcl, ymax = ucl, 
                              color = Release),  width = .1,
                position = position_dodge(.5)) +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  ylab(ylabel) +
  xlab(xlabel) +
  scale_color_manual(values=c("#007EFF", "#FF8100")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
    plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
    legend.position = "top",
    text = element_text(size=text_size)
  )

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"), 
        width = 10, height = 6, dpi = 500)


# Format estimates for report table
prefixes <- c("Gridley", "Boyds")

phi_report <- combined_region_survival %>% 
  mutate(
    estimate_se = paste0(round(estimate, digits = 2), " (", 
                         round(se, digits = 2), ")"),
    lcl = round(lcl, digits = 2),
    ucl = round(ucl, digits = 2),
    Reach = str_replace(Reach, "([\n])", "")
  ) %>% 
  select(-c("estimate", "se")) %>% 
  pivot_wider(names_from = Release, values_from = c("estimate_se",  
                                                    "lcl", "ucl"),
              names_glue = "{Release} {.value}") %>% 
  arrange(desc(rkm_start)) %>% 
  select(-c("reach_start", "reach_end", "rkm_start", 'rkm_end', 'count_start',
            'count_end'))

names_to_order <- map(prefixes, ~ names(phi_report)[grep(paste0(.x, " "), names(phi_report))]) %>% unlist
names_id <- setdiff(names(phi_report), names_to_order)

phi_report <- phi_report %>%
  select(names_id, names_to_order)

write_csv(phi_report, paste0(path, name, 
                             " Region Survival per10km.csv"))

##### Region Survival --------------------------------------------------------

###### Gridley -----------------------------------------------------------------
reach.meta <- get_receiver_GEN(all_detections)
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("FR_Gridley_Rel", "Blw_FRConf", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# reach.meta.gridley <- reach.meta.aggregate

gridley <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Gridley_Rel")

outputs <- get_mark_model(gridley, standardized = F, multiple = F)
cleanup(ask = F)

phi_gridley <- format_phi(outputs, multiple = F) %>% 
  filter(reach_end != "GoldenGateW")

# Manually change regions
phi_gridley$Region <- c("Feather River", "Lower Sacramento", "Delta", "Bay")
# region_breaks <- c(1.5, 2.5, 3.5)


###### Boyds -------------------------------------------------------------------
reach.meta <- get_receiver_GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BoydsPump", "Blw_FRConf", "Hood", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

boyds <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Boyds_Rel")

outputs <- get_mark_model(boyds, standardized = F, multiple = F)
cleanup(ask = F)

phi_boyds <- format_phi(outputs, multiple = F) %>% 
  filter(reach_end != "GoldenGateW")

# Manually change regions
phi_boyds$Region <- c("Feather River", "Lower Sacramento", "Delta", "Bay")
region_breaks <- c(1.5, 2.5, 3.5)


###### Combined Survival Estimates ---------------------------------------------

combined_region_survival <- bind_rows(phi_gridley, phi_boyds) %>% 
  mutate(
    Release = c(rep("Gridley", 4), rep("Boyds", 4))
  )

# Create the levels order 
lvls <- combined_region_survival %>% 
  select(Region, rkm_start) %>% 
  arrange(desc(rkm_start)) %>% 
  pull(Region) %>% 
  unique()

combined_region_survival <- combined_region_survival %>% 
  mutate(
    Region = factor(Region, levels = lvls)
  )

angle <- 0
hjust <- 0.5
vjust <- 0

add_breaks = TRUE
ylabel = "Survival per 10km"
xlabel = "Region"
text_size = 20
padding = 5.5

ggplot(data = combined_region_survival, mapping = aes(x = Region, y = estimate, 
                                                      group = Release)) +
  geom_point(aes(color = Release), position = position_dodge(.5), size = 3) +
  geom_errorbar(mapping = aes(x = Region, ymin = lcl, ymax = ucl, 
                              color = Release),  width = .1,
                position = position_dodge(.5)) +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  ylab(ylabel) +
  xlab(xlabel) +
  scale_color_manual(values=c("#007EFF", "#FF8100")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
    plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
    legend.position = "top",
    text = element_text(size=text_size)
  )

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"), 
        width = 10, height = 6, dpi = 500)


# Format estimates for report table
prefixes <- c("Gridley", "Boyds")

phi_report <- combined_region_survival %>% 
  mutate(
    estimate_se = paste0(round(estimate, digits = 2), " (", 
                         round(se, digits = 2), ")"),
    lcl = round(lcl, digits = 2),
    ucl = round(ucl, digits = 2),
    Reach = str_replace(Reach, "([\n])", "")
  ) %>% 
  select(-c("estimate", "se")) %>% 
  pivot_wider(names_from = Release, values_from = c("estimate_se",  
                                                    "lcl", "ucl"),
              names_glue = "{Release} {.value}") %>% 
  arrange(desc(rkm_start)) %>% 
  select(-c("reach_start", "reach_end", "rkm_start", 'rkm_end', 'count_start',
            'count_end'))

names_to_order <- map(prefixes, ~ names(phi_report)[grep(paste0(.x, " "), names(phi_report))]) %>% unlist
names_id <- setdiff(names(phi_report), names_to_order)

phi_report <- phi_report %>%
  select(names_id, names_to_order)

write_csv(phi_report, paste0(path, name, 
                             " Region Survival.csv"))
