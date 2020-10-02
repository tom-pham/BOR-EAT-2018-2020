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
# ColemanLateFall 2018/2019/2020, DeerCk_SH_Wild_2018, DeerCk_Wild_2018
# FR_Spring_2019, MillCk_Wild_2018, Mok_Fall_2018, Nimbus_Fall_2018,
# RBDD_2018, RBDD_WR_2018, SB_Spring 2018/2019, Winter_H 2018/2019/2020


library(RMark)
library(tidyverse)
library(rerddap)
library(lubridate)
library(clusterPower)
# library(RColorBrewer)
# library(cowplot)
library(RODBC)
# library(flextable)
library(leaflet)


### Load TaggedFish and ReceiverDeployments tables through ERDDAP ---------------------------------------------------------------------

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo, url = my_url)  

JSATSinfo <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, url = my_url)

# Establish ERDDAP url and database name
my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_detects', url = my_url)

# Retrieve list of all studyIDs on FED_JSATS
studyid_list <- tabledap(JSATSinfo,
                         fields = c('study_id'),
                         url = my_url,
                         distinct = TRUE
) %>% 
  filter(study_id != "2017_BeaconTag") %>% 
  pull(study_id)


# Load functions ----------------------------------------------------------

get_detections <- function(studyID) {
  # Retrieve detection data from ERDDAP
  #
  # Arguments:
  #  studyID: StudyID name to retrieve data for
  #     
  # Return:
  #  df of detection data formatted correctly, add in RKM, Region, Lat, Lon, 
  # Release RKM, Release Lat, Release Lon, format types, rename cols
  
  df <- tabledap(JSATSinfo,
                 fields = c('study_id', 'fish_id', 'receiver_general_location',
                            'time'),
                 paste0('study_id=', '"',studyID, '"'),
                 url = my_url,
                 distinct = T
  ) %>% 
    left_join(
      ReceiverDeployments %>% 
        select(receiver_general_location, receiver_general_river_km, receiver_region,
               receiver_general_latitude, receiver_general_longitude)
    ) %>% 
    left_join(
      TaggedFish %>% 
        select(fish_id, release_river_km, release_latitude, release_longitude, 
               release_location) %>% distinct()
    )
  
  # Rename columns and change column types as ERDDAP returns data all in 
  # character format
  df <- df %>% 
    rename(
      StudyID = study_id,
      FishID = fish_id,
      GEN = receiver_general_location,
      GenRKM = receiver_general_river_km,
      Region = receiver_region,
      GenLat = receiver_general_latitude,
      GenLon =receiver_general_longitude,
      RelRKM = release_river_km,
      Rel_loc = release_location
    ) %>% 
    mutate(
      GenLat = ifelse(is.na(GenLat), release_latitude, GenLat),
      GenLon = ifelse(is.na(GenLon), release_longitude, GenLon),
      GenLat = as.numeric(GenLat),
      GenLon = as.numeric(GenLon),
      GenRKM = as.numeric(GenRKM),
      RelRKM = as.numeric(RelRKM),
      time = ymd_hms(time),
      GenRKM = ifelse(is.na(GenRKM), RelRKM, GenRKM)
    )
  
  # ERDDAP by default returns a table.dap object which does not play nice with
  # maggittr (pipes) so convert to tibble
  as_tibble(df)
  
}

# get_GEN_locs <- function(detections){
#   detections %>% 
#     select(StudyID, GEN, GenRKM) %>% 
#     distinct() %>% 
#     arrange(desc(GenRKM))
# }
# 
# all_GEN_locs <- lapply(all_detections, get_GEN_locs)

aggregate_GEN <- function(detections) {
  # Replace GEN in detections df according to replacelist
  #
  # Arguments:
  #  detections: a detections df
  #     
  # Return:
  #  a detections df that has replaced list of GEN with the aggregated GEN,
  #  mean RKM, mean Lat, mean Lon. Creates reach.meta.aggregate which is the list
  #  of receiver sites with new aggregated GEN's, along with RKM, Lat, Lon
  
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <<- reach.meta
  
  # Walk through each key/pair value
  for (i in 1:length(replace_dict$replace_with)) {
    # Unlist for easier to use format
    replace_list <- unlist(replace_dict[[2]][i])
    replace_with <- unlist(replace_dict[[1]][i])
    
    # Gather receiver data for the replace_list and replace, get the mean genrkm. 
    # This will be used to replace in the detections
    replace <- reach.meta %>% 
      select(GEN, GenRKM, GenLat, GenLon, Region) %>% 
      filter(GEN %in% c(replace_list, replace_with)) %>%
      distinct() %>% 
      select(-GEN) %>% 
      group_by(Region) %>% 
      summarise_all(mean)
    
    # Replace replace_list GENs name with replace_with GEN, and replace all of 
    # their genrkm with the averaged val
    detections <- detections %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
      )
    
    # This new df shows receiver metadata and reflects the aggregation done
    reach.meta.aggregate <<- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN == "End", "End", ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region))
      ) %>% 
      distinct()
  }
  detections
}


make_EH <- function(detections) {
  # Make an encounter history df
  #
  # Arguments:
  #  detections: a detections df
  #     
  # Return:
  #  Encounter history df. A matrix of every fish tagged for a given studyID
  #  at every given receiver site (that is in reach.meta.aggregate) and whether
  #  it was present 1 or absent 0 in the detection df
  
  # Get earliest detection for each fish at each GEN
  min_detects <- detections %>% 
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
    filter(study_id == detections$StudyID[1]) %>% 
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
  EH
}

create_inp <- function(detections, EH) { 
  # Create an inp df
  #
  # Arguments:
  #  detections: a detections df
  #  EH: an encounter history df
  #     
  # Return:
  #  inp df i.e. Fish01 | 11101, a record of a fish and it's presence/absence
  #  at each given receiver location. 
  
  EH.inp <- EH %>% 
    # Collapse the encounter columns into a single column of 1's and 0's
    unite("ch", 2:(length(EH)), sep ="") %>% 
    # Use the detections df to get the StudyID assignment
    mutate(StudyID = unique(detections$StudyID))
  EH.inp
}


get.mark.model <- function(all.inp, standardized, multiple) {
  # Run a CJS Mark model
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  standardized: TRUE or FALSE, if you want outputs to be standardized to 
  #     per10km or not
  #  mutliple: TRUE or FALSE, if you have multiple studyIDs or not   
  #
  # Return:
  #  the outputs of running a CJS Mark model, df with phi and p estimates, LCI
  #  UCI, SE
  
  # For single studyID
  if (multiple == F) {
    # If standardized, set time.intervals to reach_length to get per 10km
    if (standardized) {
      all.process <- process.data(all.inp, model="CJS", begin.time=1, 
                                  time.intervals = reach_length)
    } else {
      all.process <- process.data(all.inp, model="CJS", begin.time=1)
    }
  # For multiple studyID
  } else {
    # If multiple studyIDs, set groups to "StudyID
    if (standardized) {
      all.process <- process.data(all.inp, model="CJS", begin.time=1, 
                                  time.intervals = reach_length, groups = "StudyID")
    } else {
      all.process <- process.data(all.inp, model="CJS", begin.time=1,
                                  groups = "StudyID")
    }
  }


  all.ddl <- make.design.data(all.process)
  rm(list=ls(pattern="p.t.x.y"))
  rm(list=ls(pattern="Phi.t.x.y"))
  
  if (multiple) {
    p.t.x.y <- list(formula= ~time*StudyID)
    Phi.t.x.y <- list(formula= ~time*StudyID)
  }else {
    p.t.x.y <- list(formula= ~time)
    Phi.t.x.y <- list(formula= ~time)
  }

  cml = create.model.list("CJS")
  model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) 
  outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real

}

get_cum_survival <- function(all.inp, add_release) {
  # Run a CJS Mark model for cumulative survival
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  add_release: TRUE or FALSE, if you wish to add an extra dummy row at the top
  #  to show 100% survival at the release location
  #
  # Return:
  #  Cumulative survival outputs of CJS Mark model
  
  all.process <- process.data(all.inp, model = "CJS", begin.time = 1)
  all.ddl <- make.design.data(all.process)
  
  rm(list=ls(pattern="Phi.t"))
  rm(list=ls(pattern="p.t"))
  
  p.t <- list(formula= ~time) 
  Phi.t <- list(formula= ~time)
  
  cml = create.model.list("CJS")
  
  model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl, realvcv = TRUE)
  
  reaches <- nchar(all.inp$ch[1]) - 1
  
  phi.t <- model.outputs$Phi.t.p.t$results$real$estimate[1:reaches] 
  phi.t.vcv <- model.outputs$Phi.t.p.t$results$real.vcv
  
  cum.phi <- cumprod(phi.t)
  
  # calculate standard errors for the cumulative product. 
  cum.phi.se <- deltamethod.special("cumprod", phi.t[1:reaches], 
                                    phi.t.vcv[1:(reaches),1:(reaches)])
  
  
  ### Output estimate, SE, LCI, UCI to a dataframe
  cumulative <- data.frame(cum.phi = cum.phi, cum.phi.se = cum.phi.se, 
                           LCI = expit(logit(cum.phi)-1.96*sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                           UCI = expit(logit(cum.phi)+1.96*sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))))
  
  # Round to 3 digits
  cumulative <- round(cumulative,3)
  
  # If add_release TRUE, add in the dummy row to the top which just represents
  # survival of 100% at release
  if (add_release == T) {
    cumulative <- cumulative %>% 
      add_row(
        .before = 1,
        cum.phi = 1,
        cum.phi.se = NA,
        LCI = NA, 
        UCI = NA
      ) %>% 
      mutate(
        StudyID = all.inp$StudyID[1]
      )
  }else {
    cumulative <- cumulative %>% 
      mutate(
        StudyID = all.inp$StudyID[1]
      )
  }
  
}

get.receiver.GEN <- function(all_detections) {
  # Get a list of all receiver sites and metadata for a given detections df
  #
  # Arguments:
  #  all_detections: detections df 
  #
  # Return:
  #  df of receiver sites along with RKM, Lat, Lon, Region
  
  reach.meta <- all_detections %>% 
    bind_rows() %>% 
    distinct(GEN, GenRKM, GenLat, GenLon, Region) %>% 
    # Necessary because detections files shows differing RKM, Lat, Lon for some 
    # GEN sometimes
    group_by(GEN) %>% 
    summarise(
      GenRKM = mean(GenRKM),
      GenLat = mean(GenLat),
      GenLon = mean(GenLon),
      Region = first(Region)
    ) %>% 
    arrange(desc(GenRKM))
}


format.p <- function(output, multiple){
  # Format p outputs for plotting and table outputs
  #
  # Arguments:
  #  output: output df from Mark model
  #  mutliple: TRUE/FALSE if there were multiple StudyIDs in the outputs
  #
  # Return:
  #  properly formatted df for p outputs, now ready to plot
  
  # Format for single
  if (multiple == F) {
    outputs %>%
      slice(
        # Grab bottom half of outputs which are the p values
        ((nrow(outputs)/2) +1):nrow(outputs)
      ) %>%
      select(
        -c("fixed", "note")
      ) %>%
      # Add in some metadata so values make more sense
      add_column(
        reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
        reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
        rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
        rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end),
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      mutate(reach_num = 1:n())
  } else {
    # If there are multiple StudyIDs, formatting is same idea just slightly 
    # different
    outputs %>%
      slice(
        ((nrow(outputs)/2) +1):nrow(outputs)
      ) %>%
      select(
        -c("fixed", "note")
      ) %>%
      rownames_to_column(var = "StudyID") %>%
      add_column(
        reach_start = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)], 2),
        reach_end = rep(reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))], 2),
        rkm_start = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)], 2),
        rkm_end = rep(reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))], 2)
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end),
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      rowwise() %>%
      mutate(
        StudyID = strsplit(strsplit(StudyID, "p g")[[1]][2], " ")[[1]][1]
      )
  }
}

format_phi <- function(outputs, multiple) {
  # Format phi outputs for plotting and table outputs
  #
  # Arguments:
  #  output: output df from Mark model
  #  mutliple: TRUE/FALSE if there were multiple StudyIDs in the outputs
  #
  # Return:
  #  properly formatted df for phi outputs, now ready to plot
  
  # Format for single
  if (multiple == F) {
    outputs %>%
      # Grab first half of Mark outputs which represent Phi values
      slice(1:(nrow(outputs) / 2)) %>%
      select(
        -c("fixed", "note")
      ) %>% 
      add_column(
        StudyID = studyIDs,
        reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
        reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
        rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
        rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end),
        reach_num = 1:n()
      ) 
  } else {
    outputs %>%
      slice(1:(nrow(outputs) / 2)) %>%
      select(
        -c("fixed", "note")
      ) %>% 
      rownames_to_column(var = "StudyID") %>%
      add_column(
        reach_start = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)], 2),
        reach_end = rep(reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))], 2),
        rkm_start = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)], 2),
        rkm_end = rep(reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))], 2),
        reach_num = rep(1:(nrow(reach.meta.aggregate)-1),length(all_EH))
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end)
      ) %>% 
      rowwise() %>%
      mutate(
        StudyID = strsplit(strsplit(StudyID, "Phi g")[[1]][2], " ")[[1]][1]
      )
  }

}

get.unique.detects <- function(all_aggregated){
  # Get raw number of unique fish detected at each GEN 
  #
  # Arguments:
  #  all_aggregated: df of detections that have been replaced with aggregating
  #  receiver locations
  #
  # Return:
  #  df of each GEN in a detections df and the raw number of unique fish detected
  
  all_aggregated %>% 
    bind_rows() %>% 
    select(StudyID, FishID, GEN, GenRKM) %>% 
    distinct() %>% 
    group_by(StudyID, GEN, GenRKM) %>% 
    summarise(
      count = n()
    ) %>% 
    arrange(StudyID, desc(GenRKM)) %>% 
    ungroup()
}

# DEPRECATED - does not function 100%
# get.region.breaks <- function(phi) {
#   phi %>% 
#     group_by(Region) %>% 
#     summarise(breaks = max(reach_num)) %>% 
#     mutate(breaks = breaks + .5) %>% 
#     arrange(breaks) %>% 
#     slice(1:(n()-1)) %>% 
#     pull(breaks)
# }

plot.phi <- function(phi, type, add_breaks, ylabel, xlabel, multiple, 
                     padding = 5.5) {
  # Plot phi outputs from Mark model
  #
  # Arguments:
  #  phi: phi outputs from Mark model, must be formatted with (format_phi) first
  #  type: "Reach" or "Region", dictates how the plot will be created
  #  add_breaks: TRUE/FALSE whether to add vertical line breaks to represent
  #     regions
  #  ylabel: label for y axis
  #  xlabel: label for x axis
  #  multiple: TRUE/FALSE, whether there are multiple studyIDs or not
  #  padding: leftside plot margin, default set to 5.5 good for most, but can
  #     be adjusted of the xaxis label too long and gets cut off
  #
  # Return:
  #  plot of phi with estimate and error bars representing LCI, UCI
  
  # Create the levels order 
  lvls <- phi %>% select(type) %>% distinct() %>% pull()
  
  # # Unfortunately, this method of variable assignment doesn't work when I have
  # # Multiple studyIDs, something to do with recycling rules in mutate()
  # p <- phi %>%
  #   mutate(
  #     # !! allows me to select the column by the arg I pass (Reach or Region),
  #     # := goes in hand with assignment via !!
  #     # https://github.com/tidyverse/ggplot2/releases/tag/v3.0.0
  #     !!type := factor(lvls, levels = lvls)
  #   ) %>%
  #   filter(reach_end != "GoldenGateW")
  
  # Does same thing as above in base R
  p <- phi
  p[type] <- factor(rep(lvls, length(studyIDs)), levels = lvls)
  p <- p %>% 
    # Filter out the GGE to GGW estimate, not really useful
    filter(reach_end != "GoldenGateW")

  if (type == "Reach") {
    # If plotting for Reach survival set angle of xaxis to be 45 degrees because
    # they are too long
    angle <- 45
    hjust <- 1
  }else {
    angle <- 0
    hjust <- 0.5
  }

  if (multiple == T) {
    ggplot(data = p, mapping = aes(x = get(type), y = estimate, group = StudyID)) +
      geom_point(aes(color = StudyID), position = position_dodge(.5)) +
      geom_errorbar(mapping = aes(x = get(type), ymin = lcl, ymax = ucl, 
                                  color = StudyID),  width = .1,
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
        axis.text.x = element_text(angle = angle, hjust = hjust),
        plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
        legend.position = c(.12 ,.9)
      )
      
  }  else {
    ggplot(data = p, mapping = aes(x = get(type), y = estimate)) +
      geom_point() +
      geom_errorbar(mapping = aes(x= get(type), ymin = lcl, ymax = ucl), 
                    width = .1) +
      # Conditionally add breaks 
      {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
      ylab(ylabel) +
      xlab(xlabel) +
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=.75),
        axis.text.x = element_text(angle = angle, hjust = hjust)
      )
  }
}

make.phi.table <- function(phi, standardized = T) {
  # Format phi outputs further to be ready to save as a csv
  #
  # Arguments:
  #  phi: phi outputs from Mark model, must be formatted with (format_phi) first
  #
  # Return:
  #  phi df formatted the way I want to be saved as csv
  ifelse(standardized, label <-  'Survival rate per 10km (SE)',
         label <- 'Survival rate (SE)')
  
  phi %>% 
    select(StudyID, reach_num, Reach, RKM, Region, Estimate = estimate, SE = se, 
           LCI = lcl, UCI = ucl, N= count) %>% 
    mutate(
      Reach = str_remove_all(Reach, "\n"),
      Estimate = round(Estimate, 2),
      SE = round(SE, 2),
      LCI = round(LCI, 2),
      UCI = round(UCI, 2),
      Estimate = paste0(Estimate, " (", as.character(SE), ")")
    ) %>% 
    rename(!!label := Estimate,
           'Reach #' = reach_num) %>% 
    select(-SE)
}

plot.p <- function(p) {
  # Plot p outputs from Mark model
  #
  # Arguments:
  #  p: p outputs from Mark model, must be formatted with (format.p) first
  #
  # Return:
  #  plot of p with estimate and error bars representing LCI, UCI
  p %>% 
    mutate(
      Reach = factor(reach_num, levels = p$reach_num)
    ) %>% 
    ggplot(mapping = aes(x = reach_num, y = estimate)) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = lcl, ymax = ucl)) +
    ylab("Detection probability") +
    xlab("Reach") +
    scale_x_continuous(breaks = p$reach_num) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=.75),
    )
}

format.cum.surv <- function(cum_survival_all) {
  # Format cumulative survival outputs for plotting and table outputs
  #
  # Arguments:
  #  cum_survival_all: output df from Mark model cumulative survival
  #
  # Return:
  #  properly formatted df for phi outputs, now ready to plot
  
  cum_survival_all <- cum_survival_all %>% 
    add_column(
      GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], 
                length(studyIDs)),
      RKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], 
                length(studyIDs)),
      reach_num = rep(seq(0, (nrow(reach.meta.aggregate))-1, 1), length(studyIDs))
    ) %>% 
    left_join(
      reach.meta.aggregate %>%
        select(GEN, Region) %>% 
        distinct(),
      by = c("GEN")
    ) %>% 
    mutate(
      'Survival estimate (SE)' = paste0(cum.phi, " (", as.character(cum.phi.se), ")"),
      'Reach #' = reach_num
    ) %>% filter(
      GEN != "GoldenGateW"
    )

}

plot.cum.surv <- function(cum_survival_all, add_breaks, multiple, padding = 5.5) {
  # Plot cumulative survival outputs from Mark model
  #
  # Arguments:
  #  cum_survival_all: cumulative survival outputs from Mark model, 
  #     must be formatted with (format.cum.surv) first
  #  add_breaks: TRUE/FALSE whether to add vertical line breaks to represent
  #     regions
  #  multiple: TRUE/FALSE, whether there are multiple studyIDs or not
  #  padding: leftside plot margin, default set to 5.5 good for most, but can
  #     be adjusted of the xaxis label too long and gets cut off
  #
  # Return:
  #  plot of cumulative survival with estimate and error bars representing LCI, UCI
  
  if (multiple) {
    cum_survival_all %>% 
      mutate(
        GEN = factor(GEN, levels = reach.meta.aggregate$GEN,
                     labels = paste0(reach.meta.aggregate$GEN, " (", 
                                     reach.meta.aggregate$GenRKM, ")"))
      ) %>% 
      ggplot(mapping = aes(x = GEN, y = cum.phi, group = StudyID)) +
      geom_point(size = 2, aes(color = StudyID)) +
      geom_errorbar(mapping = aes(x= GEN, ymin = LCI, ymax = UCI, 
                                  color = StudyID),  width = .1) +
      geom_line(size = 0.7, aes(color = StudyID)) +
      {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
      ylab("Cumulative survival") +
      xlab("Site (River KM)") +
      scale_y_continuous(breaks = seq(0, 1, 0.1)) +
      # NEEDS TO BE FIXED FOR IF THERE ARE MORE THEN 2 STUDYIDS
      scale_color_manual(values=c("#007EFF", "#FF8100")) +
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
        legend.position = c(.1 ,.12)
      ) 
  } else {
    cum_survival_all %>% 
      mutate(
        GEN = factor(GEN, levels = reach.meta.aggregate$GEN,
                     labels = paste0(reach.meta.aggregate$GEN, " (", 
                                     reach.meta.aggregate$GenRKM, ")"))
      ) %>% 
      ggplot(mapping = aes(x = GEN, y = cum.phi)) +
      geom_point(size = 2) +
      geom_errorbar(mapping = aes(x= GEN, ymin = LCI, ymax = UCI),  width = .1) +
      geom_line(size = 0.7, group = 1) +
      {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
      ylab("Cumulative survival") +
      xlab("Site (River KM)") +
      scale_y_continuous(breaks = seq(0, 1, 0.1)) +
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin=unit(c(.5,.5,.5,1), "cm"),
      ) 
  }
}


#### Create INP for BC-Spring-2018  -----------------------------------

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("BC-Spring-2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>% 
  filter(
    GEN != "UpperButte_RST",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

#### BC-Spring-2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)
p <- format.p(outputs, multiple = F)

plot.p(p)

# Plot detection probability
p %>% 
  mutate(
    Reach = factor(reach_num, levels = p$reach_num)
  ) %>% 
  ggplot(mapping = aes(x = reach_num, y = estimate)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = lcl, ymax = ucl)) +
  ylab("Detection probability") +
  xlab("Reach") +
  scale_x_continuous(breaks = p$reach_num) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
  )


#### BC-Spring-2018 Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))


all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

ouputs <- get.mark.model(all.inp, standardized = T, multiple = F)

phi <- format_phi(outputs, multiple = F)

unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              ungroup() %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- get.region.breaks(phi)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F)

# Create folder for studyid, save plot
dir.create(file.path(getwd(), "Outputs", "BC-Spring-2018"), showWarnings = FALSE)  
ggsave2("./Outputs/BC-Spring-2018/BC-Spring-2018 Reach Survival per10km.png", width = 6, height = 5, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, "./Outputs/BC-Spring-2018/BC-Spring-2018 Reach Survival per10km.csv")


#### BC-Spring-2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- cum_survival_2018

cum_survival_all <- format.cum.surv(cum_survival_all, final = F)

# Plot
plot.cum.surv(cum_survival_all, add_breaks = F)

cum_survival_all <- format.cum.surv(cum_survival_all, final = T)

ggsave2("./Outputs/BC-Spring-2018/BC-Spring-2018 Cumulative Survival.png", width = 10, height = 6, dpi = 500)
write_csv(cum_survival_all, "./Outputs/BC-Spring-2018/BC-Spring-2018 Cumulative Survival.csv")


#### Create INP for ColemanLateFall_2018  -----------------------------------

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("ColemanLateFall_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>% 
  filter(
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

# Filter out poor receiver sites based on receiver notes, and detection probability
reach.meta <- reach.meta %>% 
  filter(
    !(GEN %in% c("AbvColusaBr", "Colusa BC3", "Freeport", "Blw_FR_GS2",
                 "RB_Elks", "BlwChinaBend", "GCID_blw", "Blw_FremontWeir",
                 "Blw_FRConf"))
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

#### ColemanLateFall_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp)

p <- format.p(outputs)

# Plot detection probability
p %>% 
  mutate(
    Reach = factor(reach_num, levels = p$reach_num)
  ) %>% 
  ggplot(mapping = aes(x = reach_num, y = estimate)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = lcl, ymax = ucl)) +
  ylab("Detection probability") +
  xlab("Reach") +
  geom_hline(yintercept = 0.7, linetype = 'dashed', color = 'red') +
  scale_x_continuous(breaks = p$reach_num) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
  )


# Look at detection probabilities in which estimates are < 0.7
p %>% 
  filter(estimate < 0.7) %>% 
  pull(reach_start)

reach.meta <- reach.meta %>% 
  filter(
    # GEN start that have poor estimates
    !GEN %in% c("Blw_Paynes_Ck", "GCID_abv", "Colusa AC2", "Colusa BC2",
                "BlwTisdale", "Blw_FremontWeir", "Blw_FRConf"),
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

outputs <- get.mark.model(all.inp)
p <- format.p(outputs)

#### ColemanLateFall_2018 Reach survival per 10km ----------------------------------------------------------------

replace_dict <- list(replace_with = list(c("Chipps")),
                     replace_list = list(c("ChippsE", "ChippsW")))

# Make final selection of receiver sites
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BattleCk_CNFH_Rel", "BattleCk4", "Blw_Paynes_Ck", "Abv_Altube1",
               "Mill_Ck_Conf", "Abv_WoodsonBr", "Blw_IrvineFinch", "BlwOrd",
               "ButteBr", "Colusa AC2", "Colusa BC2", "Blw_Knights_GS3", 
               "Blw_Elkhorn_GS1", "TowerBridge", "Hood", "ChippsE", "ChippsW",
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

all.process <- process.data(all.inp, model="CJS", begin.time=1,time.intervals = reach_length)

#again, ? make.design.data if you are curious about this call
all.ddl <- make.design.data(all.process)

rm(list=ls(pattern="p.t.x.y"))
rm(list=ls(pattern="Phi.t.x.y"))

# ## Now, set up basic model structures
#this is setting detection probability (p) to be a function of time x year, which is typically the best for simple survival analysis
p.t.x.y <- list(formula= ~time)

#This is setting the survival model to be a function of time (which is actually reach) x year
Phi.t.x.y <- list(formula= ~time)

#create a model list of the above models
cml = create.model.list("CJS")

# # Run mark.wrapper for all model structures (the combination of all phi and p models possible). Warning, this could take awhile if you have many models
model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) #,adjust=FALSE

outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real

cleanup(ask = FALSE)

# From outputs format to include the first half of rows which are phi estimates
# Extract StudyID from the rowname jumble
phi <- outputs %>%
  slice(1:(nrow(outputs) / 2)) %>%
  select(
    -c("fixed", "note")
  ) %>% 
  add_column(
    reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
    reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
    rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
    rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>%
      distinct(),
    by = c("reach_start" = "GEN")
  ) %>% 
  mutate(
    Reach = paste0(reach_start, " to \n", reach_end),
    RKM = paste0(rkm_start, " to ", rkm_end),
    Region = case_when(
      reach_end == "Blw_Paynes_Ck" ~ "Upper Sac R",
      reach_start == "ButteBr" ~ "Upper Sac R",
      reach_start == "Hood" ~ "Delta",
      reach_start %in% c("Chipps", "BeniciaW", "GoldenGateE", "GoldenGateW") ~ "Bay",
      TRUE ~ Region
    ),
    reach_num = 1:n()
  )

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- all_aggregated %>% 
  bind_rows() %>% 
  select(StudyID, FishID, GEN, GenRKM) %>% 
  distinct() %>% 
  group_by(StudyID, GEN, GenRKM) %>% 
  summarise(
    count = n()
  ) %>% 
  arrange(StudyID, desc(GenRKM))

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- phi %>% 
  group_by(Region) %>% 
  summarise(breaks = max(reach_num)) %>% 
  mutate(breaks = breaks + .5) %>% 
  arrange(breaks) %>% 
  slice(1:(n()-1)) %>% 
  pull(breaks)

# Plot survival probability
phi %>% 
  filter(reach_end != "GoldenGateW") %>% 
  mutate(
    Reach = factor(Reach, levels = phi$Reach)
  ) %>% 
  ggplot(mapping = aes(x = Reach, y = estimate)) +
  geom_point() +
  geom_errorbar(mapping = aes(x= Reach, ymin = lcl, ymax = ucl),  width = .1) +
  ylab("Survival per 10km") +
  xlab("Reach") +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave2("./Outputs/ColemanLateFall_2018/ColemanLateFall_2018 Reach Survival per10km.png", width = 10, height = 6, dpi = 500)


phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, UCI = ucl, reach_num, count) %>% 
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

write_csv(phi_table, "./Outputs/ColemanLateFall_2018/ColemanLateFall_2018 Reach Survival per10km.csv")


#### ColemanLateFall_2018 Cumulative Survival By Year----------------------------------------------------------------

cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cleanup(ask = FALSE)

cum_survival_all <- cum_survival_all %>% 
  add_column(
    survival_to = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], length(studyIDs)),
    survival_to_rkm = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], length(studyIDs)),
    reach_num = rep(seq(1, (nrow(reach.meta.aggregate)), 1), length(studyIDs))
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>% 
      distinct(),
    by = c("survival_to" = "GEN")
  ) %>% 
  mutate(
    Region = case_when(
      survival_to == "Blw_Paynes_Ck" ~ "Upper Sac R",
      survival_to == "ButteBr" ~ "Upper Sac R",
      survival_to == "Chipps" ~ "Lower Sac R",
      survival_to %in% c("BeniciaW", "GoldenGateE", "GoldenGateW") ~ "Bay",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("Battle Ck", "Upper Sac R", "Lower Sac R", "Delta", "Bay"))
  )

region_breaks <- cum_survival_all %>%
  group_by(Region) %>%
  summarise(breaks = max(reach_num)) %>%
  mutate(breaks = breaks + .5) %>%
  slice(1:(n()-1)) %>%
  pull(breaks)

# Plot
cum_survival_all %>% 
  filter(survival_to != "GoldenGateW") %>% 
  mutate(
    survival_to = factor(survival_to, levels = reach.meta.aggregate$GEN,
                         labels = paste0(reach.meta.aggregate$GEN, " (", reach.meta.aggregate$GenRKM, ")"))
  ) %>% 
  ggplot(mapping = aes(x = survival_to, y = cum.phi)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(x= survival_to, ymin = LCI, ymax = UCI),  width = .1) +
  geom_line(size = 0.7, group = 1) +
  ylab("Cumulative survival") +
  xlab("Site (River KM)") +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin=unit(c(.5,.5,.5,1), "cm"),
  ) 


ggsave2("./Outputs/ColemanLateFall_2018/ColemanLateFall_2018 Cumulative Survival.png", width = 10, height = 6, dpi = 500)


write_csv(cum_survival_all, "./Outputs/ColemanLateFall_2018/ColemanLateFall_2018 Cumulative Survival.csv")




#### Create INP for ColemanLateFall_2020  -----------------------------------

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("ColemanLateFall_2020")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- all_detections %>% 
  bind_rows() %>% 
  distinct(GEN, GenRKM, GenLat, GenLon, Region) %>% 
  # Necessary because detections files shows differing RKM, Lat, Lon for some GEN sometimes
  group_by(GEN) %>% 
  summarise(
    GenRKM = mean(GenRKM),
    GenLat = mean(GenLat),
    GenLon = mean(GenLon),
    Region = first(Region)
  ) %>% 
  arrange(desc(GenRKM))
  


# Adjust sites, remove site above release, and Delta sites except Chipps
reach.meta <- reach.meta %>% 
  filter(
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

# Filter out poor receiver sites based on receiver notes, and detection probability
reach.meta <- reach.meta %>% 
  filter(
    !(GEN %in% c("AbvColusaBr", "Colusa BC3", "Freeport", "Blw_FR_GS2",
                 "RB_Elks", "BlwChinaBend", "GCID_blw", "Blw_FremontWeir",
                 "Blw_FRConf"))
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()


#### ColemanLateFall_2020 Detection probability ---------------------------------------------------
all.process <- process.data(all.inp, model="CJS", begin.time=1)
all.ddl <- make.design.data(all.process)

rm(list=ls(pattern="p.t.x.y"))
rm(list=ls(pattern="Phi.t.x.y"))

p.t.x.y <- list(formula= ~time)
Phi.t.x.y <- list(formula= ~time)

cml = create.model.list("CJS")

model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) 

outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real

cleanup(ask = FALSE)

p <- outputs %>%
  slice(
    ((nrow(outputs)/2) +1):nrow(outputs)
  ) %>%
  select(
    -c("fixed", "note")
  ) %>%
  add_column(
    reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
    reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
    rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
    rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
  ) %>% 
  mutate(
    Reach = paste0(reach_start, " to \n", reach_end),
    RKM = paste0(rkm_start, " to ", rkm_end),
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>%
      distinct(),
    by = c("reach_start" = "GEN")
  ) %>% 
  mutate(reach_num = 1:n())

# Plot detection probability
p %>% 
  mutate(
    Reach = factor(reach_num, levels = p$reach_num)
  ) %>% 
  # mutate(reach_num = factor(reach_num, levels(factor(seq(1, 14))))) %>% 
  ggplot(mapping = aes(x = reach_num, y = estimate)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = lcl, ymax = ucl)) +
  ylab("Detection probability") +
  xlab("Reach") +
  geom_hline(yintercept = 0.7, linetype = 'dashed', color = 'red') +
  # geom_vline(xintercept = region_breaks, linetype = "dotted") +
  # scale_y_continuous(breaks = seq(0.4, 1, .1), limits = c(0.4,1)) +
  #scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = p$reach_num) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    # axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Look at detection probabilities in which estimates are < 0.7
p %>% 
  filter(estimate < 0.7) %>% 
  pull(reach_start)

reach.meta <- reach.meta %>% 
  filter(
    # GEN start that have poor estimates
    !GEN %in% c("Blw_Paynes_Ck", "GCID_abv", "Colusa AC2", "Colusa BC2",
                "BlwTisdale", "Blw_FremontWeir", "Blw_FRConf"),
    GEN != "LSNFH",
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

all.process <- process.data(all.inp, model="CJS", begin.time=1)
all.ddl <- make.design.data(all.process)

rm(list=ls(pattern="p.t.x.y"))
rm(list=ls(pattern="Phi.t.x.y"))

p.t.x.y <- list(formula= ~time)
Phi.t.x.y <- list(formula= ~time)

cml = create.model.list("CJS")

model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) 

outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real

cleanup(ask = FALSE)

p <- outputs %>%
  slice(
    ((nrow(outputs)/2) +1):nrow(outputs)
  ) %>%
  select(
    -c("fixed", "note")
  ) %>%
  add_column(
    reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
    reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
    rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
    rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
  ) %>% 
  mutate(
    Reach = paste0(reach_start, " to \n", reach_end),
    RKM = paste0(rkm_start, " to ", rkm_end),
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>%
      distinct(),
    by = c("reach_start" = "GEN")
  ) %>% 
  mutate(reach_num = 1:n())



#### ColemanLateFall_2020 Reach survival per 10km ----------------------------------------------------------------

replace_dict <- list(replace_with = list(c("Chipps")),
                     replace_list = list(c("ChippsE", "ChippsW")))

# Make final selection of receiver sites
reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("BattleCk_CNFH_Rel", "BattleCk4", "Blw_Paynes_Ck", "Abv_Altube1",
               "Mill_Ck_Conf", "Abv_WoodsonBr", "Blw_IrvineFinch", "BlwOrd",
               "ButteBr", "Colusa AC2", "Colusa BC2", "Blw_Knights_GS3", 
               "Blw_Elkhorn_GS1", "TowerBridge", "Hood", "ChippsE", "ChippsW",
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Get river KM
KM <- reach.meta.aggregate$GenRKM

## Calculate the reach lengths. Here I devided reach lengths by 10 so that my survival estimates later will be survival per 10km
reach_length <- abs(diff(KM))/10

all.process <- process.data(all.inp, model="CJS", begin.time=1,time.intervals = reach_length)

#again, ? make.design.data if you are curious about this call
all.ddl <- make.design.data(all.process)

rm(list=ls(pattern="p.t.x.y"))
rm(list=ls(pattern="Phi.t.x.y"))

# ## Now, set up basic model structures
#this is setting detection probability (p) to be a function of time x year, which is typically the best for simple survival analysis
p.t.x.y <- list(formula= ~time)

#This is setting the survival model to be a function of time (which is actually reach) x year
Phi.t.x.y <- list(formula= ~time)

#create a model list of the above models
cml = create.model.list("CJS")

# # Run mark.wrapper for all model structures (the combination of all phi and p models possible). Warning, this could take awhile if you have many models
model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) #,adjust=FALSE

outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real

cleanup(ask = FALSE)

# From outputs format to include the first half of rows which are phi estimates
# Extract StudyID from the rowname jumble
phi <- outputs %>%
  slice(1:(nrow(outputs) / 2)) %>%
  select(
    -c("fixed", "note")
  ) %>% 
  add_column(
    reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
    reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
    rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
    rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>%
      distinct(),
    by = c("reach_start" = "GEN")
  ) %>% 
  mutate(
    Reach = paste0(reach_start, " to \n", reach_end),
    RKM = paste0(rkm_start, " to ", rkm_end),
    Region = case_when(
      reach_end == "Blw_Paynes_Ck" ~ "Upper Sac R",
      reach_start == "ButteBr" ~ "Upper Sac R",
      reach_start == "Hood" ~ "Delta",
      reach_start %in% c("Chipps", "BeniciaW", "GoldenGateE", "GoldenGateW") ~ "Bay",
      TRUE ~ Region
    ),
    reach_num = 1:n()
  )

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- all_aggregated %>% 
  bind_rows() %>% 
  select(StudyID, FishID, GEN, GenRKM) %>% 
  distinct() %>% 
  group_by(StudyID, GEN, GenRKM) %>% 
  summarise(
    count = n()
  ) %>% 
  arrange(StudyID, desc(GenRKM))

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- phi %>% 
  group_by(Region) %>% 
  summarise(breaks = max(reach_num)) %>% 
  mutate(breaks = breaks + .5) %>% 
  arrange(breaks) %>% 
  slice(1:(n()-1)) %>% 
  pull(breaks)

# Plot survival probability
phi %>% 
  filter(reach_end != "GoldenGateW") %>% 
  mutate(
    Reach = factor(Reach, levels = phi$Reach)
  ) %>% 
  ggplot(mapping = aes(x = Reach, y = estimate)) +
  geom_point() +
  geom_errorbar(mapping = aes(x= Reach, ymin = lcl, ymax = ucl),  width = .1) +
  ylab("Survival per 10km") +
  xlab("Reach") +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.75),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave2("./Outputs/ColemanLateFall_2020/ColemanLateFall_2020 Reach Survival per10km.png", width = 10, height = 6, dpi = 500)


phi_table <- phi %>% 
  select(StudyID, Reach, RKM, Region, Estimate = estimate, SE = se, LCI = lcl, UCI = ucl, reach_num, count) %>% 
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

write_csv(phi_table, "./Outputs/ColemanLateFall_2020/ColemanLateFall_2020 Reach Survival per10km.csv")


#### ColemanLateFall_2020 Cumulative Survival By Year----------------------------------------------------------------

cum_survival_2020 <- get_cum_survival((all.inp), T)
cum_survival_2020

cum_survival_all <- bind_rows(cum_survival_2020)

cleanup(ask = FALSE)

cum_survival_all <- cum_survival_all %>% 
  add_column(
    survival_to = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], length(studyIDs)),
    survival_to_rkm = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], length(studyIDs)),
    reach_num = rep(seq(1, (nrow(reach.meta.aggregate)), 1), length(studyIDs))
  ) %>% 
  left_join(
    reach.meta.aggregate %>%
      select(GEN, Region) %>% 
      distinct(),
    by = c("survival_to" = "GEN")
  ) %>% 
  mutate(
    Region = case_when(
      survival_to == "Blw_Paynes_Ck" ~ "Upper Sac R",
      survival_to == "ButteBr" ~ "Upper Sac R",
      survival_to == "Chipps" ~ "Lower Sac R",
      survival_to %in% c("BeniciaW", "GoldenGateE", "GoldenGateW") ~ "Bay",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("Battle Ck", "Upper Sac R", "Lower Sac R", "Delta", "Bay"))
  )

region_breaks <- cum_survival_all %>%
  group_by(Region) %>%
  summarise(breaks = max(reach_num)) %>%
  mutate(breaks = breaks + .5) %>%
  slice(1:(n()-1)) %>%
  pull(breaks)

# Plot
cum_survival_all %>% 
  filter(survival_to != "GoldenGateW") %>% 
  mutate(
    survival_to = factor(survival_to, levels = reach.meta.aggregate$GEN,
                         labels = paste0(reach.meta.aggregate$GEN, " (", reach.meta.aggregate$GenRKM, ")"))
  ) %>% 
  ggplot(mapping = aes(x = survival_to, y = cum.phi)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(x= survival_to, ymin = LCI, ymax = UCI),  width = .1) +
  geom_line(size = 0.7, group = 1) +
  ylab("Cumulative survival") +
  xlab("Site (River KM)") +
  geom_vline(xintercept = region_breaks, linetype = "dotted") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin=unit(c(.5,.5,.5,1), "cm"),
  ) 


ggsave2("./Outputs/ColemanLateFall_2020/ColemanLateFall_2020 Cumulative Survival.png", width = 10, height = 6, dpi = 500)


write_csv(cum_survival_all, "./Outputs/ColemanLateFall_2020/ColemanLateFall_2020 Cumulative Survival.csv")



#### Create INP for DeerCk_Wild_2018 -----------------------------------
name <- "DeerCk_Wild_2018"

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("DeerCk_Wild_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

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

#### DeerCk_Wild_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

# Identify and remove sites with low p <0.7
reach.meta <- reach.meta %>% 
  filter(
    !(GEN %in% c("AbvTisdale", "AbvColusaBr", "GCID_blw"))
  )
all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Then rerun 
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### DeerCk_Wild_2018 Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("DeerCk_RST", "DeerCk3", "Abv_WoodsonBr", "GCID_abv", "Blw_IrvineFinch", "BlwOrd", "ButteBr")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- 2.5

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
              width = 9, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                           " Reach Survival per10km.csv")))


#### DeerCk_Wild_2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
region_breaks <- 2.5
plot.cum.surv(cum_survival_all, T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 9, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv")))

#### DeerCk_Wild_2018 Region----------------------------------------------------------------
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("DeerCk_RST", "Abv_WoodsonBr", "ButteBr")
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)

all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- 1.5

# Plot survival probability
phi_per10km$Region <- c("Deer Creek", "Upper Sacramento")
plot.phi(phi_per10km, type = "Region", add_breaks = T, ylabel = "Survival per 10km",
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 4, height = 4, dpi = 500)

# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = F)
phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 4, height = 4, dpi = 500)


phi_table <- make.phi.table(phi_region)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Region Survival.csv")))


#### Create INP for Mok_Fall_2018 -----------------------------------
name <- "Mok_Fall_2018"

studyIDs <- c("Mok_Fall_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

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

#### Mok_Fall_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### Mok_Fall_2018 Reach survival per 10km ----------------------------------------------------------------
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

# region_breaks <- get.region.breaks(phi)

# Plot survival probability
plot.phi(phi, type = "Reach",  add_breaks = F, ylabel = "Survival per 10km",
         xlabel = "Reach", multiple = F)


ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 6, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### Mok_Fall_2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
plot.cum.surv(cum_survival_all %>% 
                filter(reach_num != 5), F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 6, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv")))

### Create INP for Nimbus_Fall_2018 -----------------------------------
name <- "Nimbus_Fall_2018"

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("Nimbus_Fall_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)


reach.meta <- reach.meta %>% 
  filter(
    GEN %in% c("AR_Sunrise_Ramp_Rel", "TowerBridge", "SacTrawl1",
               "SacTrawl2", "Hood", "ChippsE", "ChippsW", "BeniciaE", 
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  )

leaflet(data = reach.meta) %>% addTiles() %>%
  addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
             label = ~as.character(GEN),
             labelOptions = labelOptions(noHide = T, textOnly = TRUE))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

#### Nimbus_Fall_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

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
outputs <- get.mark.model(all.inp, F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### Nimbus_Fall_2018 Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("AR_Sunrise_Ramp_Rel", "TowerBridge", "SacTrawl1",
               "SacTrawl2", "Hood", "ChippsE", "ChippsW", "BeniciaE", 
               "BeniciaW", "GoldenGateE", "GoldenGateW")
  ) 


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
  count()

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count = n),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5, 4.5, 5.5)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks =  T, ylabel = "Survival rate per 10km", 
        xlabel = "Reach",  multiple = F)


ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### Nimbus_Fall_2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
plot.cum.surv(cum_survival_all, T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

#### Nimbus_Fall_2018 Region survival  ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("AR_Sunrise_Ramp_Rel", "TowerBridge", "Hood", "ChippsE", 
               "ChippsW","GoldenGateE", "GoldenGateW")
  ) 

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = F) %>% 
  filter(reach_num != 5) %>% 
  mutate(Region = c("American River", "Lower Sacramento", "Delta", "Bay"))

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
 summarise(count = n())

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot.phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = F)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 6, height = 6, dpi = 500)

phi_per10km_table <- make.phi.table(phi_region, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Region Survival per 10km.csv")))


# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = F) %>% 
  filter(reach_num != 5) %>% 
  mutate(Region = c("American River", "Lower Sacramento", "Delta", "Bay"))

phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 6, height = 6, dpi = 500)


phi_region_table <- make.phi.table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Region Survival.csv")))

#### Create INP for RBDD_2018 -----------------------------------
name <- "RBDD_2018"

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("RBDD_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

leaflet(data = reach.meta) %>% addTiles() %>%
  addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
             label = ~as.character(GEN),
             labelOptions = labelOptions(noHide = T, textOnly = TRUE))


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

#### RBDD_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, F)
cleanup(ask = FALSE)

p <- format.p(outputs)

# Plot detection probability
plot.p(p)

# Identify and remove sites with low p <0.7
reach.meta <- reach.meta %>% 
  filter(
    !(GEN %in% c("Blw_FremontWeir"))
  )

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Then rerun 
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### RBDD_2018 Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

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

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(7.5)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = F)

# Create folder for studyid, save plot
ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### RBDD_2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
region_breaks <- 8.5
plot.cum.surv(cum_survival_all, T)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"), 
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Cumulative Survival.csv")))

#### RBDD_2018 Region survival  ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("RBDD_Rel", "ButteBr", "Hood")
  ) 

all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = F) %>% 
  mutate(Region = c("Upper Sacramento", "Lower Sacramento"))

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
  summarise(count = n())

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5)

# Plot survival probability
plot.phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = F)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 6, height = 6, dpi = 500)

phi_per10km_table <- make.phi.table(phi_per10km, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                           " Region Survival per 10km.csv")))

# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = F)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = F) %>% 
  mutate(Region = c("Upper Sacramento", "Lower Sacramento"))

phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 6, height = 6, dpi = 500)


phi_region_table <- make.phi.table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))


##### Create INP for RBDD_WR_2018 -----------------------------------
name <- "RBDD_WR_2018"

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("RBDD_WR_2018")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    !Region %in% c("North Delta", "East Delta", "West Delta", "Yolo Bypass") |
      GEN %in% c("ChippsE", "ChippsW"),
    GEN != "I80-50_Br"
  )

leaflet(data = reach.meta) %>% addTiles() %>%
  addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
             label = ~as.character(GEN),
             labelOptions = labelOptions(noHide = T, textOnly = TRUE))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

#### RBDD_WR_2018 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, F, multiple = F)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = F)

# Plot detection probability
plot.p(p)

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
outputs <- get.mark.model(all.inp, F)
cleanup(ask = FALSE)

p <- format.p(outputs)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### RBDD_WR_2018 Reach survival per 10km ----------------------------------------------------------------
replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

reach.meta <- reach.meta %>%
  filter(
    GEN %in% c("RBDD_Rel", "Blw_Salt", "BankRobber", "Lwr_Ant_Crk", "Mill_Ck_Conf",
               "Abv_WoodsonBr", "Blw_Woodson", "GCID_abv", "GCID_blw", 
               "Blw_IrvineFinch")
  ) 


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
  addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
             label = ~as.character(GEN),
             labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(all.inp, standardized = T, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count) %>% 
              distinct(),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))


# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = F)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))

#### RBDD_WR_2018 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp), add_release = T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
plot.cum.surv(cum_survival_all, add_breaks = F)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))


### Create INP for SB_Spring_2018_2019 -----------------------------------
name <- "SB_Spring_2018_2019"

# Build the key/pair for GEN replacement
replace_dict <- list(replace_with = list(c()),
                     replace_list = list(c()))

studyIDs <- c("SB_Spring_2018", "SB_Spring_2019")

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

# Multi-year studyID so use common sites list
common_sites <- read_csv("SB_Spring_2018_2019_ReceiverSites.csv") %>% 
  filter(
    SB_Spring_2018 == "PRESENT",
    SB_Spring_2019 == "PRESENT"
  )

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

#### SB_Spring_2018_2019 Detection probability ---------------------------------------------------
outputs <- get.mark.model(all.inp, F, multiple = T)
cleanup(ask = FALSE)

p <- format.p(outputs, multiple = T)

# Plot detection probability
plot.p(p)

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
outputs <- get.mark.model(all.inp, F)
cleanup(ask = FALSE)

p <- format.p(outputs)

# Plot detection probability
plot.p(p)

dir.create(paste0("./Outputs/", name), showWarnings = FALSE)  
write_csv(p, paste0(paste0("./Outputs/", name, "/", name, 
                           " detection probability.csv")))

#### SB_Spring_2018_2019 Reach survival per 10km ----------------------------------------------------------------
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

outputs <- get.mark.model(all.inp, T, multiple = T)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = T)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count, StudyID),
            by = c("reach_end" = "GEN", "StudyID")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(4.5, 6.5, 7.5)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = T, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = T, padding = 25)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### SB_Spring_2018_2019 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "SB_Spring_2018")), T)
cum_survival_2018

cum_survival_2019 <- get_cum_survival((all.inp %>% 
                                         filter(StudyID == "SB_Spring_2019")), T)
cum_survival_2019

cum_survival_all <- bind_rows(cum_survival_2018, cum_survival_2019)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
region_breaks <- c(5.5, 7.5, 8.5)
plot.cum.surv(cum_survival_all, add_breaks = T, multiple = T, padding = 50)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    StudyID, 'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

#### SB_Spring_2018_2019 Region survival  ----------------------------------------------------------------
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

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_per10km$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
  summarise(count = n())

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot.phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = T)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 8, height = 6, dpi = 500)

phi_per10km_table <- make.phi.table(phi_per10km, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                           " Region Survival per 10km.csv")))

# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_region$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 8, height = 6, dpi = 500)


phi_region_table <- make.phi.table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))




### FOR 2019 REPORT -----------------------------------
name <- "Winter_H_2019"

studyIDs <- c("Winter_H_2019")

replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    # ColemanLateFall_2019
    # GEN %in% c("BattleCk_CNFH_Rel", "Abv_Altube1", "Mill_Ck_Conf", "Abv_WoodsonBr",
    #            "Blw_IrvineFinch", "BlwOrd", "ButteBr", "Colusa AC2", "Colusa BC2",
    #            "Blw_Knights_GS3","Blw_Elkhorn_GS1", "TowerBridge", "ChippsE", 
    #            "ChippsW","GoldenGateE","GoldenGateW")
    # SB_Spring_2019
    # GEN %in% c("SutterBypass_Weir2_RST_Rel", "Butte2", "Butte3", "Butte6",
    #            "I80-50_Br", "ChippsE", "ChippsW",
    #            "GoldenGateE","GoldenGateW")
    # Winter_H_2019
    GEN %in% c("Caldwell_Park", "Colusa AC3", "KnightsBlwRST",
               "ChippsE", "ChippsW", "GoldenGateE",
    "GoldenGateW") #"BeniciaE", "BeniciaW",
    # #CNFH_FMR
    # GEN %in% c("RBDD_Rel", "Blw_Salt", "Abv_WoodsonBr", "Blw_IrvineFinch",
    #            "BlwOrd", "ButteBr", "Colusa AC3", "Colusa BC2", "AbvTisdale", 
    #            "BlwChinaBend", "Knights_RST", "I80-50_Br", "ChippsE", "ChippsW",
    #            "GoldenGateE", "GoldenGateW")
  )

# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()



#### SB_Spring_2018_2019 Reach survival per 10km ----------------------------------------------------------------


# reach.meta <- reach.meta %>%
#   filter(
#     GEN %in% c("Caldwell_Park", "Blw_Cypress", "Blw_ClearCr", "Blw_Salt",
#                "GCID_abv", "Colusa AC3", "Colusa BC3", "BlwChinaBend", "KnightsBlwRST", 
#                "TowerBridge", "ChippsE", "ChippsW", "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
#   )



# leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(all.inp, F, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end"= "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(4.5, 6.5, 7.5)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = F, padding = 25)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### SB_Spring_2018_2019 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival(all.inp, T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
region_breaks <- c()
plot.cum.surv(cum_survival_all, add_breaks = F, multiple = F, padding = 50)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    StudyID, 'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

#### SB_Spring_2018_2019 Region survival  ----------------------------------------------------------------
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

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_per10km$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
  summarise(count = n())

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot.phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = T)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 8, height = 6, dpi = 500)

phi_per10km_table <- make.phi.table(phi_per10km, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                           " Region Survival per 10km.csv")))

# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_region$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 8, height = 6, dpi = 500)


phi_region_table <- make.phi.table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))


### FOR 2019 REPORT FR Spring -----------------------------------
name <- "FR_Spring_2019"

studyIDs <- c("FR_Spring_2019")

replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# gridley <- all_detections %>% 
#   bind_rows() %>% 
#   filter(Rel_loc == "FR_Gridley_Rel") %>% 
#   list()
# 
# boyds <- all_detections %>% 
#   bind_rows() %>% 
#   filter(Rel_loc == "FR_Boyds_Rel") %>% 
#   list()
  

# Get list of all receiver GEN
reach.meta <- get.receiver.GEN(all_detections)

reach.meta <- reach.meta %>% 
  filter(
    # Boyds
    GEN %in% c("BoydsPump", "BC_Beach",
               "Blw_FRConf", "I80-50_Br", "ChippsE", "ChippsW",
               "GoldenGateE", "GoldenGateW")
    # Gridley
    # GEN %in% c("FR_Gridley_Rel", "SunsetPumps", "BoydsPump", "BC_Beach",
    #            "Blw_FRConf", "I80-50_Br", "ChippsE", "ChippsW",
    #            "GoldenGateE", "GoldenGateW")
    # ColemanLateFall_2019
    # GEN %in% c("BattleCk_CNFH_Rel", "Abv_Altube1", "Mill_Ck_Conf", "Abv_WoodsonBr",
    #            "Blw_IrvineFinch", "BlwOrd", "ButteBr", "Colusa AC2", "Colusa BC2",
    #            "Blw_Knights_GS3","Blw_Elkhorn_GS1", "TowerBridge", "ChippsE", 
    #            "ChippsW","GoldenGateE","GoldenGateW")
    # SB_Spring_2019
    # GEN %in% c("SutterBypass_Weir2_RST_Rel", "Butte2", "Butte3", "Butte6",
    #            "I80-50_Br", "ChippsE", "ChippsW",
    #            "GoldenGateE","GoldenGateW")
    # Winter_H_2019
    # GEN %in% c("Caldwell_Park", "Blw_Cypress", "Blw_ClearCr", "Blw_Salt",
    #            "GCID_abv", "Colusa AC3", "Colusa BC3", "BlwChinaBend", "KnightsBlwRST",
    #            "ChippsE", "ChippsW", "GoldenGateE",
    #"GoldenGateW") #"BeniciaE", "BeniciaW", 
# #CNFH_FMR
# GEN %in% c("RBDD_Rel", "Blw_Salt", "Abv_WoodsonBr", "Blw_IrvineFinch",
#            "BlwOrd", "ButteBr", "Colusa AC3", "Colusa BC2", "AbvTisdale", 
#            "BlwChinaBend", "Knights_RST", "I80-50_Br", "ChippsE", "ChippsW",
#            "GoldenGateE", "GoldenGateW")
)

# leaflet(data = reach.meta) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))


all_aggregated <- lapply(all_detections, aggregate_GEN)
all_EH <- lapply(all_aggregated, make_EH)
all.inp <- pmap(list(all_detections,all_EH), create_inp) %>% 
  bind_rows()

gridley <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Gridley_Rel")

boyds <- all.inp %>% 
  left_join(TaggedFish %>% 
              select(FishID = fish_id, release_location)) %>% 
  filter(release_location == "FR_Boyds_Rel")


#### SB_Spring_2018_2019 Reach survival per 10km ----------------------------------------------------------------


# reach.meta <- reach.meta %>%
#   filter(
#     GEN %in% c("Caldwell_Park", "Blw_Cypress", "Blw_ClearCr", "Blw_Salt",
#                "GCID_abv", "Colusa AC3", "Colusa BC3", "BlwChinaBend", "KnightsBlwRST", 
#                "TowerBridge", "ChippsE", "ChippsW", "BeniciaE", "BeniciaW", "GoldenGateE", "GoldenGateW")
#   )



# leaflet(data = reach.meta.aggregate) %>% addTiles() %>%
#   addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
#              label = ~as.character(GEN),
#              labelOptions = labelOptions(noHide = T, textOnly = TRUE))

KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10

outputs <- get.mark.model(boyds, F, multiple = F)
cleanup(ask = F)

phi <- format_phi(outputs, multiple = F)

# Identify number of fish detected at each GEN, use that to filter out estimates for fish 
# no longer seen
unique_detects <- get.unique.detects(all_aggregated)

# Add in counts to phi
phi <- phi %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end"= "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(4.5, 6.5, 7.5)

# Plot survival probability
plot.phi(phi, type = "Reach", add_breaks = F, ylabel = "Survival rate per 10km",
         xlabel = "Reach", multiple = F, padding = 25)

ggsave(paste0("./Outputs/", name, "/", name, " Reach Survival per10km.png"),
       width = 10, height = 6, dpi = 500)

phi_table <- make.phi.table(phi)
write_csv(phi_table, paste0(paste0("./Outputs/", name, "/", name, 
                                   " Reach Survival per10km.csv")))


#### SB_Spring_2018_2019 Cumulative Survival By Year----------------------------------------------------------------
cum_survival_2018 <- get_cum_survival(all.inp, T)
cum_survival_2018

cum_survival_all <- bind_rows(cum_survival_2018)

cum_survival_all <- format.cum.surv(cum_survival_all)

# Plot
region_breaks <- c()
plot.cum.surv(cum_survival_all, add_breaks = F, multiple = F, padding = 50)

ggsave(paste0("./Outputs/", name, "/", name, " Cumulative Survival.png"),
       width = 10, height = 6, dpi = 500)

cum_survival_all <- cum_survival_all %>% 
  select(
    StudyID, 'Reach #', GEN, RKM, Region, 'Survival estimate (SE)', LCI, UCI
  )

write_csv(cum_survival_all, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Cumulative Survival.csv")))

#### SB_Spring_2018_2019 Region survival  ----------------------------------------------------------------
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

# Remove Mare Island Release group fish
all.inp <- all.inp %>% 
  left_join(
    TaggedFish %>% 
      select(fish_id, release_location),
    by = c("FishID" = "fish_id")
  ) %>% 
  filter(release_location != "Mare_Island_Rel")

# Per 10km analysis
KM <- reach.meta.aggregate$GenRKM
reach_length <- abs(diff(KM))/10
outputs <- get.mark.model(all.inp, standardized = T, multiple = T)
cleanup(ask = F)

phi_per10km <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_per10km$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

unique_detects <- all_aggregated %>% bind_rows() %>% 
  select(FishID, GEN) %>% 
  distinct() %>% 
  left_join(
    TaggedFish %>% select(
      FishID = fish_id, release_location
    )
  ) %>% 
  filter(release_location != "Mare_Island_Rel") %>% 
  group_by(GEN) %>% 
  summarise(count = n())

# Add in counts to phi
phi_per10km <- phi_per10km %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

region_breaks <- c(1.5, 2.5, 3.5)

# Plot survival probability
plot.phi(phi_per10km, add_breaks = T, type = "Region", xlabel = "Region", 
         ylabel = "Survival rate per 10km", multiple = T)


ggsave(paste0("./Outputs/", name, "/", name, " Region Survival per10km.png"),
       width = 8, height = 6, dpi = 500)

phi_per10km_table <- make.phi.table(phi_per10km, standardized = T)
write_csv(phi_per10km_table, paste0(paste0("./Outputs/", name, "/", name, 
                                           " Region Survival per 10km.csv")))

# Per region analysis
outputs <- get.mark.model(all.inp, standardized = F, multiple = T)
cleanup(ask = F)

phi_region <- format_phi(outputs, multiple = T) %>% 
  filter(reach_num != 5)

phi_region$Region <- rep(c("Butte Creek", "Lower Sacramento", "Delta", "Bay"), 2)

phi_region <- phi_region %>% 
  left_join(unique_detects %>% 
              select(GEN, count),
            by = c("reach_end" = "GEN")) %>% 
  mutate(count = ifelse(is.na(count), 0, count))

plot.phi(phi_region, type = "Region", add_breaks = T, ylabel = "Survival per region",
         xlabel = "Region", multiple = T)

ggsave(paste0("./Outputs/", name, "/", name, " Region Survival.png"),
       width = 8, height = 6, dpi = 500)


phi_region_table <- make.phi.table(phi_region, standardized = F)
write_csv(phi_region_table, paste0(paste0("./Outputs/", name, "/", name, 
                                          " Region Survival.csv")))

### FOR 2019 REPORT FR Spring -----------------------------------
name <- "Winter_H_2019"

studyIDs <- c("Winter_H_2019")

replace_dict <- list(replace_with = list(c("Chipps"),
                                         c("SacTrawl"),
                                         c("Benicia")),
                     replace_list = list(c("ChippsE", "ChippsW"),
                                         c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW")))

# Retreive ERDDAP data 
all_detections <- lapply(studyIDs, get_detections)

# Get list of all receiver GEN
reach.meta <- all_detections %>% 
  bind_rows() %>% 
  distinct(GEN, GenRKM, Region, GenLat, GenLon) %>% 
  # Necessary because detections files shows differing RKM, Lat, Lon for some GEN sometimes
  group_by(GEN) %>% 
  summarise(
    GenRKM = mean(GenRKM),
    GenLat = mean(GenLat),
    GenLon = mean(GenLon),
    Region = first(Region)
  ) %>% 
  filter(
    Region %in% c("Upper Sac R", "Lower Sac R", "Carquinez Strait",
                  "SF Bay") | GEN %in% c("ChippsE", "ChippsW"),
    GEN != "Caldwell_Park"
  ) %>% 
  arrange(desc(GenRKM))

replace_dict <- list(replace_with = list(c("SacTrawl"),
                                         c("Benicia"),
                                         c("Chipps"),
                                         c("Abv_Yolo"),
                                         c("YB_LibertyIsBase"),
                                         c("RioVista")),
                     replace_list = list(c("SacTrawl1", "SacTrawl2"),
                                         c("BeniciaE", "BeniciaW"),
                                         c("ChippsE", "ChippsW"),
                                         c("Abv_FremontWeir", "Blw_Knights_GS3"),
                                         c("YB_LibertyIsBase1", "YB_LibertyIsBase2"),
                                         c("RioVistaUS", "RioVistaDS")))

knights_fish <- detections %>% 
  filter(GEN %in% c("Abv_FremontWeir", "Blw_Knights_GS3")) %>% 
  mutate(date = as.Date(time)) %>% 
  group_by(FishID, GEN) %>% 
  summarise(min_date = min(date))

# knights_fish <- detections %>% 
#   filter(GEN %in% c("Knights_RST", "KnightsBlwRST", "Abv_FremontWeir")) %>%
#   mutate(date = as.Date(time)) %>% 
#   group_by(FishID) %>% 
#   summarise(min_date = min(date))

# Yolo / Sac Movement -----------------------------------------------------
detections <- all_detections[[1]]

yolo_fish <- detections %>% 
  filter(
    FishID %in% knights_fish$FishID
  ) %>% 
  mutate(
    # If a fish is found in this list of GEN give it a 1 otherwise 0
    remove = ifelse((GEN %in%  c( "Blw_FremontWeir", "Butte6", "Blw_FRConf", "Blw_FR_GS2", "Blw_Elkhorn_GS1", "TowerBridge", "I80-50_Br")), 1, 0)
  ) %>% 
  group_by(FishID) %>% 
  summarise(
    remove = sum(remove)
  ) %>% 
  mutate(
    direction = ifelse(remove == 0, "Yolo", "Sac")
  )

sac_fish <- detections %>% 
  filter(FishID %in% yolo_fish$FishID[yolo_fish$direction == "Sac"])

yolo_fish <- detections %>% 
  filter(FishID %in% yolo_fish$FishID[yolo_fish$direction == "Yolo"])

# Aggregate sites for both groups
sac_fish <- lapply(list(sac_fish), aggregate_GEN) %>% 
  bind_rows() %>% 
  mutate(Type = "Sacramento River")

yolo_fish <- lapply(list(yolo_fish), aggregate_GEN) %>% 
  bind_rows() %>% 
  mutate(Type = "Yolo Bypass")


calc_movement <- function(detections) {
  min_detects <- detections %>% 
    filter(GEN %in% c("Abv_Yolo", "Chipps")) %>% 
    group_by(Type, FishID, GEN, GenRKM) %>% 
    # Get first time a fish was detected at either of those GEN
    summarise(time = min(time)) %>% 
    group_by(Type, FishID) %>% 
    # Filter for those that have both sites cause we need both to calc movement
    filter(n() == 2) %>% 
    summarise(
      min_time = min(time),
      max_time = max(time),
      min_rkm = min(GenRKM),
      max_rkm = max(GenRKM),
    ) %>% 
    mutate(
      km_day = ifelse(Type == "Yolo Bypass", (102.3 / as.numeric(difftime(max_time, min_time, "days"))),  
                      ((max_rkm - min_rkm) / as.numeric(difftime(max_time, min_time, "days"))))
    ) 
  # %>% 
  #   group_by(Type) %>% 
  #   summarise(mean_km_day = mean(km_day),
  #             sd = sd(km_day),
  #             count = n())
  
}

yolo_movement <- calc_movement(yolo_fish)
sac_movement <- calc_movement(sac_fish)

movement <- bind_rows(yolo_movement, sac_movement)

ggplot(data = movement, mapping = aes(x = Type, y = km_day)) +
  geom_boxplot(width = .4) +
  xlab("Route") +
  ylab("Speed (km/day)") +
  theme_classic()

ggsave2("./Outputs/Section 3.1.2.9 Winter_H_2019 Yolo Bypass Movement.png", width = 5, height = 4, dpi = 500)

