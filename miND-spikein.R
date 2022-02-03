#' miND-spikein.R - Copyright TAmiRNA GmbH 2022
#' Author: Andreas B. Diendorfer
#'
#' This file processes miND metadata and read counts (RC) files to calculate
#' miRNA concentrations based on the spikeins.
#'
#' Files from a given directory are read in (matching them by filenames).
#' Graphs and CSV files containing the miRNA concentrations (final results) are
#' generated and stored in the given output folder (one subfolder for each batch).
#' 
#' For each sample the following files will be generated:
#'  - _concentrations.csv containing all miRNAs with reads > 0 and the calculated concentrations in molecules/uL
#'  - _spikein_qc.pdf a graph showing the spike in model performance and distribution of miRs in relation to the spikeins
#'
#' For one batch (=set of multiple samples) in addition those files will be generated to hold
#' data of all batch samples:
#'  - concentrations.csv
#'  - spikein_stats.csv containing model parameters and important QC information for all samples of the batch
#'
#' Tested/developed on:
#'   R version 4.1.0 (2021-05-18)
#'   Platform: x86_64-pc-linux-gnu (64-bit)
#'   Running under: Ubuntu 20.04.2 LTS

# Setup START --------------------
# The following variables must be set before running the script

# define the directory where -metadata.csv and -miRNAsTableRC.csv files are located
input_data_path <- "testdata"

# set the path to the folder where generated output files should be stored
output_path <- "output"

# Setup END -----------------------

# Load needed libraries
library(tidyverse) # tested with v1.3.1
library(gridExtra) # tested with v2.3

# check if input_data_path exists and contains at least 1 set of files (metadata and RC)
if (!dir.exists(input_data_path)) {
  stop("directory set in varibale input_data_path can not be read")
} else {
  input_data_files_metadata <- list.files(path = input_data_path, pattern = "*-metadata.csv", full.names = TRUE)
  input_data_files_rc <- list.files(path = input_data_path, pattern = "*-miRNAsTableRC.csv", full.names = TRUE)

  if (input_data_files_metadata %>% length() != input_data_files_rc %>% length() || input_data_files_metadata %>% length() < 1) {
    stop("No or invalid metadata and/or rc files found in input data path")
  }
}

# check if output_path exists and create it if not
if (!dir.exists(output_path)) {
  dir.create(output_path, showWarnings = FALSE)
}

# read all *-metadata.csv files and *-miRNAsTableRC.csv files into a tibble

for (metadata_file in input_data_files_metadata) {
  # create a sub directory in the output path
  batch_output_path <- paste(output_path, metadata_file %>% str_replace("-metadata.csv", ""), sep = "/")
  dir.create(batch_output_path, recursive = TRUE, showWarnings = FALSE)

  # get the filename of the corresponding RC file
  rc_file <- metadata_file %>% str_replace("-metadata.csv", "-miRNAsTableRC.csv")

  # check if both files exist
  if (!file.exists(metadata_file)) {
    stop(paste("medata.csv file missing: ", metadata_file))
  } else if (!file.exists(rc_file)) {
    stop(paste("miRNAsTableRC.csv file missing: ", rc_file))
  }

  # read metadata file which contains one row for each sample
  metadata_file_content <- read_csv(metadata_file, show_col_types = FALSE)

  # read rc file which contains one row per miRNA and one column per sample
  # make data longer so that we have 3 columns (miRNA, customerID (which is the sample_id) and the read count)
  rc_file_content <- read_csv(rc_file, show_col_types = FALSE) %>%
    pivot_longer(!customerID, names_to = "sample_id", values_to = "rc") %>%
    select(miRNA = customerID, customerID = sample_id, rc)

  # merge both and tidy data
  merged_batch_data <- rc_file_content %>%
    left_join(metadata_file_content, by = "customerID") %>%
    select(miRNA, customerID, rc, contains("Spike-In_")) %>%
    arrange(customerID, miRNA)

  # prepare a tibble to hold spike in results for all samples of this batch
  spikeins_stats <- tibble(
    sample_id = character(),
    spikeins_detected = numeric(),
    mirnas_in_range = numeric(),
    spikein_lower_limit = numeric(),
    spikein_upper_limit = numeric(),
    intercept = numeric(),
    slope = numeric(),
    rsq = numeric(),
    qc = character()
  )

  batch_sample_data <- tibble()

  # now loop through every sample and estimate concentrations based on the spike ins
  for (sample_id in merged_batch_data$customerID %>% unique()) {
    # filter data for this sample
    sample_data <- merged_batch_data %>%
      filter(customerID == sample_id)

    # calculate model based on spike in concentrations and read counts ------------------

    # extract spike in sequences
    sample_spikeins <- sample_data %>%
      filter(grepl("#SpikeIn", miRNA)) %>%
      select(spikein = miRNA, everything()) %>%
      mutate(spikein = str_replace(spikein, "#SpikeIn.", ""))

    # filter spikeins from the sample_data tibble
    sample_data <- sample_data %>%
      filter(!grepl("#SpikeIn", miRNA)) %>%
      select(-contains("Spike-In_"), -customerID) %>%
      filter(rc > 0) # and remove miRNAs that have 0 RC


    # extract spike in concentrations for this sample
    spikein_concentrations <- c(
      "C" = sample_spikeins[[1, "Spike-In_C_amol"]], # take the value of the first row as they are all the same for one sample
      "E" = sample_spikeins[[1, "Spike-In_E_amol"]],
      "H" = sample_spikeins[[1, "Spike-In_H_amol"]],
      "I" = sample_spikeins[[1, "Spike-In_I_amol"]],
      "K" = sample_spikeins[[1, "Spike-In_K_amol"]],
      "M" = sample_spikeins[[1, "Spike-In_M_amol"]],
      "N" = sample_spikeins[[1, "Spike-In_N_amol"]]
    )

    # remove the concentration columns from the tibble as they are no longer needed
    sample_spikeins <- sample_spikeins %>%
      select(-contains("Spike-In_"), -customerID) %>%
      arrange(spikein)

    # define some variables used for the calculations
    amol <- 602214 # number of molecules in 1 amol
    finalVolume <- 9.5 # final volume of sample + spikeInsVolume uL of the spikeins
    spikeInsVolume <- 1

    # calculate finale concentrations of spike ins in the sample
    spikein_concentrations <- spikein_concentrations * spikeInsVolume * amol / finalVolume
    sample_spikeins <- sample_spikeins %>% left_join(
      tibble(
        spikein = names(spikein_concentrations),
        concentration = spikein_concentrations
      ),
      by = "spikein"
    )

    # check if there are spikeins missing (RC of 0)
    spikein_missing <- sample_spikeins %>%
      filter(rc == 0)
    # remove them from the spike ins tibble
    sample_spikeins <- sample_spikeins %>%
      filter(rc > 0)

    if (spikein_missing %>% nrow() > 1) {
      warning(paste("Missing RC for 2 or more spike ins in sample", sample_id))
      spikeins_stats <- spikeins_stats %>% add_row(
        sample_id = sample_id,
        spikeins_detected = sample_spikeins %>% filter(rc < 0) %>% nrow(),
        mirnas_in_range = NA,
        spikein_lower_limit = NA,
        spikein_upper_limit = NA,
        intercept = NA,
        slope = NA,
        rsq = NA,
        qc = "failed (2 or more spike ins missing)"
      )
    } else {
      # calculate slope based from the spike ins, intercept is forced to 0
      fm <- lm(concentration ~ 0 + rc, data = sample_spikeins)

      intercept <- 0 # as per definition in lm
      slope <- as.matrix(coef(fm))[1]
      rsq <- summary(fm)$r.squared

      # calculate prediction intervals
      pred_int <- suppressWarnings(predict(fm, interval = "prediction")) %>% as_tibble()

      # add prediction intervals back to the spike ins tibble
      sample_spikeins <- sample_spikeins %>%
        bind_cols(pred_int)

      # spike ins range
      spikein_upper_limit <- max(sample_spikeins$rc)
      spikein_lower_limit <- min(sample_spikeins$rc)

      # check for miRNAs in range
      mirnas_in_range <- sample_data %>%
        filter(rc >= spikein_lower_limit & rc <= spikein_upper_limit)

      # calculate the ratio of miRNAs in range to all detected miRNAs
      mirnas_in_range_ratio <- (mirnas_in_range %>% nrow()) / (sample_data %>% nrow())

      # some basic QC checks
      qc_result <- c()
      if (sample_spikeins %>% nrow() != length(spikein_concentrations)) {
        qc_result <- c(qc_result, "warning (not all spikeins detected)")
      }
      if (is.na(slope)) {
        qc_result <- c(qc_result, "FAILED (model failed)", )
      }
      if (mirnas_in_range_ratio < 0.5) {
        qc_result <- c(qc_result, "warning (less than 50% of miRNAs in spikeins range)")
      }
      if (rsq < 0.95) {
        qc_result <- c(qc_result, "FAILED (R squared less than 0.95)")
      }

      if (qc_result %>% length() == 0) {
        qc_result <- "OK"
      } else {
        qc_result <- paste(qc_result, collapse = " & ")
      }

      # save results of spikeins for this sample back to the collection tibble
      spikeins_stats <- spikeins_stats %>% add_row(
        sample_id = sample_id,
        spikeins_detected = sample_spikeins %>% nrow(),
        mirnas_in_range = mirnas_in_range %>% nrow(),
        spikein_lower_limit = spikein_lower_limit,
        spikein_upper_limit = spikein_upper_limit,
        intercept = intercept,
        slope = slope,
        rsq = rsq,
        qc = qc_result
      )

      # now calculate the miRNA concentrations based on the spikeins model
      mirnas_predicted <- predict(fm, data.frame(rc = sample_data$rc), interval = "p") %>%
        as_tibble() %>%
        mutate(interval = round(fit - lwr, 2)) %>%
        mutate(concentration = round(fit, 2)) %>%
        select(concentration, interval)

      # and add them back to the sample_data tibble
      sample_data <- sample_data %>%
        bind_cols(mirnas_predicted)

      # add range information for each miRNA to know if the rc of the miR was in
      # range of the spikeins
      sample_data <- sample_data %>%
        mutate(range = ifelse(rc < spikein_lower_limit, "too low", "in range")) %>%
        mutate(range = ifelse(rc > spikein_upper_limit, "too high", range))

      # save mirna concentrations to csv
      sample_data %>%
        select(-interval) %>%
        write_csv2(paste0(batch_output_path, "/", sample_id, "_concentrations.csv"))

      # add data to collection tibble for the whole batch
      batch_sample_data <- batch_sample_data %>%
        mutate(customerID = sample_id) %>%
        select(customerID, everything()) %>%
        bind_rows(sample_data)

      sample_spikein_data <- sample_data %>%
        mutate(isSpikein = "no") %>%
        bind_rows(sample_spikeins %>%
          mutate(isSpikein = "yes") %>%
          select(miRNA = spikein, everything()))


      # generate plots and save them with the sample filename

      p1 <- ggplot(sample_spikeins, aes(rc, concentration)) +
        geom_point() +
        geom_line(aes(x = rc, y = fit), color = "blue", linetype = "dashed", size = 0.2) +
        scale_x_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)) +
        scale_y_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)) +
        annotation_logticks(size = 0.1) +
        labs(
          title = "Spike-in calibrator fit",
          subtitle = paste0("R-squared = ", round(rsq, 4)),
          x = "Reads",
          y = "Copies [molecules/uL]"
        )

      p2 <- ggplot(sample_spikein_data, aes(y = rc, x = paste0(sample_id), color = range)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = spikein_lower_limit, alpha = .1, fill = "black") +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = spikein_upper_limit, ymax = Inf, alpha = .1, fill = "black") +
        geom_point(data = function(x) {
          filter(x, isSpikein != "yes")
        }, aes(x = paste0(sample_id), y = rc), alpha = 0.5, position = position_jitter(width = .3)) +
        geom_hline(data = function(x) {
          filter(x, isSpikein == "yes")
        }, aes(yintercept = rc), alpha = 0.5) +
        geom_text(data = function(x) {
          filter(x, isSpikein == "yes")
        }, aes(y = rc, x = Inf, label = paste0(rc, " (", miRNA, ")")), alpha = 0.9, color = "black", vjust = 0, hjust = 1, size = 3) +
        scale_y_log10() +
        expand_limits(y = 0.01) +
        labs(
          title = "Reads distribution",
          subtitle = paste0(mirnas_in_range_ratio %>% round(4) %>% "*"(100), "% miRNAs in range"),
          x = "Sample",
          y = "Reads",
          color = "Quantification"
        ) +
        theme(
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 0),
          legend.text = element_text(size = 10)
        )

      # arrange and save both plots
      g <- arrangeGrob(p1, p2, layout_matrix = rbind(c(1, 1, 2)), nrow = 1)
      ggsave(file = paste0(batch_output_path, "/", sample_id, "_spikein_qc.pdf"), g, scale = 1.5)
    }
  }

  # save spikein_stats into a csv
  spikeins_stats %>%
    write_csv2(paste0(batch_output_path, "/spikein_stats.csv"))

  # save concentrations of all batch samples into a csv
  batch_sample_data %>%
    write_csv2(paste0(batch_output_path, "/concentrations.csv"))
}
