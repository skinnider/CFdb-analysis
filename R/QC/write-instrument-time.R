# Write the instrument time associated with each file in the combined analysis.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser()
parser$add_argument('--base_dir', type = 'character', 
                    default = "~/projects/rrg-ljfoster-ab/skinnim/CFdb")
parser$add_argument('--rawtools_dir', type = 'character', 
                    default = "data/QC/RawTools")
args = parser$parse_args()

library(tidyverse)
library(magrittr)

# first, read list of all files
files = read.csv("~/git/CFdb-searches/data/files.csv") %>%
  # expand full filepath
  mutate(filepath = file.path(args$base_dir, Accession, 
                              gsub("^.*\\\\", "", File)))

# fix filepaths
files %<>%
  mutate(Replicate = fct_recode(Replicate,
                                'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis') %>% 
           as.character(),
         filepath = ifelse(Author == 'SaeLee',
                           file.path(dirname(filepath), 
                                     Replicate,
                                     basename(filepath)),
                           filepath))

# iterate over filepaths
filepaths = files$filepath

# create empty column
files$Time = NA
for (idx in seq_along(filepaths)) {
  filepath = filepaths[idx]
  message("[", idx, "/", length(filepaths), "] working on file: ", filepath)
  
  # ignore files that don't exist
  if (!file.exists(filepath)) {
    warning("  ignoring file ", filepath)  
    next
  }
  
  # get instrument time
  if (endsWith(tolower(filepath), '.raw')) {
    # read instrument time from RawTools output
    rawtools_file = file.path(args$rawtools_dir, paste0(basename(filepath), 
                                                        '_Metrics.txt.gz'))
    if (!file.exists(rawtools_file)) {
      warning("file does not exist: ", rawtools_file)
      next
    }
    files$Time[idx] = readLines(rawtools_file) %>%
      extract(grepl("TotalAnalysisTime", .)) %>%
      gsub("^.*:", "", .) %>%
      as.numeric()
  } else if (endsWith(tolower(filepath), '.mzxml')) {
    # parse mZXML plain text directly
    con = file(filepath, open = 'r')
    start_time = end_time = NULL
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      if (grepl("startTime", line)) {
        start_time = line %>%
          gsub("^.*startTime=\"", "", .) %>%
          gsub("\".*$", "", .) %>%
          gsub("^PT|S$", "", .)
      }
      if (grepl("endTime", line)) {
        end_time = line %>%
          gsub("^.*endTime=\"", "", .) %>%
          gsub("\".*$", "", .) %>%
          gsub("^PT|S$", "", .)
      }
      if (!is.null(start_time) && !is.null(end_time)) {
        files$Time[idx] = (as.numeric(end_time) - as.numeric(start_time)) / 60
        break
      }
    }
    close(con)
  } else if (endsWith(tolower(filepath), '.wiff') |
             endsWith(tolower(filepath), '.d')) {
    warning('not yet implemented: Bruker/AB Sciex (', filepath, ')')
  }
}

# print missing values
filter(files, is.na(Time)) %>% 
  dplyr::count(Author, Year, Accession)

# fill in a few missing values by manual review
#' McBride2017, wiff files: 90 min gradient
#' Aryal2017, wiff files: 120 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD006694'] = 90
files$Time[is.na(files$Time) & files$Accession == 'MSV000081206'] = 120
#' v2: 
#' Spaniol2022 PXD023443 = 21 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD023443'] = 21
#' Lee2021 PXD021451 = 125 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD021451'] = 125
#' Lee2021 PXD021454 = 125 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD021454'] = 125
#' Lee2021 PXD021452 = 90 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD021452'] = 125
#' Kerr2020 PXD013809 = 90 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD013809'] = 90
#' Skinnider2021 PXD022309 = 150 min gradient
files$Time[is.na(files$Time) & files$Accession == 'PXD022309'] = 150

# are any values still missing?
str(filter(files, is.na(Time)))

# print summary
message("total instrument time: ", sum(files$Time, na.rm = T), " min")

# write output table
write.csv(files, 'data/QC/instrument-time.csv', row.names = F)

# remove Kerr2020 rep1
files %>% 
  filter(!(Author == 'Kerr' & grepl('rep1', Replicate))) %$% 
  sum(Time, na.rm = TRUE) %>% 
  message("total instrument time: ", ., " min ")
