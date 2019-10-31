library(tidyverse)

load("data/snr/2016125_ICEPRO_CAST_002_160504_163523_URC.tsv.RData")

df <- cops$LuZ %>%
  as_tibble() %>%
  add_column(depth = cops$Depth, .before = 1) %>%
  pivot_longer(-depth, names_to = "wavelength", values_to = "luz") %>%
  mutate(wavelength = parse_number(wavelength)) %>%
  filter(wavelength %in% c(443, 510, 683))

luz_fitted <- cops$LuZ.fitted %>%
  as_tibble() %>%
  add_column(depth = cops$depth.fitted, .before = 1) %>%
  pivot_longer(-depth, names_to = "wavelength", values_to = "luz") %>%
  mutate(wavelength = parse_number(wavelength)) %>%
  filter(wavelength %in% c(443, 510, 683))

write_csv(df, "data/snr_luz_raw.csv")
write_csv(luz_fitted, "data/snr_luz_fitted.csv")

get_sd <- function(df) {
  z.discr <- c(seq(000.00, 001.00, by = 0.01),
               seq(001.02, 002.00, by = 0.02),
               seq(002.05, 005.00, by = 0.05),
               seq(005.10, 010.00, by = 0.10),
               seq(010.20, 020.00, by = 0.20),
               seq(020.50, 050.00, by = 0.50),
               seq(051.00, 100.00, by = 1.00),
               seq(102.00, 200.00, by = 2.00))
  sdev <- array(NA, length(z.discr))
  for (i.z in 1:(length(z.discr)-1)) {
    idx <- which((df$depth >= z.discr[i.z]) & (df$depth < z.discr[i.z + 1]))
    sdev[i.z] <- sd(df$value[idx])
  }
  res <- data.frame(z.discr, sdev)
}

# work on the HDD - set the WD and get the files list
files.list <- list.files("data/snr/", pattern = "*.RData")

###***********************
### main loop on all files
###***********************
for (i.file in 1:length(files.list)) {
  myfile <- files.list[i.file]
  # myfile <- files.list[2]
  load(paste("data/snr/", myfile, sep = "")) # return a variable (list) called "cops" with everything in it

  # get date and time from filename
  filename.split <- stringr::str_split(myfile, "_")
  i.date <- filename.split[[1]][5]
  i.time <- filename.split[[1]][6]

  i.year  <- substr(i.date, 0, 2)
  i.month <- substr(i.date, 3, 4)
  i.day   <- substr(i.date, 5, 6)
  i.hour  <- substr(i.time, 0, 2)
  i.min   <- substr(i.time, 3, 4)
  i.sec   <- substr(i.time, 5, 6)

  my.date.txt <- paste(i.year, "/", i.month, "/", i.day, "-", i.hour, ":", i.min, ":", i.sec, sep = "")
  my.date <- strptime(my.date.txt, "%y/%m/%d-%H:%M:%S")

  # get the wavelengths, all first, then the red ones
  wl <- cops$LuZ.waves
  wl.keep <- as.character(wl) #c("589", "625", "665", "683", "694", "710")

  # define the depth offset between LuZ depth and ref. depth (historically not applied to the measured data!)
  offset <- cops$delta.capteur.optics["LuZ"]

  # define the measured LuZ data frame, especially tidy it!
  LuZ.mes <- as.data.frame(cops$LuZ)[wl.keep] %>%
    mutate(depth = (cops$Depth) + offset) %>%
    gather(lambda, value, - depth) %>%
    mutate(lambda = parse_number(lambda))

  # group by depth and lambda, then compute sdev of the grouped data
  # on the fitted data depth grid
  LuZ.sd <- LuZ.mes %>%
    group_by(lambda) %>%
    nest() %>%
    mutate(sdev = map(data, get_sd)) %>%
    select(- data) %>%
    unnest() %>%
    rename(depth = z.discr) %>%
    select(depth, lambda, sdev) %>%
    mutate(depth = signif(depth,3))

  # get the fitted data in the data frame also
  # note: I have to add the "signif" line otherwise
  # I get issue later whn joining, for some reason
  # it is tempering with the numbers representation!?!?!?
  LuZ.fit <- as.data.frame(cops$LuZ.fitted)[wl.keep] %>%
    mutate(depth = (cops$depth.fitted)) %>%
    gather(lambda, value, - depth) %>%
    mutate(lambda = parse_number(lambda)) %>%
    mutate(depth = signif(depth,3))

  # join the fitted values df and the sdev df to compute the relative sdev
  LuZ.sd.norm <- LuZ.sd %>%
    full_join(LuZ.fit, by = c("lambda", "depth")) %>%
    mutate(rel_sdev = 100.0 * sdev / value) %>%
    select(-sdev, -value) %>%
    mutate(datetime = as.POSIXct(my.date))
  #drop_na()

  df <- LuZ.sd.norm %>%
    mutate(date = as.Date(datetime))

  export_file <- paste0("data/snr_cops_", unique(df$date), ".csv")
  write_csv(df, export_file)

}
