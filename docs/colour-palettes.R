# colour palettes and scales

pal_coral_proxies <- c(#"Sr/Ca & d18O" = "#F21A00", #  "#33a02c",
  "SrCa & d18O" =  "#a6d854", 
  "SrCa" = "#fc8d62",
  #"d18O" = "#67a9cf",
  "d18O" = "#8da0cb",
  "ERSSTv5" = "#33a02c",
  "Sr/Ca" = "#ef8a62"
  )


lbl_coral_proxies <- list("SrCa & d18O" = expression("Sr/Ca & delta^18*O"),
                          "SrCa" = expression("Sr/Ca"),
                          "d18O" = expression(delta^18*O))


lbl_SNR_components <- list(S_clim=expression(S[clim]),
                 S_proxy=expression(S[proxy]),
                 S_noise=expression(S[noise]),
                 "S_resids (noise)" = expression(S[resids]*" (noise)"),    
                 "S_fttd (climate)" = expression(S[fttd]*" (climate)"), 
                 "S_instrumental" = expression(S[instrumental]),
                 "S_ERSSTv5" = expression(S[ERSSTv5]),
                 "S_CMIP5" = expression(S[CMIP5]))



lbl_SNR_coral_proxies <- list("SrCa & d18O" = expression("Sr/Ca & delta^18*O"),
                              "SrCa" = expression("Sr/Ca"),
                              "d18O" = expression(delta^18*O),
                              S_clim=expression("Proxy corrected ("*S[clim]*")"),
                              S_proxy=expression("Proxy uncorrected ("*S[proxy]*")"),
                              S_noise=expression("Proxy noise"),
                              S_resids=expression("Calibration residuals"),
                              "S_resids (noise)" = expression(S[resids]*" (noise)"),    
                              "S_fttd (climate)" = expression(S[fttd]*" (climate)"), 
                              "S_instrumental" = expression(S[instrumental]),
                              "S_ERSSTv5" = expression("Instrumental (ERSSTv5)"),
                              "S_CMIP5" = expression("CMIP5"))


pal_coral_taxa <- c("Porites lobata", "Porites lutea", "Porites sp.",
                    "Diploastrea heliopora", "Orbicella faveolata",
                    "Pseudodiploria strigosa", "Siderastrea siderea")


pal_SNR_components <- c(
  "S_proxy" = "#7570b3", 
  "S_noise" = "#d95f02",
  "S_resids (noise)" = "#d95f02",
  "S_resids" = "#33a02c",
  "S_clim" = "#1b9e77", 
  "S_fttd (climate)" = "#1b9e77", 
  "S_instrumental" = "#e7298a",
  "S_ERSSTv5" = "#e7298a",
  "S_CMIP5" = "Gold"
  )


lty_SNR_components <- c(
  "S_clim" = 1, 
  "S_proxy" = 1, 
  "S_noise" = 1,
  "S_instrumental" = 2,
  "S_ERSSTv5" = 2,
  "S_CMIP5" = 3
)


Spec_SrCa_d18O_labeller <- ggplot2::as_labeller(
  c(SrCa="Sr/Ca", d18O = "delta^18*O",
    S_clim="S[clim]", S_proxy="S[proxy]",
    S_noise="S[noise]"),
  default = label_parsed)


