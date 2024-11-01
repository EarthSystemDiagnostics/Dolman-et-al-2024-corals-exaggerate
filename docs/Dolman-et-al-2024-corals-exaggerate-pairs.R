## ----setup, include=FALSE----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE,
                      cache = TRUE,
                      fig.width = 8, fig.height = 5,
                      dev.args = list(png = list(type = "cairo")))


## ----load_packages-----------------------------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(progress)
library(ccaPP)


## ----load_github_packages----------------------------------------------------------------------------------------------------
if(!require(PaleoSpec)){
    remotes::install_github("EarthSystemDiagnostics/paleospec")
    library(PaleoSpec)
}


## ----load_additional_functions-----------------------------------------------------------------------------------------------
source("colour-palettes.R")
source("functions.R")


## ----------------------------------------------------------------------------------------------------------------------------
# equalise the variance of all records in a cluster
eq_var <- FALSE

# use SSA gap filled version of data
use_ssa <- TRUE

# maximum pairwise distance allowed
max_dist <- 100

# minimum correlation in OISST allowed
min_rho <- 0.99

# which SST product to use when using only one
SST_prod_to_use <- "SST_OISST2_dis"

# use cached bootstrap?
use_cache = TRUE

n_boot <- 100


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_sub0 <- readRDS("../data/coral_pairs_sub0_mon.rds")

coral_pairs_sub0 <- coral_pairs_sub0 %>% 
  mutate(paleoData_variableName = factor(paleoData_variableName,
                                         ordered = TRUE,
                                         levels = c("SrCa", "d18O")))
coral_pairs_sub0 <- coral_pairs_sub0 %>% 
  mutate(genus = gsub("([A-Za-z]+).*", "\\1", paleoData_archiveSpecies)) %>% 
  mutate(overlap_years = round(last_joint_year.dec - first_joint_year.dec, 2))


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_sub <- coral_pairs_sub0 

if (use_ssa == TRUE){
  coral_pairs_sub$proxy <- coral_pairs_sub$ssa_values
} else {
  coral_pairs_sub$proxy <- coral_pairs_sub$paleoData_values
}


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_sub <- coral_pairs_sub %>%
  group_by(paleoData_TSid) %>%
  mutate(longest_instance = overlap_years == max(overlap_years))

coral_pairs_once <- coral_pairs_sub %>%
  filter(longest_instance == TRUE)


coral_pairs_meta <- coral_pairs_sub %>%
  ungroup() %>%
  select(
    pair, dist_km, overlap_years, member,
    paleoData_TSid, paleoData_variableName,
    hasResolution_nominal, paleoData_archiveSpecies, genus,
    geo_description, geo_longitude, geo_latitude
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = c(member),
    values_from = c(
      paleoData_TSid, paleoData_variableName, hasResolution_nominal,
      paleoData_archiveSpecies, genus, geo_description,
      geo_longitude, geo_latitude
    )
  )

stopifnot(
  length((coral_pairs_meta$pair)) == length(unique(coral_pairs_meta$pair))
  )


## ----------------------------------------------------------------------------------------------------------------------------
OISST_ann <- readRDS("../data/OISST2_ch2k_remap_int.rds") %>% 
  group_by(paleoData_ch2kCoreCode, year) %>% 
  summarise(SST_OISST2_dis = mean(SST_OISST2_dis, na.rm = TRUE))

OISST_ann_pairs <- coral_pairs_sub %>% 
  select(pair, dist_km, overlap_years, resolution, paleoData_ch2kCoreCode,
         paleoData_variableName, paleoData_TSid, member) %>% 
  distinct() %>% 
  left_join(., OISST_ann, relationship = "many-to-many") %>% 
  ungroup() %>% 
  select(pair, paleoData_variableName, dist_km, overlap_years, resolution,
         member, year, SST_OISST2_dis) %>% 
  pivot_wider(names_from = member, values_from = SST_OISST2_dis) %>% 
  group_by(pair) %>% 
   # mutate(a2 = pracma::detrend(a),
   #        b2 = pracma::detrend(b)
   #        ) %>% 
     mutate(a2 = residuals(lm(a~year)),
          b2 = residuals(lm(b~year))
          )

OISST_pair_rho <- OISST_ann_pairs %>% 
  group_by(pair, paleoData_variableName, dist_km, overlap_years, resolution) %>% 
  summarise(rho = as.numeric(cor(a2, b2)))

OISST_pair_rho %>% 
  filter(#dist_km <= 100,
         rho >= min_rho)


## ----------------------------------------------------------------------------------------------------------------------------
OISST_pair_rho %>% 
 filter(resolution == "mon") %>%
  ggplot(aes(x = dist_km, y = rho)) +
  geom_hline(yintercept = c(min_rho), lty = 2, colour = "black") +
  geom_point(alpha = 0.5) +
  expand_limits(y = c(0,1)) +
  #geom_smooth(method = "lm") +
  labs(y = "Pearson correlation coefficient") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #scale_x_sqrt(labels = scales::comma) +
  scale_y_continuous(breaks = c(seq(0, 0.99, length.out = 4))) +
  labs(x = "Pairwise distance [km]")


## ----------------------------------------------------------------------------------------------------------------------------
pairs_to_use <- OISST_pair_rho %>% 
  filter(resolution %in% c("mon"),
         rho >= min_rho)


## ----------------------------------------------------------------------------------------------------------------------------
f_breaks <- c(1/100, 1/10, 1/2, 1/1.1, 1.1/1)
f_names <- as.character(expression(1/100, 1/10, 1/2, 1/1.1, 1.1/1))

fbands <- tibble(
  f_upr = tail(f_breaks, -1),
  f_lwr = head(f_breaks, -1),
  f_band = paste0("[", head(f_names, -1), ", ", tail(f_names, -1), "]")
  ) %>% 
  mutate(f_band = factor(f_band, ordered = TRUE, levels = f_band))

fbands$f_band_name <- factor(c("Decadal to\nCentennial",
                               "Interannual\nto Decadal",
                               "Interannual",
                               "Annual"),
                             ordered = TRUE,
                             levels = c("Decadal to\nCentennial",
                                        "Interannual\nto Decadal", 
                                        "Interannual",
                                        "Annual"))



## ----------------------------------------------------------------------------------------------------------------------------
genus_count <- coral_pairs_sub %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  select(paleoData_variableName, paleoData_TSid,
         paleoData_archiveSpecies, genus) %>% 
  distinct() %>% 
  group_by(genus) %>% 
  summarise(n= n()) %>% 
  arrange(desc(n))

knitr::kable(genus_count)


## ----inset_figures_2---------------------------------------------------------------------------------------------------------
sub_mon <- coral_pairs_sub %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  filter(paleoData_TSid == "FL17DTO01_SrCa") 

sub_ann <- coral_pairs_sub %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  filter(paleoData_TSid == "FL17DTO01_SrCa") %>% 
  group_by(year, paleoData_TSid) %>% 
  summarise_if(is.numeric, mean) 

lm_ERSST <- sub_mon %>% 
  lm(proxy~offset(-0.06*SST_ERSST5_dis), data = .)

cfs_ERSST <- coef(lm_ERSST)


lm_OISST <- sub_mon %>% 
  lm(proxy~offset(-0.06*SST_OISST2_dis), data = .)

cfs_OISST <- coef(lm_OISST)


## ----inset_figures_big-------------------------------------------------------------------------------------------------------
ins_mon_long <- sub_mon %>% 
  filter(year > 1979, year < 2001,
         complete.cases(SST_OISST2_dis)) %>% 
  ggplot(aes(x = year.dec, y = proxy)) +
  geom_line(aes(colour = "SrCa")) +
  geom_line(aes(y = SST_OISST2_dis * -0.06 + cfs_OISST[1], colour = "OISSTv2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)
        ) +
  scale_color_manual(values = c(pal_coral_proxies, OISSTv2 = "#66a61e")) +
  labs(x = "Year CE", y = "Monthly Sr/Ca [mmol/mol]", colour = "", linetype = "")


ins_ann_long <- sub_ann %>% 
  filter(year > 1979, year < 2001,
         complete.cases(SST_OISST2_dis)) %>% 
  ggplot(aes(x = year.dec, y = proxy)) +
  geom_line(aes(colour = "SrCa")) +
  geom_line(aes(y = SST_OISST2_dis * -0.06 + cfs_OISST[1],
                colour = "OISSTv2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)
        ) +
  scale_color_manual(values = c(pal_coral_proxies, OISSTv2 = "#66a61e")) +
  labs(x = "Year CE", y = "Annual mean Sr/Ca [mmol/mol]",
       colour = "", linetype = "")


patchwork::wrap_plots(ins_mon_long, ins_ann_long, ncol = 1)


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_once_SST_long <- coral_pairs_once %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  filter(resolution == "mon") %>% 
  #filter(complete.cases(SST * proxy)) %>% 
  group_by(paleoData_variableName, paleoData_ch2kCoreCode, paleoData_TSid,
           genus, paleoData_archiveSpecies, year
           ) %>% 
  select(-pair) %>% 
  distinct() %>% 
  pivot_longer(cols = starts_with("SST_"),
               names_to = "SST_product",
               values_to = "SST") 

month_SD <- coral_pairs_once_SST_long %>% 
  group_by(SST_product, paleoData_variableName, paleoData_ch2kCoreCode,
           paleoData_TSid, genus, paleoData_archiveSpecies, year
           ) %>% 
  summarise(n = sum(is.nan(SST * proxy)==FALSE),
            sd_SST = sd(SST, na.rm = TRUE),
            sd_SST = sd(SST, na.rm = TRUE),
            sd_proxy = sd(proxy, na.rm = TRUE)) %>% 
  filter(n >= 12) %>% 
  select(-year) %>% 
  summarise_all(.funs = list(mean=mean, sd=sd), na.rm = TRUE) %>% 
  mutate(across(ends_with("_sd"), ~./sqrt(n_mean)))


SESD <- function(S, n){
  S / sqrt(2*n - 2)
}

coral_pairs_sub_ann <- coral_pairs_once_SST_long %>% 
  filter(complete.cases(proxy, SST)) %>% 
  group_by(SST_product, paleoData_variableName,
           paleoData_ch2kCoreCode, paleoData_TSid,
           paleoData_archiveSpecies, genus, year) %>% 
  #select(-pair) %>% 
  distinct() %>% 
  mutate(n = sum(is.nan(SST * proxy)==FALSE)) %>% 
  filter(n >= 12) %>% 
  summarise_if(is.numeric, mean)


ann_SD <- coral_pairs_sub_ann %>% 
  select(-year) %>% 
  summarise(sd_SST = sd(SST, na.rm = TRUE),
            sd_proxy = sd(proxy, na.rm = TRUE),
            n = n())%>% 
  mutate(se_sd_SST = SESD(sd_SST, n),
         se_sd_proxy = SESD(sd_proxy, n))


dec_SD <- coral_pairs_sub_ann %>% 
  mutate(decade = round(year, -1)) %>% 
  select(-year) %>% 
  group_by(SST_product, paleoData_variableName, paleoData_ch2kCoreCode,
           paleoData_TSid, paleoData_archiveSpecies, genus, decade) %>% 
  select(-decade) %>% 
  group_by(SST_product, paleoData_variableName, paleoData_ch2kCoreCode,
           paleoData_TSid, paleoData_archiveSpecies, genus) %>% 
  summarise(sd_SST = sd(SST, na.rm = TRUE),
            sd_proxy = sd(proxy, na.rm = TRUE),
            n = n())%>% 
  mutate(se_sd_SST = SESD(sd_SST, n),
         se_sd_proxy = SESD(sd_proxy, n))
  


## ----fig.width=5, fig.height = 3.5-------------------------------------------------------------------------------------------
#"FL17DTO01_SrCa"
highlt_month <- month_SD %>% 
  filter(paleoData_ch2kCoreCode == "FL17DTO01",
         paleoData_variableName == "SrCa")%>% 
  filter(SST_product %in% SST_prod_to_use)

fig_sd_mon <- month_SD %>% 
  filter(SST_product %in% SST_prod_to_use,
         complete.cases(sd_SST_mean)) %>% 
  filter(paleoData_TSid %in% c("AL16YUC01_SrCa", "AL16PUR02_SrCa") == FALSE) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  ggplot(aes(x = sd_SST_mean, y = sd_proxy_mean
             )) +
  geom_linerange(aes(ymin = sd_proxy_mean - sd_proxy_sd, 
                      ymax = sd_proxy_mean + sd_proxy_sd), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST_mean - sd_SST_sd,
                     xmax = sd_SST_mean + sd_SST_sd), alpha = 0.5) +
  geom_point(alpha = 0.5) +
  # highlight inset
  geom_linerange(data = highlt_month, colour = "Red",
                 aes(ymin = sd_proxy_mean - sd_proxy_sd, 
                      ymax = sd_proxy_mean + sd_proxy_sd)) + 
  geom_linerange(data = highlt_month, colour = "Red",
                 aes(xmin = sd_SST_mean - sd_SST_sd,
                     xmax = sd_SST_mean + sd_SST_sd)) +
  geom_point(data = highlt_month, colour = "Red", size = 3) +
  geom_abline(intercept = 0, slope = c(0.06), lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  expand_limits(x = c(0, 2.8), y = c(0, 0.06*2.8)) +
  theme_bw() +
  coord_fixed(ratio = 1/0.06) +
  theme(panel.grid = element_blank()) +
  labs(x = "SD annual cycle OISSTv2 [°C]", 
       y = "SD annual cycle Sr/Ca [mmol/mol]",
       colour = "Genus") 

#fig_sd_mon


## ----------------------------------------------------------------------------------------------------------------------------
highlt_ann <- ann_SD %>% 
  filter(paleoData_ch2kCoreCode == "FL17DTO01",
         paleoData_variableName == "SrCa")%>% 
  filter(SST_product %in% SST_prod_to_use)

fig_sd_ann <- ann_SD %>% 
  filter(SST_product %in% SST_prod_to_use) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  ggplot(aes(x = sd_SST, y = sd_proxy, group = paleoData_ch2kCoreCode)) +
  geom_linerange(aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST), alpha = 0.5) +
  geom_point(alpha = 0.5)+
  # highlight inset
  geom_linerange(data = highlt_ann, colour = "Red",
                 aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy)) + 
  geom_linerange(data = highlt_ann, colour = "Red",
                 aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST)) +
  geom_point(data = highlt_ann, colour = "Red", size = 3) +

  geom_abline(intercept = 0, slope = 0.06, lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  theme_bw() +
  coord_fixed(ratio = 1/0.06) +
  expand_limits(x = c(0, 1.1), y = c(0, 0.06*1.1)) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "SD annual mean OISSTv2 [°C]",
       y = "SD annual mean Sr/Ca [mmol/mol]") +
  scale_color_brewer("", type = "qual", palette = "Set2")#+
  #facet_wrap(~SST_product)

#fig_sd_ann


## ----------------------------------------------------------------------------------------------------------------------------
fig_sd_dec <- dec_SD %>%
  filter(SST_product %in% SST_prod_to_use) %>% 
  filter(paleoData_TSid %in% c("AL16YUC01_SrCa", "AL16PUR02_SrCa") == FALSE) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  ggplot(aes(x = sd_SST, y = sd_proxy, group = paleoData_ch2kCoreCode)) +
  geom_linerange(aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST), alpha = 0.5) +
  geom_point(alpha = 0.5)+
  # highlight inset
  geom_linerange(data = highlt_ann, colour = "Red",
                 aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy)) + 
  geom_linerange(data = highlt_ann, colour = "Red",
                 aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST)) +
  geom_point(data = highlt_ann, colour = "Red", size = 3) +

  geom_abline(intercept = 0, slope = 0.06, lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  scale_color_discrete("") +
  theme_bw() +
  coord_fixed(ratio = 1/0.06#, ylim = c(0, 0.08)
              ) +
  expand_limits(x = c(0, 0.15*1/0.06), y = c(0, 2 * 0.06)) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "SD decadal mean OISSTv2 [°C]",
       y = "SD decadal mean Sr/Ca [mmol/mol]") +
  scale_color_brewer("", type = "qual", palette = "Set2")#+
  #facet_wrap(~SST_product)

fig_sd_dec


## ----fig_show_issue, fig.width=8---------------------------------------------------------------------------------------------
l <- fig_sd_mon+
  labs(tag = "b")  + 
  patchwork::inset_element(
  ins_mon_long + 
    theme(legend.position = "none", 
          axis.text.y = element_blank()
          ) +
    expand_limits(x = 1980) +
     scale_x_continuous(breaks = c(1980, 1990, 2000)) +
    labs (x = "", y = ""), 
  0.55, 0.075, 0.975, 0.42, align_to = "plot") 

r <- fig_sd_ann + 
  labs(tag = "c")  + 
  patchwork::inset_element(
  ins_ann_long + 
    theme(legend.position = "top",
          legend.box.spacing = unit(0.01, "mm"),
          legend.key.spacing.y = unit(0.01, "mm"),
          legend.justification = "right",
          axis.text.y = element_blank()
          ) +
    expand_limits(x = 1980) +
    scale_x_continuous(breaks = c(1980, 1990, 2000)) +
    labs (x = "", y = ""),
  0.55, 0.075, 0.975, 0.625, align_to = "plot")+ 
  guides(colour = guide_legend(ncol = 1)) 


rr <- fig_sd_dec + 
  labs(tag = "d") + 
  guides(colour = guide_legend(ncol = 1)) 



fig_var_ill <- patchwork::wrap_plots(l, r) 
fig_var_ill


## ----fig.width=5, fig.height = 3.5-------------------------------------------------------------------------------------------
fig_sd_mon_Genus <- month_SD %>% 
  filter(SST_product %in% SST_prod_to_use,
         complete.cases(sd_SST_mean)) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  ggplot(aes(x = sd_SST_mean, y = sd_proxy_mean, colour = genus
             )) +
  geom_linerange(aes(ymin = sd_proxy_mean - sd_proxy_sd, 
                      ymax = sd_proxy_mean + sd_proxy_sd), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST_mean - sd_SST_sd,
                     xmax = sd_SST_mean + sd_SST_sd), alpha = 0.5) +
  geom_point(alpha = 1) +
  geom_abline(intercept = 0, slope = c(0.06), lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  expand_limits(y = c(0, 0.06*3), x = c(0, 3)) +
  theme_bw() +
  coord_fixed(ratio = 1/0.06) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "SD annual cycle OISSTv2 [°C]",
       y = "SD annual cycle Sr/Ca [mmol/mol]", colour = "Genus") +
  scale_color_brewer("", type = "qual", palette = "Set2") 


## ----------------------------------------------------------------------------------------------------------------------------
fig_sd_ann_Genus <- ann_SD %>% 
  filter(SST_product %in% SST_prod_to_use) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  ggplot(aes(x = sd_SST, y = sd_proxy, colour = genus)) +
  geom_linerange(aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST), alpha = 0.5) +
  geom_point(alpha = 1)+
  geom_abline(intercept = 0, slope = 0.06, lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  scale_color_discrete("") +
  theme_bw() +
  coord_fixed(ratio = 1/0.06, ylim = c(0, 1.1*0.06), xlim = c(0, 1.1)) +
  #expand_limits(x = c(0, 0.085*1/0.06), ) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(face = "italic")) +
  labs(x = "SD annual mean OISSTv2 [°C]", 
       y = "SD annual mean Sr/Ca [mmol/mol]", 
       colour = "Genus") +
  scale_color_brewer("", type = "qual", palette = "Set2")



## ----fig.width=5, fig.height = 3.5-------------------------------------------------------------------------------------------
fig_sd_mon_d18O_Genus <- month_SD %>% 
  filter(SST_product %in% SST_prod_to_use, 
         complete.cases(sd_SST_mean)) %>% 
  filter(paleoData_variableName == "d18O") %>% 
  ggplot(aes(x = sd_SST_mean, y = sd_proxy_mean, colour = genus)) +
  geom_linerange(aes(ymin = sd_proxy_mean - sd_proxy_sd, 
                      ymax = sd_proxy_mean + sd_proxy_sd), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST_mean - sd_SST_sd,
                     xmax = sd_SST_mean + sd_SST_sd), alpha = 0.5) +
  geom_point(alpha = 1) +
  geom_abline(intercept = 0, slope = c(0.22), lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  expand_limits(y = c(0, 0.22*2.5), x = c(0, 2.5)) +
  theme_bw() +
  coord_fixed(ratio = 1/0.22) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "SD annual cycle OISSTv2 [°C]",
       y = expression("SD annual cycle "*delta^18*O*" [‰]"),
       colour = "Genus") +
  scale_color_brewer("", type = "qual", palette = "Set2")


## ----------------------------------------------------------------------------------------------------------------------------
fig_sd_ann_d18O_Genus <- ann_SD %>% 
  filter(SST_product %in% SST_prod_to_use) %>% 
  filter(paleoData_variableName == "d18O") %>% 
  ggplot(aes(x = sd_SST, y = sd_proxy, colour = genus)) +
  geom_linerange(aes(ymin = sd_proxy - se_sd_proxy, 
                      ymax = sd_proxy + se_sd_proxy), alpha = 0.5) + 
  geom_linerange(aes(xmin = sd_SST - se_sd_SST,
                     xmax = sd_SST + se_sd_SST), alpha = 0.5) +
  geom_point(alpha = 1)+
  geom_abline(intercept = 0, slope = 0.22, lty = 4, #lwd = 1,
                                alpha = 1, colour = "Darkblue") +
  scale_color_discrete("") +
  theme_bw() +
  coord_fixed(ratio = 1/0.22) +
  expand_limits(x = c(0, 1.3), y = c(0, 0.22*1.3)) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "SD annual mean OISSTv2 [°C]",
       y = expression("SD annual mean "*delta^18*O*" [‰]"),
       colour = "Genus") +
  scale_color_brewer("", type = "qual", palette = "Set2")


## ----fig_show_issue_Genus, fig.width=9, fig.height=8, warning=FALSE----------------------------------------------------------
patchwork::wrap_plots(fig_sd_mon_Genus, fig_sd_ann_Genus, 
                      fig_sd_mon_d18O_Genus, fig_sd_ann_d18O_Genus, 
                      guides = "auto") +
  patchwork::plot_annotation(tag_levels = "a")




## ----------------------------------------------------------------------------------------------------------------------------
#bin_width <- 1
coral_pairs_binned <- coral_pairs_sub %>% 
  mutate(month_CE = year * 12 + month,
         bin_width = ifelse(resolution == "mon", 1, 12)) %>% 
  group_by(pair, paleoData_variableName) %>% 
  mutate(min.age = min(month_CE), max.age = max(month_CE)) %>%
  arrange(month_CE) %>%
  group_by(resolution, pair, dist_km,
           paleoData_variableName, paleoData_TSid) %>%
  reframe(
    BinTimeseries(month_CE, proxy,
                  bin.width = bin_width[1],
                  strt.time = min.age[1], end.time = max.age[1]
                  )
    ) %>%
  group_by(pair, paleoData_variableName, paleoData_TSid) %>%
  mutate(n.pts = sum(is.na(mean.value) == FALSE)) %>%
  rename(age = time) %>%
  filter(n.pts > 0) %>%
  ungroup() %>% 
  group_by(pair, paleoData_variableName) %>% 
  mutate(p_NA = sum(is.na(mean.value)) / n())


## ----------------------------------------------------------------------------------------------------------------------------
# check for number of gaps
coral_pairs_binned %>% 
  select(paleoData_TSid, resolution, p_NA) %>% 
  distinct() %>% 
  ggplot(aes(x = p_NA)) +
  geom_histogram() +
  facet_wrap(~resolution, scales = "free")


## ----------------------------------------------------------------------------------------------------------------------------
# Make each pair a matrix of time x record
coral_pairs_binned_wide <- coral_pairs_binned %>%
  filter(p_NA < 0.1) %>% 
  ungroup() %>%
  select(-n.bin, -n.pts, -mean.time) %>%
  group_by(resolution, bin.width, pair, paleoData_variableName) %>%
  do({
    m = pivot_wider(., names_from = paleoData_TSid,
                    values_from = mean.value,
                    names_prefix = "record_") %>%
     arrange(age) %>%
     select(starts_with("record_")) %>%
     as.matrix(.)
    
     colnames(m) <- NULL
  
    tibble(
      pair_size = ncol(m),
      n_gaps = sum(is.na(m)),
      p_gaps = n_gaps / length(m),
      m = list(m)
    )
  })


stopifnot(all(sapply(coral_pairs_binned_wide$m, ncol) == 2))


## ----------------------------------------------------------------------------------------------------------------------------
## Get spec of climate at locations
ERSST_pairs_binned <- coral_pairs_sub %>%
  mutate(month_CE = year * 12 + month,
         bin_width = ifelse(resolution == "mon", 1, 12)) %>% 
  group_by(paleoData_variableName, resolution, pair, 
           dist_km, paleoData_TSid) %>%
  mutate(min.age = min(month_CE), max.age = max(month_CE)) %>%
  arrange(paleoData_variableName, pair, paleoData_TSid, month_CE) %>%
  group_by(paleoData_variableName, resolution, pair, 
           dist_km, paleoData_TSid) %>%
  rename(OISSTv2 = SST_OISST2_dis, 
         ERSSTv5 = SST_ERSST5_dis) %>% 
  pivot_longer(cols = c(OISSTv2, ERSSTv5),
               names_to = "SST_product", 
               values_to = "SST") %>% 
  filter(complete.cases(SST)) %>% 
  group_by(paleoData_variableName, resolution, pair, dist_km,
           paleoData_TSid, SST_product) %>%
  filter(n() > 9) %>% 
  reframe(
    BinTimeseries(month_CE, SST,
                  bin.width = bin_width[1],
                  strt.time = min.age[1], end.time = max.age[1]
                  )
    ) %>%
  group_by(paleoData_variableName, resolution, pair,
           paleoData_TSid, SST_product) %>%
  mutate(n.pts = sum(is.na(mean.value) == FALSE)) %>%
  rename(age = time) %>%
  filter(n.pts > 0) %>%
  ungroup()


ERSST_pairs_binned_wide <- ERSST_pairs_binned %>%
  #filter(p_NA < 0.2) %>%
  ungroup() %>%
  select(-n.bin, -n.pts, -mean.time) %>%
  group_by(bin.width, resolution, paleoData_variableName, pair, dist_km,
           paleoData_TSid, SST_product) %>%
  do({
    dat <- .
    m = pivot_wider(dat, names_from = paleoData_TSid,
                    values_from = mean.value,
                    names_prefix = "record_") %>%
     arrange(age) %>%
     select(starts_with("record_")) %>%
     as.matrix(.)

     colnames(m) <- NULL

     m <- PaleoSpec:::TrimNA(m)
    tibble(
      pair_size = ncol(m),
      n_gaps = sum(is.na(m)),
      p_gaps = n_gaps / length(m),
      m = list(m)
    )
  })



## ----------------------------------------------------------------------------------------------------------------------------
ERSST_specs <- ERSST_pairs_binned_wide %>%
  filter(p_gaps < 0.2#,
         #paleoData_TSid != "RE19GBR04_SrCa"
         ) %>%
  group_by(bin.width, resolution, paleoData_variableName, pair, dist_km,
           paleoData_TSid, SST_product) %>%
   reframe(
      Spec2DF(
           SpecACF(x = m[[1]], bin.width = bin.width[1]/12, k = 5, nw = 3)
        )
   ) %>%
  mutate(spec_id = paleoData_TSid) %>%
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
ERSST_specs_reg <- ERSST_specs %>%
  group_by(spec_id, bin.width, pair, dist_km, paleoData_variableName,
           resolution, SST_product) %>%
  filter(complete.cases(spec)) %>% 
  filter(n() > 9) %>% 
  do({
    min_f <- min(.$freq)

    strt_f <-  0.01
    maxf <- 6
    rfreq = seq(strt_f, maxf, strt_f)

    sr <- as.spec(list(freq = .$freq, spec = .$spec, dof = .$dof))

    si <- PaleoSpec::SpecInterpolate(sr,
                                     freqRef = rfreq,
                                     check = FALSE)

    si <- FilterSpecLog(si)

   tibble(freq = si$freq[si$freq >= min_f], spec = si$spec[si$freq >= min_f],
          dof = si$dof[si$freq >= min_f])
  })

ERSST_specs_reg_mean <- ERSST_specs_reg %>% 
  filter(dist_km <= max_dist) %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  group_by(paleoData_variableName, SST_product, pair, freq) %>%
  summarise(n = sum(is.na(spec) == FALSE),
            se = sd(spec, na.rm = TRUE) / sqrt(n),
            spec = mean(spec, na.rm = TRUE),
            .groups = "keep") %>%
  group_by(paleoData_variableName, SST_product, freq) %>%
  summarise(n = sum(is.na(spec) == FALSE),
            se = sd(spec, na.rm = TRUE) / sqrt(n),
            spec = mean(spec, na.rm = TRUE),
            .groups = "keep") %>%
  mutate(lim.1 = spec + 2*se, lim.2 = spec - 2*se) %>%
  mutate(X97.5. = spec + 2*se, X2.5. = spec - 2*se,
         X81.4. = spec + se, X15.9. = spec - se) %>%
  mutate(spec_id = "S_instrumental") %>% 
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs <- coral_pairs_binned_wide %>% 
  filter(p_gaps < 0.2) %>% 
  #filter(resolution == "mon") %>% 
  group_by(resolution, bin.width, pair, pair_size, 
           paleoData_variableName) %>% 
  reframe(
      Spec2DF(SNRStack(x = m[[1]], 
                       prefilter = FALSE,
                       logsmooth = FALSE,
                       k = 5, nw = 3,
                       equalise_var = eq_var,
                       bin_width = bin.width[1]/12)) 
  ) %>% 
  as_spec_df() %>% 
  left_join(., coral_pairs_meta)


## ----------------------------------------------------------------------------------------------------------------------------
# interpolated each to same freq axis
coral_specs_reg <- coral_specs %>% 
 group_by(spec_id, resolution, bin.width, pair, dist_km, 
          paleoData_variableName,
          paleoData_TSid_a, paleoData_TSid_b) %>% 
  do({
    min_f <- min(.$freq)
    
    strt_f <-  0.01
    if (.$resolution[1] == "ann"){
      maxf <- 1/2
    } else {
     maxf <- 6  
    }
    
    rfreq = seq(strt_f, maxf, strt_f)
    
    sr <- as.spec(list(freq = .$freq, spec = .$spec, dof = .$dof))
    si <- PaleoSpec::SpecInterpolate(sr,
                                     freqRef = rfreq,
                                     check = FALSE)
    
    if (.$spec_id[[1]] != "SignalNoise"){
    si <- FilterSpecLog(si)  
    }
    
    tibble(freq = si$freq[si$freq >= min_f], 
           spec = si$spec[si$freq >= min_f], 
           dof = si$dof[si$freq >= min_f])
  })


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs_reg_mean <- coral_specs_reg %>%
  filter(pair %in% pairs_to_use$pair) %>% 
  group_by(spec_id, paleoData_variableName, freq) %>%
  summarise(mean_dist_km = mean(dist_km),
            n_pairs = sum(is.na(spec)==FALSE),
            se = sd(spec, na.rm = TRUE) / sqrt(n_pairs),
            mad_e = 1.2533 * mad(spec, na.rm = TRUE) / sqrt(n_pairs),
            spec_mean = mean(spec, na.rm = TRUE),
            spec_median = median(spec, na.rm = TRUE),
            .groups = "keep") %>%
  pivot_longer(starts_with("spec_m"), 
               values_to = "spec", 
               names_to = "statistic") %>%
  mutate(#spec_id = name,
         lim.2 = ifelse(statistic == "spec_mean",
                        spec - se, spec - mad_e),
         lim.1 = ifelse(statistic == "spec_mean",
                        spec + se, spec + mad_e)) %>%
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs_reg_mean_degC <- coral_specs_reg_mean %>% 
  as_tibble() %>% 
  #select(spec_id, paleoData_variableName, freq, spec) %>% 
  pivot_longer(cols = c(spec, starts_with("lim."), mad_e, se)) %>% 
  mutate(value = ifelse(spec_id != "SignalNoise",
                        ifelse(paleoData_variableName == "d18O", 
                               value * 1/0.22^2,
                               value * 1/0.06^2), value)
         ) %>% 
  pivot_wider() %>% 
  as_spec_df() 


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs_distsub <- coral_specs_reg %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  ungroup() %>% 
  select(pair, paleoData_variableName, dist_km,
         paleoData_TSid_a, paleoData_TSid_b, spec_id, freq, spec) %>% 
  rename(a = paleoData_TSid_a, b = paleoData_TSid_b) %>% 
  pivot_longer(cols = c(a, b), names_to = "member",
               values_to = "paleoData_TSid")


## ----------------------------------------------------------------------------------------------------------------------------
BootAndSummarise <- function(dat) {
  
  boot_sample_data <- dat %>%
    ungroup() %>%
    select(paleoData_TSid) %>%
    reframe(used_records = unique(c(paleoData_TSid))) %>%
    reframe(
      boot_record_sample = sample(used_records, n(), replace = TRUE),
      i = 1:n()
    ) %>%
    left_join(., dat, by = c("boot_record_sample" = "paleoData_TSid"),
              relationship = "many-to-many")

  boot_sample_summary <- boot_sample_data %>%
    group_by(spec_id, paleoData_variableName, freq) %>%
    summarise(
      spec_mean = .Internal(mean(spec)),
      spec_median = ccaPP::fastMedian(spec),
      .groups = "keep"
    ) %>%
    pivot_longer(starts_with("spec_m"), values_to = "spec",
                 names_to = "statistic")
 
  return(boot_sample_summary)
}


## ----load_boot, cache.lazy = FALSE-------------------------------------------------------------------------------------------
if (file.exists("..data/coral_specs_reg_mean_boot.RData") == FALSE | file.exists("..data/coral_specs_reg_mean_boot_degC.RData") == FALSE |
    use_cache == FALSE) {
  if (file.exists("../data/coral_pair_specs_boot.RData") == FALSE | 
      use_cache == FALSE) {
  n_reps <- n_boot
  pb <- progress_bar$new(total = n_reps)

  foo <- function(x) {
    pb$tick(0)
    y <- BootAndSummarise(x)
    pb$tick()
    return(y)
  }
  system.time({
    boot_spec <- tibble(
      boot_rep = 1:n_reps
    ) %>%
      group_by(boot_rep) %>%
      reframe(
        # suppressMessages(
        foo(coral_specs_distsub)
        # )
      )

    saveRDS(boot_spec, file = "../data/coral_pair_specs_boot.RData")
  })
} else {
  boot_spec <- readRDS(file = "../data/coral_pair_specs_boot.RData")
}

  
coral_specs_reg_mean_boot <- boot_spec %>% 
  group_by(spec_id, paleoData_variableName, statistic, freq) %>% 
  hamstr:::summarise_q(., var = spec)

coral_specs_reg_mean_boot_degC <- coral_specs_reg_mean_boot %>%
  as_tibble() %>% 
  select(spec_id, paleoData_variableName, statistic,
         freq, mean, contains("%")) %>% 
  pivot_longer(cols = c(mean, contains("%"))) %>%
  mutate(value = ifelse(spec_id != "SignalNoise",
                        ifelse(paleoData_variableName == "d18O",
                        value * 1/0.22^2,
                        value * 1/0.06^2), value)
         ) %>%
  pivot_wider() %>%
  as_spec_df() %>% 
  left_join(., 
            (coral_specs_reg_mean_degC))

saveRDS(coral_specs_reg_mean_boot,
        "../data/coral_specs_reg_mean_boot.RData")

saveRDS(coral_specs_reg_mean_boot_degC,
        "../data/coral_specs_reg_mean_boot_degC.RData")

} else {
  
  coral_specs_reg_mean_boot <- 
    readRDS("../data/coral_specs_reg_mean_boot.RData")
 
   coral_specs_reg_mean_boot_degC <- 
     readRDS("../data/coral_specs_reg_mean_boot_degC.RData")
  
}


## ----------------------------------------------------------------------------------------------------------------------------
pairs_per_freq <- coral_specs_reg %>%
  filter(pair %in% pairs_to_use$pair) %>%
  filter(complete.cases(spec)) %>% 
  group_by(spec_id, paleoData_variableName, freq) %>%
  summarise(mean_dist_km = mean(dist_km),
            n_pairs = n(),
            n_records = n_distinct(c(paleoData_TSid_a, paleoData_TSid_b)),
            .groups = "keep") %>% 
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
fig_no_pairs <- pairs_per_freq %>% 
  filter(spec_id == "S_clim") %>% 
  ggplot(aes(x = freq)) +
  geom_line(aes(y = n_pairs, linetype = "No. pairs")) +
  geom_line(aes(y = n_records, linetype = "No. unique records")) +
  scale_x_log10("Frequency [1/years]", 
                sec.axis = sec_axis(transform = ~1/.,
                                    name = "Timescale [years]")) +
  #scale_y_continuous(breaks = c(3, 10, 20, 30, 40)) +
  expand_limits(y = c(0, 70)) +
  facet_wrap(~paleoData_variableName, labeller = Spec_SrCa_d18O_labeller) +
  annotation_logticks(sides = "bt") +
  theme_bw() +
  labs(x = "Frequency", y = "No. of pairs", colour = "", linetype = "") +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right")

fig_no_pairs


## ----------------------------------------------------------------------------------------------------------------------------
PlotSNRStackBoot <- function(coral_dat, SST_dat, 
                             components = c("S_clim", "S_proxy", "S_noise"),
                             f_range = NULL,
                             statistic = c("spec_mean", "spec_median"),
                             plot_SST = TRUE) {
  statistic_type <- match.arg(statistic)

  if (is.null(f_range)) {
    f_range <- range(coral_dat$freq)
  }

  coral_dat_sub <- coral_dat %>%
    filter(statistic == statistic_type) %>%
    filter(
      freq >= f_range[1],
      freq <= f_range[2]
    ) %>%
    mutate(type = spec_id, spec = spec) %>%
    filter(type %in% components) %>%
    mutate(spec_id = type) %>%
    as_spec_df()

  SST_dat_sub <- SST_dat %>%
    filter(
      freq >= f_range[1],
      freq <= f_range[2]
    )

  p <- coral_dat_sub %>%
    gg_spec(., min.colours = 1, quantiles = TRUE, conf = FALSE,
            force.lims = TRUE, time_unit = "years") +
    facet_wrap(~paleoData_variableName,
      labeller = Spec_SrCa_d18O_labeller # , scales = "free_y"
    ) +
    scale_linetype_manual(values = c(lty_SNR_components, lty_coral_proxies)) +
    labs(colour = "", fill = "", x = "Frequency [1/years]", 
         y = expression(PSD ~ "[K"^
      {
        2
      } * "yr]")) +
    scale_colour_manual(
      values = c(pal_coral_proxies, pal_SNR_components),
      breaks = names(c(pal_SNR_components, pal_coral_proxies)),
      labels = lbl_SNR_coral_proxies,
      aesthetics = c("colour", "fill")
    )

  if (plot_SST) {
    p <- p +
      geom_ribbon(
        data = SST_dat_sub, aes(
          x = freq, ymax = X81.4., ymin = X15.9.,
          fill = SST_product,
          group = paste(SST_product, spec_id)
        ),
        alpha = 1 / 3, colour = NA
      ) +
      geom_line(data = SST_dat_sub, 
                aes(x = freq, y = spec, colour = SST_product,
                    group = paste(SST_product, spec_id)), lty = 2)
  }
  return(p)
}


## ----------------------------------------------------------------------------------------------------------------------------
PlotSNRStackBoot_byVar <- function(
    coral_dat, SST_dat, 
    components = c("S_clim", "S_proxy", "S_noise"),
    f_range = NULL,
    statistic = c("spec_mean", "spec_median"),
    plot_SST = TRUE) {
  statistic_type <- match.arg(statistic)

  if (is.null(f_range)) {
    f_range <- range(coral_dat$freq)
  }

  coral_dat_sub <- coral_dat %>%
    filter(statistic == statistic_type) %>%
    filter(
      freq >= f_range[1],
      freq <= f_range[2]
    ) %>%
    mutate(type = spec_id, spec = mean) %>%
    filter(type %in% components) %>%
    mutate(spec_id = paleoData_variableName) %>%
    as_spec_df()

  SST_dat_sub <- SST_dat %>%
    filter(
      freq >= f_range[1],
      freq <= f_range[2]
    ) %>%
    crossing(., type = c("S_clim", "S_proxy"))

  p <- coral_dat_sub %>%
    gg_spec(., min.colours = 1, quantiles = TRUE, conf = FALSE,
            force.lims = TRUE, time_unit = "years") +
    facet_wrap(~ factor(type, levels = rev(c("S_clim", "S_noise", "S_proxy"))),
      labeller = Spec_SrCa_d18O_labeller, as.table = FALSE
    ) +
    scale_linetype_manual(
      values = c(lty_SNR_components, lty_coral_proxies, "SST" = 2),
      breaks = names(c(lty_SNR_components, lty_coral_proxies, "SST" = 2)),
      labels = names(c(lty_SNR_components, lty_coral_proxies, "SST" = 2))
    ) +
    labs(
      colour = "", fill = "", y = expression(PSD ~ "[K"^
        {
          2
        } * "yr]"),
      linetype = ""
    ) +
    scale_colour_manual(
      values = c(pal_coral_proxies, pal_SNR_components),
      breaks = names(c(pal_SNR_components, pal_coral_proxies)),
      labels = lbl_SNR_coral_proxies,
      aesthetics = c("colour", "fill")
    )

  if (plot_SST) {
    p <- p +
      geom_ribbon(
        data = SST_dat_sub, aes(
          x = freq, ymax = X81.4., ymin = X15.9.,
          fill = paste(paleoData_variableName),
          # group = paste(SST_product, spec_id, paleoData_variableName)
          # fill = "SST",
          group = paste(SST_product, spec_id, paleoData_variableName)
        ),
        alpha = 1 / 3, colour = NA
      ) +
      geom_line(data = SST_dat_sub, aes(
        x = freq, y = spec,
        colour = paste(paleoData_variableName),
        linetype = "SST",
        # colour = "SST",
        group = paste(SST_product, spec_id, paleoData_variableName)
      ))
  }
  return(p)
}


## ----warning=FALSE-----------------------------------------------------------------------------------------------------------
fig_SNR_full_Freq <- PlotSNRStackBoot(
  coral_specs_reg_mean_boot_degC,
  ERSST_specs_reg_mean
  )
fig_SNR_full_Freq


## ----warning=FALSE-----------------------------------------------------------------------------------------------------------
fig_S_Clim_Instr_degC <- PlotSNRStackBoot(coral_specs_reg_mean_boot_degC,
                                          ERSST_specs_reg_mean,
                                          statistic = "spec_mean",
                                          components = c("S_clim", "S_proxy"),
                                          f_range = c(1 / 100, 1 / 2)
                                          )


## ----warning=FALSE-----------------------------------------------------------------------------------------------------------
fig_S_components_degC <- PlotSNRStackBoot_byVar(coral_specs_reg_mean_boot_degC,
                                                ERSST_specs_reg_mean,
                                                f_range = c(1/100, 1/2),
                                                plot_SST = FALSE) +
  coord_cartesian(ylim = c(0.01, 2))



## ----fig_Clim_Instr_degC_No_Comp, fig.height=8, warning=FALSE----------------------------------------------------------------
patchwork::wrap_plots(
  fig_S_Clim_Instr_degC + 
    coord_cartesian(xlim = c(0.01, 1/2), ylim = c(1e-02, 3e00))+ 
    theme(axis.text.x.bottom = element_blank(),
          axis.title.x.bottom = element_blank()) +
  labs(tag = "a"),
  fig_no_pairs + 
    coord_cartesian(xlim = c(0.01, 1/2))+ 
    theme(axis.text.x.top = element_blank(),
          axis.title.x.top = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.spacing = unit(0, "mm"),
          legend.margin = unit(0, "mm")) +
  labs(tag = "b"),
  fig_S_components_degC,
  ncol = 1, heights = c(4, 2, 4)) +
  labs(tag = "c")


## ----------------------------------------------------------------------------------------------------------------------------
GetPairVarFband <- function(pair_specs){
  pair_specs %>% 
  filter(spec_id != "SignalNoise") %>% 
  group_by(pair, dist_km) %>% 
  mutate(min_f = min(freq)) %>% 
  ungroup() %>% 
  crossing(., fbands) %>%
  filter(min_f <= f_lwr) %>% 
  group_by(pair, dist_km, paleoData_variableName,
           paleoData_TSid_a, paleoData_TSid_b, spec_id,
           f_band, f_band_name) %>% 
  do({
   sp1 <- DF2Spec(.)
   var <- GetVarFromSpectra(sp1,
                            f = c(unique(.$f_lwr), unique(.$f_upr)), bw = 1)
   data.frame(var_proxy = var$var, dof = var$dof)
  }) %>% 
  pivot_wider(names_from = spec_id, values_from = c(var_proxy, dof)) %>% 
  group_by(pair, dist_km, f_band, f_band_name, paleoData_variableName) %>% 
  mutate(var_ratio = var_proxy_S_clim / var_proxy_S_proxy,
         SNR = var_proxy_S_clim / var_proxy_S_noise) 
}

SummariseVarFband <- function(var_fband){
  var_fband %>%
    group_by(f_band, f_band_name, paleoData_variableName) %>% 
    summarise(var_ratio_mu = mean(var_ratio, na.rm = TRUE),
            var_proxy_S_clim = mean(var_proxy_S_clim, na.rm = TRUE),
            var_proxy_S_proxy = mean(var_proxy_S_proxy, na.rm = TRUE),
            SNR = mean(SNR, na.rm = TRUE),
            dof_S_clim = sum(dof_S_clim),
            dof_S_proxy = sum(dof_S_proxy),
            .groups = "keep") 
}


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs_var_fband <- coral_specs %>%
  filter(resolution == "mon") %>% 
  GetPairVarFband(.)


coral_specs_var_fband_mean <- coral_specs_var_fband %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  SummariseVarFband(.) 


## ----------------------------------------------------------------------------------------------------------------------------
BootAndSummariseVarFband <- function(dat) {
 
  boot_sample_data <- dat %>%
    ungroup() %>%
    select(paleoData_TSid) %>%
    reframe(used_records = unique(c(paleoData_TSid))) %>%
    reframe(
      boot_record_sample = sample(used_records, n(), replace = TRUE),
      i = 1:n()
    ) %>%
    left_join(., dat, by = c("boot_record_sample" = "paleoData_TSid"),
              relationship = "many-to-many")

  boot_sample_summary <- SummariseVarFband(boot_sample_data)

  return(boot_sample_summary)
}


## ----------------------------------------------------------------------------------------------------------------------------
coral_specs_varsub <- coral_specs_var_fband %>% 
   filter(pair %in% pairs_to_use$pair) %>% 
  ungroup() %>% 
  rename(a = paleoData_TSid_a, b = paleoData_TSid_b) %>% 
  pivot_longer(cols = c(a, b),
               names_to = "member", 
               values_to = "paleoData_TSid")


## ----------------------------------------------------------------------------------------------------------------------------
if (file.exists("../data/coral_pair_var_boot.RData") == FALSE |
    use_cache == FALSE) {
  n_reps <- n_boot
  pb <- progress_bar$new(total = n_reps)

  foo2 <- function(x) {
    pb$tick(0)
    y <- BootAndSummariseVarFband(x)
    pb$tick()
    return(y)
  }

  system.time({
    boot_var_fband <- tibble(
      boot_rep = 1:n_reps
    ) %>%
      group_by(boot_rep) %>%
      reframe(
        foo2(coral_specs_varsub)
        )

    saveRDS(boot_var_fband,
      file = "../data/coral_pair_var_boot.RData"
    )
  })
} else {
  boot_var_fband <- readRDS(
    file = "../data/coral_pair_var_boot.RData"
  )
}





## ----------------------------------------------------------------------------------------------------------------------------
boot_var_fband_summary <- boot_var_fband %>% 
  group_by(f_band, f_band_name, paleoData_variableName) %>% 
  summarise_q_2(., var = 1/var_ratio_mu)


## ----fig_inflation_factor, fig.width = 5-------------------------------------------------------------------------------------
boot_var_fband_summary %>% 
  ggplot(aes(x = f_band_name, y = `50%`, colour = paleoData_variableName)) +
  geom_linerange(aes(ymax = `97.5%`, ymin = `2.5%`),
                 alpha = 1, position = position_dodge(width = 1/4)) +
  geom_pointrange(aes(ymax = `84.1%`, ymin = `15.9%`),
                 lwd = 1.2, alpha = 1, position = position_dodge(width = 1/4)) +
  #geom_point() +
  geom_hline(yintercept = 1, lty = 2) +
  theme_bw() +
  expand_limits(y = c(1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank())+
  #facet_grid(.~., scales = "free_y") +
  labs(x = "Timescale", y = "Inflation factor") +
  scale_colour_manual("", values = pal_coral_proxies,
                      labels = lbl_SNR_coral_proxies) +
  annotation_logticks(sides = "l") +
  scale_y_log10(sec.axis = sec_axis(name = "Corrected / Raw Variance", 
                                    trans = ~1/.,
                                    breaks = c(0.033, 0.1, 0.33, 1))) +
  theme(legend.position = "top", panel.grid.major.x = element_blank())



## ----------------------------------------------------------------------------------------------------------------------------
boot_var_fband_summary 


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_binned_SST <- coral_pairs_sub %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
   mutate(month_CE = year * 12 + month,
         bin_width = ifelse(resolution == "mon", 1, 12)) %>% 
  group_by(resolution, pair, dist_km, paleoData_variableName) %>% 
  mutate(min.age = min(month_CE), max.age = max(month_CE)) %>%
  pivot_longer(cols = c(SST_OISST2_dis, SST_ERSST5_dis),
               names_to = "SST_product", values_to = "SST") %>% 
  group_by(resolution, pair, dist_km, paleoData_variableName,
           paleoData_TSid, SST_product) %>%
  arrange(month_CE) %>%
  reframe( 
    BinTimeseries(month_CE, SST,
                  bin.width = bin_width[1],
                  strt.time = min.age[1], end.time = max.age[1]
                  )
    ) %>%
  group_by(resolution, pair, dist_km, paleoData_variableName,
           paleoData_TSid, SST_product) %>%
  mutate(n.pts = sum(is.na(mean.value) == FALSE)) %>%
  rename(age = time) %>%
  filter(n.pts > 0) %>%
  ungroup() %>% 
  group_by(resolution, pair, paleoData_variableName) %>% 
  mutate(p_NA = sum(is.na(mean.value)) / n())

coral_pairs_binned_wide_SST <- coral_pairs_binned_SST %>%
  filter(p_NA < 0.2) %>% 
  ungroup() %>%
  select(-n.bin, -n.pts, -mean.time) %>%
  group_by(resolution, bin.width, pair, dist_km, 
           paleoData_variableName, SST_product) %>%
  do({
    m = pivot_wider(., names_from = paleoData_TSid,
                    values_from = mean.value,
                    names_prefix = "record_") %>%
     arrange(age) %>%
     select(starts_with("record_")) %>%
     as.matrix(.)
    
     colnames(m) <- NULL
    
    tibble(
      pair_size = ncol(m),
      n_gaps = sum(is.na(m)),
      p_gaps = n_gaps / length(m),
      m = list(m)
    )
  })

coral_pairs_binned_wide_SST %>% 
  group_by(paleoData_variableName, SST_product) %>% 
  summarise(n_pairs = n_distinct(pair))


## ----------------------------------------------------------------------------------------------------------------------------
SST_specs <- coral_pairs_binned_wide_SST %>% 
  #filter(p_gaps < 0.2) %>% 
  group_by(bin.width, pair, resolution, dist_km, pair_size, 
           paleoData_variableName, SST_product) %>% 
  reframe(
      Spec2DF(SNRStack(x = m[[1]], 
                       k = 5, nw = 4,
                        logsmooth = FALSE,
                        equalise_var = eq_var,
                          bin_width = bin.width[1]/12)) 
  ) %>% 
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
SST_specs_reg <- SST_specs %>% 
  group_by(spec_id, bin.width, pair, resolution, dist_km, 
           paleoData_variableName, SST_product) %>% 
  do({
    dat <- .
    min_f <- min(dat$freq)
    
    strt_f <- if (dat$SST_product[1] == "SST_OISST2_dis") 0.1 else {0.01}
    maxf <- 6
    rfreq = seq(strt_f, maxf, strt_f/10)
    
    sr <- as.spec(list(freq = dat$freq, spec = dat$spec, dof = dat$dof))
    
    si <- PaleoSpec::SpecInterpolate(sr,
                                     freqRef = rfreq,
                                     check = FALSE)
    
   tibble(freq = si$freq[si$freq >= min_f], spec = si$spec[si$freq >= min_f],
          dof = si$dof[si$freq >= min_f])
  }) 

SST_specs_reg_mean <- SST_specs_reg %>% 
  # remove zero distance pairs
  
  filter(dist_km > 0) %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  group_by(spec_id, paleoData_variableName, SST_product, freq) %>% 
  summarise(n = sum(is.na(spec) == FALSE),
            se = sd(spec, na.rm = TRUE) / sqrt(n),
            spec = mean(spec, na.rm = TRUE)
               ) %>% 
  mutate(lim.1 = spec + se, lim.2 = spec - se) %>% 
  as_spec_df()


## ----------------------------------------------------------------------------------------------------------------------------
c_specs_var <- coral_specs %>% 
  filter(spec > 0, freq >= 0.01) %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  filter(spec_id %in% c("S_clim", "S_proxy")) %>% 
  group_by(bin.width, pair, pair_size, paleoData_variableName, spec_id) %>% 
  do({
   sp1 <- DF2Spec(.)
   
   #print(min(sp1$freq))
   var <- PaleoSpec::GetVarFromSpectra(sp1, f = c(min(sp1$freq), 1/2))$var
   
   tibble(var_proxy = var, n_freq = length(sp1$freq))
   
  }) %>% 
  mutate(sd_proxy = sqrt(var_proxy))

c_specs_var_SST <- SST_specs %>% 
  filter(spec > 0, freq >= 0.01) %>% 
  filter(pair %in% pairs_to_use$pair) %>% 
  filter(spec_id %in% c("S_clim", "S_proxy")) %>% 
  group_by(bin.width, pair, pair_size, paleoData_variableName,
           SST_product, spec_id) %>% 
  do({
   sp1 <- DF2Spec(.)
   
   # print(.$pair[1])
   # print(.$SST_product[1])
   # print(min(sp1$freq))
   var <- PaleoSpec::GetVarFromSpectra(sp1, f = c(min(sp1$freq), 1/2))$var
   
   tibble(var_SST = var)
   
  }) %>% 
  mutate(sd_SST = sqrt(var_SST))

c_specs_var_degC <- c_specs_var %>% 
  mutate(var_proxy = ifelse(spec_id != "SignalNoise",
                            ifelse(paleoData_variableName == "d18O", 
                                   var_proxy * 1/0.22^2,
                                   var_proxy * 1/0.06^2), var_proxy)
         )

c_specs_var_comb <- left_join(c_specs_var_degC, 
                              filter(c_specs_var_SST,
                                     spec_id == "S_proxy") %>% 
                                ungroup() %>% 
                                select(-spec_id) ,
                              relationship = "many-to-many"
          )


comp_names <- tibble(
  spec_id = c("S_clim", "S_noise", "S_proxy", "S_stack"),
  name = factor(c("Corrected", "Noise", "Uncorrected", "Stack"), ordered = TRUE, 
                levels = c("Uncorrected", "Corrected", "Noise", "Stack"))
)


## ----fig.width=7, fig.height=5-----------------------------------------------------------------------------------------------
c_specs_var_comb %>% 
  left_join(., comp_names) %>% 
  mutate(slope = ifelse(paleoData_variableName == "d18O", 0.22, 0.06)) %>% 
  filter(complete.cases(var_SST),
         spec_id %in% c("S_proxy", "S_clim"), 
         SST_product == "SST_OISST2_dis") %>% 
  ggplot(aes(x = sqrt(var_SST), y = sqrt(var_proxy),
             colour = paleoData_variableName)) +
  geom_point() +
  # geom_abline(aes(slope = slope, intercept = 0),
  #             lty = 2) +
  geom_abline(slope = 1, intercept = 0,
              lty = 2) +
  theme_bw() +
  #scale_x_sqrt() +
  #scale_y_sqrt() +
  scale_colour_manual(values = pal_coral_proxies,
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "SD SST [°C]", y = "SD Sr/Ca [mmol/mol]",
       colour = "") +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        #strip.background = element_blank(),
        legend.position = "none"#,
        #strip.text = element_blank()
        ) +
  coord_equal() +
  expand_limits(x = c(0, 1), y = c(0,1)) +
  labs(x = "SD OISSTv2 [°C]", y = "SD scaled coral proxy [°C]", colour = "") +
  facet_grid(name ~ paleoData_variableName, labeller = Spec_SrCa_d18O_labeller)


## ----------------------------------------------------------------------------------------------------------------------------
c_specs_var_comb %>% 
  left_join(., comp_names) %>% 
  mutate(slope = ifelse(paleoData_variableName == "d18O", 0.22, 0.06)) %>% 
  filter(complete.cases(var_SST),
         spec_id %in% c("S_proxy", "S_clim"), 
         SST_product == "SST_ERSST5_dis") %>% 
  ggplot(aes(x = sqrt(var_SST), y = sqrt(var_proxy), 
             colour = paleoData_variableName)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,
              lty = 2) +
  theme_bw() +
  scale_colour_manual(values = pal_coral_proxies, 
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "SD SST [°C]", y = "SD Sr/Ca [mmol/mol]",
       colour = "") +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none"#,
        ) +
  coord_equal() +
  expand_limits(x = c(0, 1), y = c(0,1)) +
  labs(x = "SD ERSSTv5 [°C]", y = "SD scaled coral proxy [°C]", colour = "") +
  facet_grid(name ~ paleoData_variableName, labeller = Spec_SrCa_d18O_labeller)


## ----------------------------------------------------------------------------------------------------------------------------
fig_SNR_dist_SrCa <-  coral_specs_var_fband %>% 
  filter(SNR > 0) %>% 
  filter(paleoData_variableName == "SrCa") %>% 
  filter(f_band_name != "Decadal to\nCentennial") %>% 
  ggplot(aes(x = dist_km+0.1, y = (SNR),
             colour = f_band_name, fill = f_band_name)
         ) +
  geom_point() +
  geom_vline(xintercept = max_dist, linetype = 2) +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(breaks = c(0, 1, 10, 100, 1000)) +
  geom_smooth(method = "lm") +
  annotation_logticks(sides = "lb") +
  labs(y = "Signal to Noise ratio", 
       fill = "", colour = "",
       x = "Mean pairwise distance [km+0.1]") +
  scale_color_viridis_d(aesthetics = c("fill", "colour"), end = 0.9) +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) 

fig_SNR_dist_SrCa


## ----------------------------------------------------------------------------------------------------------------------------
coral_pairs_sub_lm <- coral_pairs_once %>%
  filter(pair %in% pairs_to_use$pair) %>% 
  select(year, month, starts_with("paleoData_"),
         SST_OISST2_dis, SST_ERSST5_dis, proxy) %>% 
  pivot_longer(cols = c("SST_ERSST5_dis", "SST_OISST2_dis"),
               names_to = "SST_product", values_to = "SST") %>% 
  group_by(paleoData_variableName, paleoData_TSid, SST_product) %>%
  mutate(n_sst = sum(is.na(SST) == FALSE)) %>%
  filter(n_sst > 10) %>%
  do({
    dat <- .

    lm1 <- lm(proxy~SST, na.action = na.exclude, data = dat)

    tibble(dat, fttd_proxy_lm = fitted(lm1), resid_proxy_lm = residuals(lm1))

  })


## ----------------------------------------------------------------------------------------------------------------------------
regr_resids <- coral_pairs_sub_lm %>% 
  summarise_if(is.numeric, sd, na.rm = TRUE) %>% 
  mutate(resid_proxy_lm = ifelse(paleoData_variableName == "SrCa",
                                 resid_proxy_lm*1/0.06,
                                 resid_proxy_lm * 1/0.22))


## ----------------------------------------------------------------------------------------------------------------------------
regr_resids %>% 
  group_by(paleoData_variableName) %>% 
  hamstr:::summarise_q(., resid_proxy_lm)


## ----------------------------------------------------------------------------------------------------------------------------
## Bin timeseries to insert gaps
bin.width <- 1
coral_pairs_sub_lm_binned <- coral_pairs_sub_lm %>%
  mutate(month_CE = year * 12 + month) %>%
  group_by(paleoData_variableName, paleoData_TSid,
           paleoData_ch2kCoreCode, SST_product) %>%
  mutate(mu_dt = round(median(diff(month_CE)))) %>%
  arrange(month_CE) %>%
  select(paleoData_variableName, paleoData_TSid, 
         paleoData_ch2kCoreCode, SST_product,
         month_CE, mu_dt, resid_proxy_lm, resid_proxy_lm,
         proxy, fttd_proxy_lm) %>%
  pivot_longer(cols = c(resid_proxy_lm,
                        resid_proxy_lm, proxy,
                        fttd_proxy_lm), names_to = "type") %>%
  group_by(paleoData_variableName, paleoData_ch2kCoreCode,
           paleoData_TSid, SST_product, type) %>%
  reframe(
    #BinTimeseries(month_CE, value, bin.width = 2*mu_dt[1])
    BinTimeseries(month_CE, value, bin.width = bin.width)
  )


## ----------------------------------------------------------------------------------------------------------------------------
resid.specs <- coral_pairs_sub_lm_binned %>%
  group_by(paleoData_variableName, paleoData_TSid,
           paleoData_ch2kCoreCode, SST_product, type) %>%
  arrange(time) %>%
  reframe(Spec2DF(FilterSpecLog(
    FilterSpec(SpecACF(mean.value, bin.width = unique(bin.width)/12),
               spans = c(3))
    ))
    ) %>%
  mutate(spec_id = type) %>%
  mutate(n_neg = sum(spec <= 0)) %>%
  ungroup() %>%
  select(-type)


## ----------------------------------------------------------------------------------------------------------------------------
resid.specs_w <- resid.specs %>%
  group_by(paleoData_TSid, SST_product) %>%
  filter(rank(freq) > 2) %>%
  select(starts_with("paleoData"), SST_product, spec_id, freq, spec) %>%
  pivot_wider(names_from = spec_id, values_from = spec)%>%
  mutate(SNR = fttd_proxy_lm / resid_proxy_lm) %>%
  mutate(spec = SNR) %>%
  mutate(spec_id = paleoData_variableName)


## ----------------------------------------------------------------------------------------------------------------------------
resid.specs_w_reg <- resid.specs_w %>%
  mutate(dof = 1000) %>%
  group_by(paleoData_variableName, paleoData_TSid,
           paleoData_ch2kCoreCode, SST_product) %>%
  do({
    minf <- 0.01 # min(bin.ind.spec.df$freq)
    maxf <- max(resid.specs_w$freq)
    rfreq = seq(minf, maxf, minf*2)
    SNR <- PaleoSpec::SpecInterpolate(
      list(freq = .$freq, spec = .$SNR, dof = .$dof),
                                     freqRef = rfreq)

    S_proxy <- PaleoSpec::SpecInterpolate(
      list(freq = .$freq, spec = .$proxy, dof = .$dof),
                                     freqRef = rfreq)

    S_resids <- PaleoSpec::SpecInterpolate(
      list(freq = .$freq, spec = .$resid_proxy_lm, dof = .$dof),
                                                   freqRef = rfreq)

    S_fttd <- PaleoSpec::SpecInterpolate(
      list(freq = .$freq, spec = .$fttd_proxy_lm, dof = .$dof),
                                                   freqRef = rfreq)

    tibble(freq = SNR$freq, SNR = SNR$spec,
           S_proxy = S_proxy$spec,
           S_resids=S_resids$spec,
           S_fttd = S_fttd$spec)
  }) %>%
  mutate(S_fttd_plus_resids = S_fttd + S_resids) %>%
  pivot_longer(cols = starts_with("S_"), values_to = "spec")


resid.specs_w_summary <- resid.specs_w_reg  %>%
  mutate(spec = ifelse(paleoData_variableName == "SrCa",
                       spec*1/0.06^2, spec * 1/0.22^2)) %>%
  group_by(paleoData_variableName, SST_product, name, freq) %>%
  summarise(n = sum(is.na(spec)==FALSE),
            se = sd(spec, na.rm = TRUE) / sqrt(n),
            mad_e = 1.2533 * mad(spec, na.rm = TRUE) / sqrt(n),
            spec_mean = mean(spec, na.rm = TRUE),
            spec_median = median(spec, na.rm = TRUE)) %>%
  pivot_longer(starts_with("spec_"), values_to = "spec",
               names_to = "statistic") %>%
  mutate(spec_id = name,
         `15.9%` = ifelse(statistic == "spec_mean", spec - se, spec - mad_e),
         `84.1%` = ifelse(statistic == "spec_mean", spec + se, spec + mad_e),
         `97.5%` = ifelse(statistic == "spec_mean", spec + 2*se, spec + 2*mad_e),
         `2.5%` = ifelse(statistic == "spec_mean", spec - 2*se, spec - 2*mad_e))


## ----------------------------------------------------------------------------------------------------------------------------
resids_spec <- resid.specs_w_summary %>%
  filter(spec_id == "S_resids", statistic == "spec_mean", freq > 1/100) %>% 
  mutate(spec_id = SST_product)


## ----fig_S_noise_resids, fig.width=8-----------------------------------------------------------------------------------------
regr_labs <- list(S_noise = expression("Proxy noise (" * S[noise] * ")"),
     SST_ERSST5_dis = expression("ERSSTv5 residuals"),
     SST_OISST2_dis = expression("OISSTv2 residuals"))

fig_SNR_vs_resids <- coral_specs_reg_mean_boot_degC %>%
  filter(spec_id == "S_noise", statistic == "spec_mean") %>%
  bind_rows(., resids_spec) %>%
  filter(`2.5%` > 0) %>% 
  as_spec_df() %>%
  gg_spec(., quantiles = TRUE, conf = FALSE, time_unit = "years") +
  facet_wrap(~paleoData_variableName,
              labeller = Spec_SrCa_d18O_labeller, scales = "free_y") +
  scale_colour_manual(values = c("S_noise" = "#d95f02",
                                 "SST_ERSST5_dis" = "#b2df8a",
                                 "SST_OISST2_dis" = "#1f78b4"),
                      breaks = names(regr_labs),
                      labels = regr_labs,
                      aesthetics = c("colour", "fill"))+
  labs(colour = "", fill = "", y = expression(PSD~"[K"^{2}*"yr]"))

fig_SNR_vs_resids


## ----------------------------------------------------------------------------------------------------------------------------
ex_dat <- readRDS("../data/FE18RUS01_SrCa.RDS")


## ----------------------------------------------------------------------------------------------------------------------------
# Empirically estimated noise spectrum
SrCa_err_spec <- coral_specs_reg_mean_boot_degC %>% 
  filter(paleoData_variableName == "SrCa", 
         spec_id == "S_noise", statistic == "spec_mean") %>% 
  as_spec_df() %>% 
  DF2Spec()

# Power-law noise spectrum
pl_spec <- tibble(
  freq = seq(1/1000, 6, 1/1000),
  spec = exp(-2.65) * freq^-0.599) %>%  
  as_spec_df(.) 

# Splice the spectra together and interpolated to regular frequency axis
spliced_spec <- bind_rows(
  pl_spec %>% 
    filter(freq < 1/100),
  SrCa_err_spec %>% 
    Spec2DF() %>% 
    filter(freq >= 1/100)
) %>% 
  as_tibble() %>% 
  select(freq, spec) %>% 
  mutate(dof = 20) %>% 
  as_spec_df() %>% 
  do({
    sp1 <- DF2Spec(.)
    
    sp2 <- PaleoSpec::SpecInterpolate(sp1, seq(1/1000, 6, 1/1000))
    
    Spec2DF(sp2)
  }) %>% 
  DF2Spec()

gg_spec(spliced_spec)


## ----------------------------------------------------------------------------------------------------------------------------
delta_t_ex <- mean(diff(ex_dat$year.dec))
tau_ex <- diff(range(ex_dat$year.dec))


sigsq_spec_dec <- GetVarFromSpectra(spliced_spec,
                                    f = c(1/tau_ex, 1/20))$var

sigsq_spec_all_white <- GetVarFromSpectra(spliced_spec,
                                          f = c(1/tau_ex, 1/(2*delta_t_ex)))$var / (10*1/delta_t_ex)


## ----------------------------------------------------------------------------------------------------------------------------
ex_dat_dec <- ex_dat %>% 
  mutate(decade = round(year.dec, -1)) %>% 
  group_by(decade) %>% 
  summarise_if(is.numeric, mean) %>% 
  mutate(SST_proxy = (paleoData_values - mean(paleoData_values)) / -1/0.06) %>% 
  select(decade, SST_proxy) %>% 
  mutate(sigma_auto = sqrt(sigsq_spec_dec), 
         sigma_inde = sqrt(sigsq_spec_all_white)) %>% 
  pivot_longer(starts_with("sigma_"),
               names_to = "noise_type", values_to = "sigma")

noise_type_labeller_2 <- ggplot2::as_labeller(
  c(
    sigma_meas = "Measurement error only",
    sigma_auto = "Autocorrelated noise (this study)",
    sigma_inde = "Independent noise"))


ex_dat_dec %>% 
  filter(noise_type != "sigma_meas") %>% 
  ggplot(aes(x = decade, y = SST_proxy)) +
  geom_line(aes(colour = noise_type)) +
  geom_ribbon(aes(
    ymax = SST_proxy + sigma,
    ymin = SST_proxy - sigma,
    fill = noise_type
    ), alpha = 0.5)+
  geom_ribbon(aes(
    ymax = SST_proxy + 2*sigma,
    ymin = SST_proxy - 2*sigma,
    fill = noise_type
    ), alpha = 0.25) +
  facet_grid(noise_type~., labeller = noise_type_labeller_2, as.table = FALSE) +
  scale_colour_manual("Noise type",
                      breaks = rev(c(
                        "sigma_auto",
                        "sigma_inde")),
                      values = c(sigma_inde = "#5ab4ac",
                                 sigma_auto = "#d8b365"),
                      label =  rev(c(
                        sigma_auto = "Autocorrelated (this study)",
                        sigma_inde = "Independent")),
                      aesthetics = c("fill", "colour")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")+
  labs(x = "Decade [CE]", y = "Decadal temperature anomaly [°C]")

