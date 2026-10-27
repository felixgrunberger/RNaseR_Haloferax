# Figure 2: Growth curves

library(here)
source(here("scripts/helpers.R"))

get_data2 <- function(which_data){
  example_data <- read_xlsx(here("data/growth/Neue Wachstumskurven.xlsx"),
                            sheet = which_data, skip = 1, col_names = F)
  
  example_data_wide <- example_data[c(1:16,18)] %>%
    janitor::clean_names() %>%
    dplyr::rename(strain = 1, time = x18) %>%
    fill(strain) %>%
    pivot_longer(cols = -c(strain, time)) %>%
    dplyr::rename(replicate = name,
                  od = value)  %>%
    dplyr::mutate(replicate = as.factor(as.numeric(str_remove_all(replicate, "x"))-1)) %>%
    dplyr::mutate(bio = ifelse(replicate %in% 1:5, 1,
                               ifelse(replicate %in% 6:10, 2,
                                      ifelse(replicate %in% 11:15, 3,NA)))) %>%
    dplyr::mutate(tech = ifelse(replicate %in% c(1,6,11), 1,
                                ifelse(replicate %in% c(2,7,12), 2,
                                       ifelse(replicate %in% c(3,8,13), 3,
                                              ifelse(replicate %in% c(4,9,14), 4,
                                                     ifelse(replicate %in% c(5,10,15), 5,NA)))))) %>% 
    drop_na() %>%
    dplyr::mutate(condition = which_data)
}

# data ----
st   <- get_data2("Standard")
st2  <- get_data2("Standard-2")
gr30 <- get_data2("30°C")
gr50 <- get_data2("50°C")
mm   <- get_data2("Minimal Medium")
all <- bind_rows(st, st2, gr30,gr50, mm)
# remove outliers ----
outlier_table <- data.table(replicate = c("15","3","5","14","6","4",
                                          "10", "9",
                                          "10",
                                          "9", "7", "10","8",
                                          "1", "10", "11",
                                          "8", "10", "9","12",
                                          "15", "14", "4", "5",
                                          "5","15","6",
                                          "15", "4", "6", "11", "5", "1","3", "2", "14", "13"),
                            condition = c(rep("Standard-2",6),
                                          rep("Standard-2",2),
                                          rep("30°C",1),
                                          rep("50°C",4),
                                          rep("Minimal Medium",3),
                                          rep("Standard",4),
                                          rep("50°C",4),
                                          rep("Minimal Medium",3),
                                          rep("Standard",10)),
                            strain = c(rep("H119",6),
                                       rep("∆Rnase R",2),
                                       rep("∆Rnase R",1),
                                       rep("∆Rnase R",4),
                                       rep("∆Rnase R",3),
                                       rep("∆Rnase R",4),
                                       rep("H119",4),
                                       rep("H119",3),
                                       rep("H119",10))) %>%
  dplyr::mutate(remove = "yes")
all_f <- all %>%
  left_join(outlier_table) %>%
  dplyr::filter(is.na(remove))
ribbon_data <- all_f %>%
  group_by(condition, strain, time) %>%
  summarise(median_od = mean(od),
            sd = sd(od)) %>%
  group_by(condition, time) %>%
  mutate(ymin = min(median_od),
         ymax = max(median_od))

# Create the plot
dir.create(here("figures"), showWarnings = FALSE)
pdf(here("figures/Fig2_growth.pdf"), height = 4,width = 8)
p <- ggplot(data = ribbon_data %>%
         dplyr::filter(time < 4000,
                       condition!="Standard-2") %>%
         dplyr::mutate(condition = factor(condition,
                                          levels = c("Standard", "Minimal Medium",
                                                     "30°C", "50°C"))), 
       aes(x = time / 60, 
       y = median_od, 
       color = strain)) +
  facet_grid(cols = vars(condition)) +
  geom_ribbon(aes(x = time/60, ymin = median_od-sd, ymax = median_od+sd, fill = strain), inherit.aes = FALSE, alpha = .1) +
  geom_ribbon(aes(x = time/60, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.5, fill = "grey") +
  geom_line(linewidth = 1) +
  theme_pubclean() +
  xlab("time (h)") +
  ylab("Optical density (a.u.)") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#B7AEE3", "#A4A4A4")) +
  scale_color_manual(values = c("#B7AEE3", "#A4A4A4")) 
print(p)
dev.off()
