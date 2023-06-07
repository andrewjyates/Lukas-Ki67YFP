

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

packages <- c("tidyverse", "ggplot2", "readr", "forcats", "cowplot")

tmp=lapply(packages, require, character.only = TRUE)


#Gathering all parameters: 
  
#1. From literature 
#Order: lower CI, value, upper CI 

ggplot_scientific_notation_axes_labels <- function(l) {
  # 'prettyNum' puts things into scientific notation when digit strings get past a certain length,
  # but keeps things  >0.0001 and < 10^7 in non-sci notation
  l <- prettyNum(l)
  # … but for axis labels, we want to switch to powers of 10 for anything >=10^3 or < 0.01
  # (though you decide)
  for(i in 1:length(l)){
    number=eval(parse(text=l)[i])
    if(!is.na(number) & (abs(number)<0.001 | abs(number)>=1000)) l[i]=format(number, scientific=TRUE)
  }
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  l <- gsub("e0","e",l)  
  l <- gsub("e-0","e-",l) 
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # convert 0.000.. * 10^x to 0
  l <- gsub("\\'0[\\.0]*\\'\\%\\*\\%", "", l)
  l <- gsub("10\\^0", "0", l)
  # return this as an expression
  parse(text=l)
}

plotwidth=12


literature_decay <- list(
"CD4" = data.frame("upper.97.5." =c( NA, 42, 35.7, 90.9,   45.5,25.6, 28),
                  "lower.2.5." = c( NA, 28, 12.35,  52.6,    8.3, 16.7,18),
                  "value" =      c( 47, 34, 21.73, 66.7,  22.7, 20,22),
                  "receptor" = rep("CD4", 7),
                  "Method" = c("Den Braber 2012", "Hogan 2015", "van Hoeven 2017 (RTE)",  "van Hoeven 2017 (MN)","Rane 2018 (WT)","Rane 2018 (Chimeras)", "Rane 2022")
),
"CD8" = data.frame("upper.97.5." = c(NA,64, NA, 90.9,   38.5, 66.7,46),
                   "lower.2.5." = c(NA, 41, NA,  58.8,    10.3,  40, 34),
                   "value" =       c(80,55, NA,  71.4,    15.87, 52.3,40),
                   "receptor" = rep("CD8", 7),
                   "Method" = c("Den Braber 2012", "Hogan 2015", "van Hoeven 2017 (RTE)",  "van Hoeven 2017 (MN)","Rane 2018 (WT)","Rane 2018 (Chimeras)", "Rane 2022")
))
  # CD4 RTE  0.046 (CI = 0.028–0.081), MN 0.015 (CI = 0.011–0.019) 
# CD8 homog  0.014 (CI = 0.011–0.017)

literature_decay<- rbind(literature_decay[[1]], literature_decay[[2]])

#the 2.5% percentile is the lower limit and the 97.5% percentile is the upper limit of the 95% CI.
CIs <- function(y){
  map(y, function(x)
  map(names(x[["par"]]),  function(a) c("upper" = quantile(x[["bootstrap"]][,a], 0.975),
  "lower" = quantile(x[["bootstrap"]][,a], 0.025),
  "value" = x[["par"]][[a]]
  )) %>% setNames(names(x[["par"]]))
  )
}

shapiros <- function(x){
  map(x, ~shapiro.test(.[["residuals"]]))
}

# CIs_mu <- function(y){
#   map(y, function(x)
#     c("upper" = 1/quantile(x[["bootstrap"]][,"mu"], 0.975),
#       "lower" = 1/quantile(x[["bootstrap"]][,"mu"], 0.025),
#       "value" = 1/x[["par"]][["mu"]]
#     )
#   )
# }

CIs_mu <- function(y){
  map(y, function(x)
    c("upper" = quantile(1/(x[["bootstrap"]][,"mu"]), 0.975),
      "lower" = quantile(1/(x[["bootstrap"]][,"mu"]), 0.025),
      "value" = 1/x[["par"]][["mu"]]
    )
  )
}

#2. From separate models 

#2a. Simple model 
load("../fitting/DataObjects/simple_abs.RData")
load("../fitting/DataObjects/simple_frac.RData")


#2b. Flow model 
# for the flow model with fixed thymic inputs, p was non zero only for CD4
load("../fitting/DataObjects/flow_frac_pzero.RData")
load("../fitting/DataObjects/flow_frac_pfree.RData")

# combine the correct fits into one list
flow_frac = list(CD4 = flow_frac_pfree[["CD4"]], CD8 = flow_frac_pzero[["CD8"]])


#3. From simultaneous fits 

#3a. Simple model 
load("../fitting/DataObjects/simple_sim_abs.Rdata")
load("../fitting/DataObjects/simple_sim_frac.Rdata")

#3b. Flow model 
load("../fitting/DataObjects/flow_sim_frac_pzero.Rdata") # for flow sim models, p=0 was the best fit

#####
all <- list("simple_abs" = simple_abs, # total YFP, p=0 
            "simple_frac" = simple_frac,#, # fraction YFP, p=0
             "flow_frac" = flow_frac,  # Ki67 model with fixed input
             "simple_sim_abs" = simple_sim_abs, # total YFP, estimate thymus at same time
             "simple_sim_frac" = simple_sim_frac, # frac YFP, estimate thymus at same time
             "flow_sim_frac" = flow_sim_frac_pzero)

all_shapiros <- map(all, shapiros)
all_CIs <- map(all, CIs)
all_CIs <- map_depth(all_CIs, 2, ~data.frame(.))
all_CIs <- map(all_CIs, ~plyr::rbind.fill(cbind(.[["CD4"]], "receptor" = rep("CD4", nrow(.[["CD4"]]))),
                               cbind(.[["CD8"]], "receptor" = rep("CD8", nrow(.[["CD8"]]))))
)

all_CIs_csv <- map2(all_CIs, names(all_CIs), ~ {.x$model = .y;.x}) %>%
  map(~{.$tmp <- rep(c("upper", "lower","value"), 2);.}) %>%
  map(~pivot_longer(., 1:(ncol(.)-3), names_to = "parameter", values_to = "value") %>%
  pivot_wider( names_from = tmp, values_from = value)
  )

all_CIs_csv <- do.call(rbind, all_CIs_csv)
rownames(all_CIs_csv) <- NULL

write.csv(all_CIs_csv, file = "../fitting/parameters/parameter_and_CIs_final.csv")


# Make plot of residence time estimates 

mu_CIs <- map(all, CIs_mu)
mu_CIs <- map(mu_CIs, ~rbind(c(.[["CD4"]], "receptor" = "CD4"),
                               c(.[["CD8"]], "receptor" = "CD8"))
) %>% 
  map(~data.frame(.)) 

mu_CIs <- map2(mu_CIs, names(mu_CIs), ~ cbind(.x, "Method" = c(.y,.y)))
mu_CIs <- do.call(rbind, mu_CIs)

literature_fitting <- rbind(literature_decay, mu_CIs)




num.lit.estimates=7
for (i in 1:3){
literature_fitting[[i]] <- as.double(literature_fitting[[i]])
}

literature_fitting$panel <- c(rep("Literature", num.lit.estimates*2), rep("Fitting thymus and\nperiphery separately", 6), 
                              rep("Fitting thymus and\nperiphery simultaneously", 6))
literature_fitting$panel = factor(literature_fitting$panel,
                                  levels=c("Literature", "Fitting thymus and\nperiphery separately","Fitting thymus and\nperiphery simultaneously"))
literature_fitting$Method <- c(literature_fitting$Method[1:(num.lit.estimates*2)],
                               rep("YFP cell count", 2), rep("YFP cell frequency",2),  rep("YFP cell frequency",2), rep("YFP cell count", 2), rep("YFP cell frequency",2),  rep("YFP cell frequency",2)) 
literature_fitting$colors <- factor(c(rep("Literature", num.lit.estimates*2), rep("Fitting on YFP only (p=0)", 4), rep("Fitting on Ki67 and YFP", 2), rep("Fitting on YFP only (p=0)", 4), rep("Fitting on Ki67 and YFP", 2)))

literature_fitting$colors <- factor(literature_fitting$colors, levels = c("Literature", "Fitting on YFP only (p=0)", "Fitting on Ki67 and YFP"))


lit="EB5424"
blue="1F78B4"
green="63B546"

color.labels=c("Literature", "Fitting on YFP only (p=0)", "Fitting on Ki67 and YFP (estimating p)")

# get the sizes of the caps on errorbars right! ggplot scales them by the number of items
# in each group. I have hacked this fix:
literature_fitting <- literature_fitting %>%
  group_by(Method) %>%
  mutate(
    width = 0.2
  )

literature_fitting$width[literature_fitting$Method=="YFP cell count"] = 0.1
literature_fitting$width[literature_fitting$Method=="YFP cell frequency"] = 0.2

  
parameter_plot <- ggplot(literature_fitting, aes(x  = fct_inorder(Method)))+
  geom_point(aes(y = value, color = colors, shape = colors, size=colors), position=position_dodge(width=0.5))+
  scale_color_manual(name = "colors",
                     labels = color.labels,
                     values = c("Black", "orangered2","royalblue3")) + 
  scale_shape_manual(name = "colors",
                     labels = color.labels,
                     values = c(15,18,19)) +
  scale_size_manual(name = "colors",
                    labels = color.labels,
                    values = c(4,6,4)) +
  guides(color = guide_legend(title = NULL, override.aes = list(linetype = 0))) + 
  guides(shape = guide_legend(title = NULL)) + #, override.aes = list(linetype = 0))) + 
  guides(size = guide_legend(title = NULL)) + #, override.aes = list(linetype = 0))) +
  geom_errorbar(aes(ymin = lower.2.5., ymax = upper.97.5., color = colors, width=width), 
                linewidth = 1.2, position=position_dodge(width=0.5))+
  #facet_wrap(~receptor + panel, scales = "free")+ 
  facet_grid(receptor ~ panel, scales = "free")+ 
  labs(x = "", y = bquote("Mean residence time," ~ mu^-1 ~"(days)"), color = "")+ 
  # scale_color_brewer(palette = "Paired")+
  scale_x_discrete(guide = guide_axis(angle = 45)) + #, labels = label_wrap(10))+
  theme(panel.spacing.y = unit(3, "lines")) + 
  theme(panel.spacing.x = unit(1, "lines")) + 
  theme(plot.title = element_text(hjust = 0.5, family = "Times", size = 18),
        text = element_text(family = "Helvetica", size = 18),
        legend.position = "bottom",
        #legent.title=element_text(NULL),
        panel.grid= element_line(color = "grey",
                                  linewidth = 0.75),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 16),
        legend.key = element_rect(fill = "transparent"),
        legend.text=element_text(NULL),
        panel.background = element_rect(fill = "#FBFCFC"),
        strip.background = element_rect(color = "black", fill = "#FBFCFC"),
        strip.text = element_text(face = "bold", size = 18))+
  scale_y_continuous(trans = 'log10', labels=ggplot_scientific_notation_axes_labels, minor_breaks= unique(as.numeric(1:10 %o% 10 ^ (0:3)))) + 
  coord_cartesian(ylim=c(9,280)) +
  theme(plot.margin = margin(1,1,1,1.2, "cm")) # top right bottom lftt

  #scale_x_discrete(guide = guide_axis(angle = 90), labels=c("Mohri et al. (2001)", "Ribeiro et al. (2002)",
   #"Den Braber et al. (2012)", "Hogan et al. (2015)",
   #"Rane et al. (2018)", "Rane et al. (2022)",
  # rep(c("Cell counts", "Frequencies"), 4)))


plotwidth=12

pdf("../figures/fig_5.pdf", width=plotwidth, height=plotwidth/1.2)
print(parameter_plot)
dev.off()

pdf("../../Manuscript/PaperFigures/Figure-5.pdf", width=plotwidth, height=plotwidth/1.2)
print(parameter_plot)
dev.off()


































