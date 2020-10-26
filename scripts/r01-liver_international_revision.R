# Revision for Liver International

library(tidyverse)
library(here)
library(metafor)
library(meta)
library(glue)
library(openxlsx)
library(bernr)


# Setup functions for later
cm_to_inch <- 2.54

get_oddsratios <- function(metareg.object) {
  bind_cols(term = rownames(metareg.object$b),
            OR = metareg.object$b[,1], 
            ci.lb = metareg.object$ci.lb, 
            ci.ub = metareg.object$ci.ub) %>% 
    mutate_if(is.numeric, exp) %>% 
    bind_cols(pval = metareg.object$pval) %>% 
    mutate(pval = round(pval, 4))
}


# Read data
cirr_data <- read_rds(here("data", "02-cirrhosis.rds"))
fib_data <- read_rds(here("data", "02-fibrosis_apri_fib4.rds"))

# Modifications for nicer output
cirr_data <- cirr_data %>% 
  mutate(author = case_when(author == "Lemoine blood" ~ "Lemoine blood donors", 
                            author == "Lemoine comm" ~ "Lemoine community", 
                            TRUE ~ author)) %>% 
  mutate(cirr_test = case_when(cirr_test == "apri" ~ "APRI", 
                               cirr_test == "fibroscan" ~ "Transient elastography",
                               cirr_test == "fibrotest" ~ "Fibrotest", 
                               TRUE ~ as.character(cirr_test)), 
         cirr_test = factor(cirr_test, levels = c("Transient elastography", 
                                                  "APRI", 
                                                  "Fibrotest")))

# Main model without Aberra
m_cirrhosis <- metaprop(event = cirrhosis, n_fibrosis_assessed, 
                        data = cirr_data, comb.fixed = FALSE, 
                        studlab = glue("{author} ({year})"), 
                        method = "Inverse", 
                        subset = author != "Aberra")



m_cirrhosis %>% 
  forest(xlab = "Proportion with cirrhosis (%)", 
         pscale = 100, digits = 1, 
         leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
         leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
         smlab = "", 
         rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
         digits.addcols.left = 0, xlim = c(0, 30),
         digits.tau2 = 2, print.by = FALSE, col.by  = "black",
         overall = TRUE,
         overall.hetstat = FALSE)


update(m_cirrhosis, byvar = studypop2)

# Subgroup model without Aberra
cirr_bystudypop <- function() {
  update(m_cirrhosis, byvar = studypop2) %>% 
    forest(xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2, print.by = FALSE, col.by  = "black",
           overall = FALSE,
           overall.hetstat = FALSE, 
           sortvar = author)
}


# Write subgroup model
png("graph/r01-cirrhosis_studypop_sens_noaberra.png",
    res = 300, width = 24, height = 16, unit = "cm")
cirr_bystudypop()
dev.off()

pdf("graph/r01-cirrhosis_studypop_sens_noaberra.pdf",
    width = 24/cm_to_inch, height = 16/cm_to_inch)
cirr_bystudypop()
dev.off()


metareg(m_cirrhosis, studypop2 + group + cirr_test) %>% 
  get_oddsratios() %>% 
  filter(term != "intrcpt") %>% 
  mutate(across(OR:ci.ub, round, 2),
         pval = bernr::nice_p(pval, 2))







# Main model without Vinikoor 2018 and Kilonzo (numbers 21, 11 for cirrhosis) (TE only performed in less than 60% of population)
m_cirrhosis <- metaprop(event = cirrhosis, n_fibrosis_assessed, 
                        data = cirr_data, comb.fixed = FALSE, 
                        studlab = glue("{author} ({year})"), 
                        method = "Inverse", 
                        subset = !(study_nr %in% c(21, 11)))

m_cirrhosis %>% 
  forest(xlab = "Proportion with cirrhosis (%)", 
         pscale = 100, digits = 1, 
         leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
         leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
         smlab = "", 
         rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
         digits.addcols.left = 0, xlim = c(0, 30),
         digits.tau2 = 2, print.by = FALSE, col.by  = "black",
         overall = TRUE,
         overall.hetstat = FALSE, 
         sortvar = author)


update(m_cirrhosis, byvar = studypop2) %>% 
  forest(xlab = "Proportion with cirrhosis (%)", 
         pscale = 100, digits = 1, 
         leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
         leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
         smlab = "", 
         rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
         digits.addcols.left = 0, xlim = c(0, 30),
         digits.tau2 = 2, print.by = FALSE, col.by  = "black",
         overall = TRUE,
         overall.hetstat = FALSE, 
         sortvar = author)

# Subgroup model without Kilonzo and Vinikoor 2018
cirr_bystudypop <- function() {
  update(m_cirrhosis, byvar = studypop2) %>% 
    forest(xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2, print.by = FALSE, col.by  = "black",
           overall = FALSE,
           overall.hetstat = FALSE, 
           sortvar = author)
}


# Write subgroup model
png("graph/r01-cirrhosis_studypop_sens_without_suff_coverage.png",
    res = 300, width = 24, height = 16, unit = "cm")
cirr_bystudypop()
dev.off()

pdf("graph/r01-cirrhosis_studypop_sens_without_suff_coverage.pdf",
    width = 24/cm_to_inch, height = 16/cm_to_inch)
cirr_bystudypop()
dev.off()


metareg(m_cirrhosis, studypop2 + group + cirr_test) %>% 
  get_oddsratios() %>% 
  filter(term != "intrcpt") %>% 
  mutate(across(OR:ci.ub, round, 2),
         pval = bernr::nice_p(pval, 2))





# Setting analysis restricted to mono-studies
m_cirrhosis <- metaprop(event = cirrhosis, n_fibrosis_assessed, 
                        data = cirr_data, comb.fixed = FALSE, 
                        studlab = glue("{author} ({year})"), 
                        method = "Inverse", 
                        subset = group == "mono")



cirr_bystudypop_mono <- function() {
  update(m_cirrhosis, byvar = studypop2) %>% 
    forest(xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2, print.by = FALSE, col.by  = "black",
           overall = FALSE,
           overall.hetstat = FALSE, 
           sortvar = author)
}


# Write subgroup model
png("graph/r01-cirrhosis_studypop_sens_onlymono.png",
    res = 300, width = 24, height = 14, unit = "cm")
cirr_bystudypop_mono()
dev.off()

pdf("graph/r01-cirrhosis_studypop_sens_onlymono.pdf",
    width = 24/cm_to_inch, height = 14/cm_to_inch)
cirr_bystudypop_mono()
dev.off()
