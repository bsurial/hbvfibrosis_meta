  library(tidyverse)
  library(here)
  library(metafor)
  library(meta)
  library(glue)
  library(openxlsx)
  library(bernr)
  
  
  data <- read_rds(here("processed", "01-hbv_fibrosis.rds"))
  
  
  # To get odds ratios between the groups, we have to exponentiate the logit-transformed estimates
  
  get_oddsratios <- function(metareg.object) {
    bind_cols(term = rownames(metareg.object$b),
              OR = metareg.object$b[,1], 
              ci.lb = metareg.object$ci.lb, 
              ci.ub = metareg.object$ci.ub) %>% 
      mutate_if(is.numeric, exp) %>% 
      bind_cols(pval = metareg.object$pval) %>% 
      mutate(pval = round(pval, 4))
  }
  
  
  # Preparation -------------------------------------------------------------
  
  studytable <- data %>% 
    select(group, author, year, journal, country, study_type, population_type_simple, age_group, institution_simple, n_studypop, n_chb, n_coinfected) %>% 
    arrange(desc(group), year)
  studytable
  
  write_rds(studytable, here("tabs", "02-studytable.rds"))
  
  
  # Study summary for excel file
  df <- data %>%
    mutate(study = glue("{author}, {journal} {year}")) %>% 
    select(group, study_nr, study, country, setting, population_type_simple, study_type, institution_simple, setting_simple,
           n_studypop, n_chb, n_coinfected, median_age_hbvpatients, n_chb_female, n_sig_fibrosis_fibroscan,n_sig_fibrosis_apri, 
           n_sig_fibrosis_fib4, n_sig_fibrosis_fibrotest, n_cirrhosis_fibroscan, n_cirrhosis_apri, 
           n_cirrhosis_fib4, n_cirrhosis_fibrotest) %>% 
    distinct() %>% 
    arrange(study_nr)
  
  write.xlsx(df, 'tabs/02-studysummaries.xlsx')
  
  
  
  # Create dataset with fibrosis values, and use fibroscan over APRI over Fib-4 
  fib_apri_fib4 <- data %>% 
    gather(n_sig_fibrosis_fibroscan,n_sig_fibrosis_apri, 
           n_sig_fibrosis_fib4, n_sig_fibrosis_fibrotest,
           key = fib_test, value = fibrosis) %>% 
    mutate(fib_test = str_remove(fib_test, "n_sig_fibrosis_")) %>% 
    mutate(fib_test = factor(fib_test, levels = c(
      "fibroscan", "apri", "fib4", "fibrotest"
    ))) %>% 
    filter(!is.na(fibrosis)) %>% 
    group_by(group, author, year, country) %>% 
    arrange(group, author, year, country, fib_test) %>% 
    slice(1) %>% 
    ungroup()
  
  
  # Create dataset with fibrosis values, and use fibroscan over Fib-4 over APRI 
  fib_fib4_apri <- data %>% 
    gather(n_sig_fibrosis_fibroscan,n_sig_fibrosis_apri, 
           n_sig_fibrosis_fib4, n_sig_fibrosis_fibrotest,
           key = fib_test, value = fibrosis) %>% 
    mutate(fib_test = str_remove(fib_test, "n_sig_fibrosis_")) %>% 
    mutate(fib_test = factor(fib_test, levels = c(
      "fibroscan", "fib4", "apri", "fibrotest"
    ))) %>% 
    filter(!is.na(fibrosis)) %>% 
    group_by(group, author, year, country) %>% 
    arrange(group, author, year, country, fib_test) %>% 
    slice(1) %>% 
    ungroup()
    
  # Create dataset with cirrhosis values, and use fibroscan over APRI
  cirr_data <- data %>% 
    gather(n_cirrhosis_fibroscan, n_cirrhosis_apri, 
           n_cirrhosis_fib4, n_cirrhosis_fibrotest, 
           key = "cirr_test", value = "cirrhosis") %>% 
    mutate(cirr_test = str_remove(cirr_test, "n_cirrhosis_"),
           cirr_test = factor(cirr_test, levels = 
                                c("fibroscan", "apri", 
                                  "fibrotest", "fib4"))) %>% 
    group_by(group, author, country, year) %>% 
    arrange(group, author, country, year, cirr_test) %>% 
    filter(!is.na(cirrhosis)) %>% 
    slice(1) %>% 
    ungroup()
    
  
  write_rds(fib_apri_fib4, here("data", "02-fibrosis_apri_fib4.rds"))
  write_rds(cirr_data, here("data", "02-cirrhosis.rds"))
  
  # Overall fibrosis ------------------------------------------------------------
  
  m_fibrosis <- metaprop(event = fibrosis, n_fibrosis_assessed, 
                         data = fib_apri_fib4, comb.fixed = FALSE, 
                         studlab = glue("{author} ({year})"), 
                         method = "Inverse")
  
  summary(m_fibrosis)
  
  sink(file = "results/02-meta_fibrosis_overall.txt")
  m_fibrosis
  sink()
  
  fib_overall <- function() {
    forest(m_fibrosis, xlab = "Proportion with significant fibrosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 40),
           digits.tau2 = 2)
  }
  
  cm_to_inch <- 2.54
  
  png("graph/02-fibrosis_overall.png",
      res = 300, width = 30, height = 16, unit = "cm")
  fib_overall()
  dev.off()
  
  pdf("graph/02-fibrosis_overall.pdf",
      width = 30/cm_to_inch, height = 16/cm_to_inch)
  fib_overall()
  dev.off()
   
  
  # subgroup fibrosis -------------------------------------------------------
  
  
  metareg(m_fibrosis, fib_test)
  metareg(m_fibrosis, group)
  metareg(m_fibrosis, country)
  metareg(m_fibrosis, setting_simple)
  metareg(m_fibrosis, population_type_simple)
  metareg(m_fibrosis, institution_simple)
  metareg(m_fibrosis, prop_chb)
  metareg(m_fibrosis, prop_chb_female)
  metareg(m_fibrosis, median_age_hbvpatients)
  metareg(m_fibrosis, prop_hbe_ag_tested_positive)
  metareg(m_fibrosis, studypop2)
  
  
  
  
  # Institution seems significant
  # -----------------------------
  m_fib_byinst <- update(m_fibrosis, byvar = institution_simple)
  
  
  fib_byinst <- function(){
   m_fib_byinst %>% 
      forest(xlab = "Proportion with significant fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 40),
             digits.tau2 = 2, print.byvar = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  sink(file = here("results", "02-meta_fibrosis_byinst.txt"))
  m_fib_byinst
  sink()
  
  png("graph/02-fibrosis_byinstitution.png",
      res = 300, width = 30, height = 22, unit = "cm")
  fib_byinst()
  dev.off()
  
  pdf("graph/02-fibrosis_byinstitution.pdf",
      width = 30/cm_to_inch, height = 22/cm_to_inch)
  fib_byinst()
  dev.off()
  
  
  # Studypop2 catpures better whether its teartiary or primary care / screening
  # ----------------------------------------------------------------------------
  m_fib_bycaretype <- update(m_fibrosis, byvar = studypop2)
  
  
  fib_bycaretype <- function(){
    m_fib_bycaretype %>% 
      forest(xlab = "Proportion with significant fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 40),
             digits.tau2 = 2, print.byvar = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  sink(file = here("results", "02-meta_fibrosis_bycaretype.txt"))
  m_fib_bycaretype
  sink()
  
  png("graph/02-fibrosis_bycaretype.png",
      res = 300, width = 30, height = 22, unit = "cm")
  fib_bycaretype()
  dev.off()
  
  pdf("graph/02-fibrosis_bycaretype.pdf",
      width = 30/cm_to_inch, height = 22/cm_to_inch)
  fib_bycaretype()
  dev.off()
  
  # by Infection group
  # ------------------
  
  m_fib_bygroup <- update(m_fibrosis, byvar = if_else(group == "mono", 
                                     "HBV infection", 
                                     "HIV/HBV co-infection")) 
  
  summary(m_fib_bygroup)
  
  sink(file = here("results", "02-meta_fibrosis_bygroup.txt"))
  m_fib_bygroup
  sink()
  
  fib_bygroup <- function() {
    m_fib_bygroup %>% 
      forest(xlab = "Proportion with significant fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 40),
             digits.tau2 = 2, print.byvar = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
    
  
  
  
  
  
  png("graph/02-fibrosis_bygroup.png",
      res = 300, width = 28, height = 20, unit = "cm")
  fib_bygroup()
  dev.off()
  
  pdf("graph/02-fibrosis_bygroup.pdf",
      width = 28/cm_to_inch, height = 20/cm_to_inch)
  fib_bygroup()
  dev.off()
  
  
  # by test
  #---------
  m_fib_bytest <- update(m_fibrosis, byvar = fib_test) 
  
  sink(file = here("results", "02-meta_fibrosis_bytest.txt"))
  m_fib_bytest
  sink()
  
  
  fib_bytest <- function(){
    m_fib_bytest %>% 
      forest(xlab = "Proportion with significant fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 40),
             digits.tau2 = 2, print.byvar = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  
  png("graph/02-fibrosis_bytest.png",
      res = 300, width = 28, height = 23, unit = "cm")
  fib_bytest()
  dev.off()
  
  pdf("graph/02-fibrosis_bytest.pdf",
      width = 28/cm_to_inch, height = 23/cm_to_inch)
  fib_bytest()
  dev.off()
  
  
  # Only fibroscan
  # --------------
  
  m_fib_fibroscan <- update(m_fibrosis, subset = fib_test == "fibroscan")
  
  fib_fibroscan <- function() {
    m_fib_fibroscan %>% 
      forest(xlab = "Proportion with significant fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nsign. fibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 40),
             digits.tau2 = 2, print.byvar = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  sink(file = "results/02-meta_fibrosis_fibroscan.txt")
  m_fib_fibroscan
  sink()
  
  
  png("graph/02-fibrosis_fibroscan.png",
      res = 300, width = 28, height = 11, unit = "cm")
  fib_fibroscan()
  dev.off()
  
  pdf("graph/02-fibrosis_fibroscan.pdf",
      width = 28/cm_to_inch, height = 11/cm_to_inch)
  fib_fibroscan()
  dev.off()
  
  
  # Cirrhosis overall -------------------------------------------------------
  
  m_cirrhosis <- metaprop(event = cirrhosis, n_fibrosis_assessed, 
                          data = cirr_data, comb.fixed = FALSE, 
                          studlab = glue("{author} ({year})"), 
                          method = "Inverse")
  
  summary(m_cirrhosis)
  
  sink(file = "results/02-meta_cirrhosis_overall.txt")
  m_cirrhosis
  sink()
  
  cirr_overall <- function(){
    forest(m_cirrhosis, xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2)
  }
  
  
  png("graph/02-cirrhosis_overall.png",
      res = 300, width = 26, height = 14, unit = "cm")
  cirr_overall()
  dev.off()
  
  pdf("graph/02-cirrhosis_overall.pdf",
      width = 26/cm_to_inch, height = 14/cm_to_inch)
  cirr_overall()
  dev.off()
  
  
  # subgroup cirrhosis  -------------------------------------------------------------------

  
  metareg(m_cirrhosis, group) %>% get_oddsratios()
  metareg(m_cirrhosis, cirr_test) %>% get_oddsratios()
  metareg(m_cirrhosis, studypop2) %>% get_oddsratios()
  
  
  metareg(m_cirrhosis, studypop2 + group + cirr_test) %>% 
    get_oddsratios()
  
  metareg(m_cirrhosis, prop_chb)
  metareg(m_cirrhosis, prop_chb_female) 
  metareg(m_cirrhosis, prop_chb_female) %>% bubble(xlab = "Proportion of women", ylab = "Proportion of cirrhosis (logit transformed)")
  metareg(m_cirrhosis, median_age_hbvpatients) %>% bubble(xlab = "Median age (years)", ylab = "Proportion of cirrhosis (logit transformed)")
  x <- metareg(m_cirrhosis, prop_hbe_ag_tested_positive)
  
  
  png("graph/02-bubble_cirrhosis_hbe.png",
      res = 300, width = 20, height = 14, unit = "cm")
  bubble(x, xlab = "Proportion of positive HBeAg", ylab = "Proportion of cirrhosis (logit transformed)",
         studlab = TRUE, pos.studlab = 4, offset = 1.8, xlim = c(0, 0.26))
  dev.off()
  
  pdf("graph/02-bubble_cirrhosis_hbe.pdf",
      width = 20/cm_to_inch, height = 14/cm_to_inch)
  bubble(x, xlab = "Proportion of positive HBeAg", ylab = "Proportion of cirrhosis (logit transformed)",
         studlab = TRUE, pos.studlab = 4, offset = 1.8, xlim = c(0, 0.26))
  dev.off()
  

  # by institution
  #---------------
  
  m_cirrhosis_byinst <- update(m_cirrhosis, byvar = institution_simple)
  
  sink(file = "results/02-meta_cirrhosis_byinst.txt")
  m_cirrhosis_byinst
  sink()
  
  cirr_byinst <- function() {
    m_cirrhosis_byinst %>% 
    forest(xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2, print.by = FALSE, col.by  = "black",
           overall = FALSE,
           overall.hetstat = FALSE)
  }
  
  png("graph/02-cirrhosis_byinstitution.png",
      res = 300, width = 28, height = 20, unit = "cm")
  cirr_byinst()
  dev.off()
  
  pdf("graph/02-cirrhosis_byinstitution.pdf",
      width = 28/cm_to_inch, height = 20/cm_to_inch)
  cirr_byinst()
  dev.off()
  
  
  # by care type (studypop2 = tertiary vs. primary care)
  #-----------------------------------------------------
  
  m_cirrhosis_bycaretype <- update(m_cirrhosis, byvar = studypop2)
  
  sink(file = "results/02-meta_cirrhosis_bycaretype.txt")
  m_cirrhosis_bycaretype
  sink()
  
  cirr_bycaretype <- function() {
    m_cirrhosis_bycaretype %>% 
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
  
  png("graph/02-cirrhosis_bycaretype.png",
      res = 300, width = 28, height = 20, unit = "cm")
  cirr_bycaretype()
  dev.off()
  
  pdf("graph/02-cirrhosis_bycaretype.pdf",
      width = 26/cm_to_inch, height = 17/cm_to_inch)
  cirr_bycaretype()
  dev.off()
  
  
  
  
  
  # by infection group
  #-------------------
  
  m_cirr_bygroup <- update(m_cirrhosis, byvar = if_else(group == "mono", 
                                                        "HBV infection", 
                                                        "HIV/HBV co-infection"))
  
  sink(file = "results/02-meta_cirrhosis_bygroup.txt")
  m_cirr_bygroup
  sink()
  
  
  cirr_bygroup <- function() {
    m_cirr_bygroup %>% 
      forest(xlab = "Proportion with cirrhosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black", 
             sortvar = author)
  }
  
  png("graph/02-cirrhosis_bygroup.png",
      res = 300, width = 28, height = 18, unit = "cm")
  cirr_bygroup()
  dev.off()
  
  pdf("graph/02-cirrhosis_bygroup.pdf",
      width = 28/cm_to_inch, height = 18/cm_to_inch)
  cirr_bygroup()
  dev.off()
  
  
  
  
  # by cirrhosis test
  #-------------------
  
  m_cirr_bytest <- update(m_cirrhosis, byvar = cirr_test)
  
  sink(file = "results/02-meta_cirrhosis_bytest.txt")
  m_cirr_bytest
  sink()
  
  
  cirr_bytest <- function() {
    m_cirr_bytest %>% 
      forest(xlab = "Proportion with cirrhosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black", 
             sortvar = author)
  }
  
  png("graph/02-cirrhosis_bytest.png",
      res = 300, width = 28, height = 19, unit = "cm")
  cirr_bytest()
  dev.off()
  
  pdf("graph/02-cirrhosis_bytest.pdf",
      width = 28/cm_to_inch, height = 19/cm_to_inch)
  cirr_bytest()
  dev.off()
  
  
  
  # Cirrhosis by fibroscan
  # ----------------------
  
  m_cirr_fibroscan <- update(m_cirrhosis, subset = cirr_test == "fibroscan")
  
  cirr_fibroscan <- function() {
  m_cirr_fibroscan %>% 
    forest(xlab = "Proportion with cirrhosis (%)", 
           pscale = 100, digits = 1, 
           leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
           leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
           smlab = "", 
           rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
           digits.addcols.left = 0, xlim = c(0, 30),
           digits.tau2 = 2, print.by = FALSE, col.by  = "black",
           overall = FALSE,
           overall.hetstat = FALSE)
  }
  
  sink(file = "results/02-meta_cirrhosis_fibroscan.txt")
  m_cirr_fibroscan
  sink()
  
  
  png("graph/02-cirrhosis_fibroscan.png",
      res = 300, width = 28, height = 11, unit = "cm")
  cirr_fibroscan()
  dev.off()
  
  pdf("graph/02-cirrhosis_fibroscan.pdf",
      width = 28/cm_to_inch, height = 11/cm_to_inch)
  cirr_fibroscan()
  dev.off()
  
  
  # Analysis with FIBROSCAN ONLY
  # ----------------------------
  
  # cirrhosis by institution, fibroscan only
  #-----
  
  m_cirr_inst_fibroscan <- update(m_cirrhosis, subset = cirr_test == "fibroscan",
                             byvar = institution_simple )
  
  cirr_inst_fibroscan <- function() {
    m_cirr_inst_fibroscan %>% 
      forest(xlab = "Proportion with cirrhosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  pdf("graph/02-cirrhosis_institution_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 16/cm_to_inch)
  cirr_inst_fibroscan()
  dev.off()
  
  png("graph/02-cirrhosis_institution_fibroscan_only.png",
      width = 26, height = 16, unit = "cm", res = 300)
  cirr_inst_fibroscan()
  dev.off()
  
  
  # cirrhosis by care type, fibroscan only
  #-----
  
  m_cirr_caretype_fibroscan <- update(m_cirrhosis, subset = cirr_test == "fibroscan",
                                      byvar = studypop2 )
  
  cirr_caretype_fibroscan <- function() {
    m_cirr_caretype_fibroscan %>% 
      forest(xlab = "Proportion with cirrhosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  pdf("graph/02-cirrhosis_caretype_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 16/cm_to_inch)
  cirr_caretype_fibroscan()
  dev.off()
  
  png("graph/02-cirrhosis_caretype_fibroscan_only.png",
      width = 26, height = 16, unit = "cm", res = 300)
  cirr_caretype_fibroscan()
  dev.off()
  
  
  
  # Cirrhosis by infection group, fibroscan only
  #-----
  
  m_cirr_group_fibroscan <- update(m_cirrhosis, subset = cirr_test == "fibroscan",
                                  byvar = group )
  
  cirr_group_fibroscan <- function() {
    m_cirr_group_fibroscan %>% 
      forest(xlab = "Proportion with cirrhosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "cirrhosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\ncirrhosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  
  pdf("graph/02-cirrhosis_group_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 14/cm_to_inch)
  cirr_group_fibroscan()
  dev.off()
  
  png("graph/02-cirrhosis_group_fibroscan_only.png",
      width = 26, height = 14, unit = "cm", res = 300)
  cirr_group_fibroscan()
  dev.off()
  
  
  # Fibrosis by institution, fibroscan only
  #-----
  
  m_fib_inst_fibroscan <- update(m_fibrosis, subset = fib_test == "fibroscan",
                                  byvar = institution_simple )
  
  fib_inst_fibroscan <- function() {
    m_fib_inst_fibroscan %>% 
      forest(xlab = "Proportion with fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nfibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  pdf("graph/02-fibrosis_institution_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 16/cm_to_inch)
  fib_inst_fibroscan()
  dev.off()
  
  png("graph/02-fibrosis_institution_fibroscan_only.png",
      width = 26, height = 16, unit = "cm", res = 300)
  cirr_inst_fibroscan()
  dev.off()
  
  
  # Fibrosis by caretype, fibroscan only
  #-----
  
  m_fib_caretype_fibroscan <- update(m_fibrosis, subset = fib_test == "fibroscan",
                                     byvar = studypop2)
  
  fib_caretype_fibroscan <- function() {
    m_fib_caretype_fibroscan %>% 
      forest(xlab = "Proportion with fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nfibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  pdf("graph/02-fibrosis_caretype_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 16/cm_to_inch)
  fib_caretype_fibroscan()
  dev.off()
  
  png("graph/02-fibrosis_caretype_fibroscan_only.png",
      width = 26, height = 16, unit = "cm", res = 300)
  fib_caretype_fibroscan()
  dev.off()
  
  
  
  # Fibrosis by infection group, fibroscan only
  #-----
  
  m_fib_inst_fibroscan <- update(m_fibrosis, subset = fib_test == "fibroscan",
                                 byvar = group )
  
  fib_group_fibroscan <- function() {
    m_fib_inst_fibroscan %>% 
      forest(xlab = "Proportion with fibrosis (%)", 
             pscale = 100, digits = 1, 
             leftcols = c("studlab", "country", "fibrosis", "n_fibrosis_assessed"),
             leftlabs = c("Author (year)", "Country", "N with\nfibrosis", "N study\npopulation"),
             smlab = "", 
             rightlabs = c("Proportion (%)", "[95% CI]", "Weight"),
             digits.addcols.left = 0, xlim = c(0, 30),
             digits.tau2 = 2, print.by = FALSE, col.by  = "black",
             overall = FALSE,
             overall.hetstat = FALSE)
  }
  
  pdf("graph/02-fibrosis_group_fibroscan_only.pdf",
      width = 26/cm_to_inch, height = 14/cm_to_inch)
  fib_group_fibroscan()
  dev.off()
  
  png("graph/02-fibrosis_group_fibroscan_only.png",
      width = 26, height = 14, unit = "cm", res = 300)
  fib_group_fibroscan()
  dev.off()
  
  
  funnel(m_cirrhosis)
  dmetar::eggers.test(m_cirrhosis)
  
  
  
  
  
  
  # Table for meta regression
  
  sex <- metareg(m_cirrhosis, prop_chb_female) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt")
  
  age <- metareg(m_cirrhosis, median_age_hbvpatients) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt")
  
  test <- metareg(m_cirrhosis, cirr_test) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt")
  
  pop2 <- metareg(m_cirrhosis, studypop2) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt")
  
  groups <- metareg(m_cirrhosis, group) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt")
  
  uni <- bind_rows(sex, age, test, pop2, groups) %>% 
    mutate(pval = nice_p(pval, 2)) %>% 
    mutate(across(where(is.numeric), comma, digits = 2, trailing = TRUE)) %>% 
    mutate(OR = glue("{OR} ({ci.lb} to {ci.ub})")) %>% 
    select(term, "univariable OR" = OR, "p-value" = pval)
  
  
  metareg(m_cirrhosis, cirr_test)
  
  multi <- metareg(m_cirrhosis, studypop2 + group + cirr_test) %>% 
    get_oddsratios() %>% 
    filter(term != "intrcpt") %>% 
    mutate(pval = nice_p(pval, 2)) %>% 
    mutate(across(where(is.numeric), comma, digits = 2, trailing = TRUE)) %>% 
    mutate(OR = glue("{OR} ({ci.lb} to {ci.ub})")) %>% 
    select(term, "adjusted OR" = OR, "p-value" = pval)
  
  
  reg_table <- uni %>% 
    left_join(multi, by = "term")
  
  ad <- tibble(term = c("fibroscan", "Tertiary Care", "coinfected"), "univariable OR" = 1, "p-value.x" = NA, "adjusted OR" = 1, "p-value.y" = NA) %>% 
    mutate(across(where(is.numeric), as.character))
  
  reg_table <- reg_table %>% 
    add_row(ad)
  
  # # Uncomment to update
  # WriteXLS::WriteXLS(reg_table, here("results", "02-metaregression_table.xlsx"))
