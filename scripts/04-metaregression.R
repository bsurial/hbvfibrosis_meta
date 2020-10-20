# Run 02-meta.R first

metareg(m_fibrosis, prop_chb_female) %>% 
  bubble(xlab = "Proportion of female patients", 
         ylab = "Proportion of patients with sig. fibrosis (logit)",
         main = "Metaregression: Sig. fibrosis ~ female sex")


metareg(m_fibrosis, median_age_hbvpatients) %>% 
  bubble(xlab = "Median age of study populations", 
         ylab = "Proportion of patients with sig. fibrosis (logit)",
         main = "Metaregression: Sig. fibrosis ~ Median age")

metareg(m_fibrosis, prop_hbe_ag_tested_positive) %>% 
  bubble(xlab = "Proportion of HBeAg positive", 
         ylab = "Proportion of patients with sig. fibrosis (logit)",
         main = "Metaregression: Sig. Fibrosis ~ HBeAg")



metareg(m_cirrhosis, prop_chb_female) %>% 
  bubble(xlab = "Proportion of female patients", 
         ylab = "Proportion of patients with cirrhosis (logit)",
         main = "Metaregression: Cirrhosis ~ female sex")

metareg(m_cirrhosis, median_age_hbvpatients) %>% 
  bubble(xlab = "Median age of study populations", 
         ylab = "Proportion of patients with cirrhosis (logit)",
         main = "Metaregression: Cirrhosis ~ Median age")

metareg(m_cirrhosis, prop_hbe_ag_tested_positive) %>% 
  bubble(xlab = "Proportion of HBeAg positive", 
         ylab = "Proportion of patients with cirrhosis (logit)",
         main = "Metaregression: Cirrhosis ~ HBeAg")

