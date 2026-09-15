#!/usr/bin/env Rscript
# ============================================================================
# R2.5: directional alignment of cis-pQTL effects with therapeutic mechanism
#       for all 23 pQTL-supported, launched target-indication pairs
# ============================================================================
# Nothing is excluded. Every pair is shown with the reasoning exposed so a
# reader can audit each call.
#
# Expected sign is derived a priori from three inputs, none of which depend on
# the observed beta:
#
#   1. analyte_class - what the plasma measurement represents relative to the
#      therapeutic axis:
#        "effector"      secreted bioactive protein; higher protein = more
#                        pathway activity (PCSK9, F2, VWF, SERPINA1, CSF3,
#                        IL12B, APOB)
#        "shed_receptor" membrane receptor; the plasma species is the shed
#                        ectodomain, which acts as a ligand decoy, so higher
#                        protein = LESS pathway activity (CSF3R, IL4R, FLT3,
#                        KIT)
#        "enzyme"        intracellular enzyme; higher protein = more activity
#                        (EGLN1)
#   2. drug_direction - whether the therapy raises or lowers pathway activity
#   3. trait_orientation - whether the matched GWAS trait runs in the same
#      direction as the indication ("direct") or the opposite ("inverse", e.g.
#      neutrophil count for neutropenia)
#
# Derivation:
#   higher protein -> more pathway activity?   effector/enzyme = yes, shed = no
#   more activity  -> more disease?            inhibitor = yes, agonist = no
#   => expected sign on a disease-scale trait, then flipped if the matched
#      trait is inversely oriented.
# ============================================================================

suppressPackageStartupMessages({library(tidyverse); library(openxlsx)})

parts <- readRDS("/tmp/r25_parts.rds")
lau <- parts$lau; cis <- parts$cis; trans <- parts$trans

ann <- tribble(
 ~gene,~indication_mesh_term,~analyte_class,~drug_direction,~modality,~example_drug,~trait_orientation,
 "APOB","Hyperlipoproteinemia Type II","effector","lowers","ASO knockdown","Mipomersen","direct",
 "CSF3","Neutropenia","effector","raises","Recombinant protein","Filgrastim","inverse",
 "CSF3R","Leukopenia","shed_receptor","raises","Receptor agonist (ligand)","Filgrastim/pegfilgrastim","inverse",
 "CSF3R","Neutropenia","shed_receptor","raises","Receptor agonist (ligand)","Filgrastim/pegfilgrastim","inverse",
 "EGLN1","Anemia","enzyme","lowers","Small molecule","Roxadustat, daprodustat","inverse",
 "F2","Thrombosis","effector","lowers","Small molecule","Dabigatran","direct",
 "F2","Venous Thrombosis","effector","lowers","Small molecule","Dabigatran","direct",
 "FLT3","Thrombocytopenia","shed_receptor","lowers","Small molecule","Fedratinib, gilteritinib, midostaurin","inverse",
 "FLT3","Thrombocytosis","shed_receptor","lowers","Small molecule","Midostaurin, quizartinib","direct",
 "IL12B","Colitis, Ulcerative","effector","lowers","mAb","Ustekinumab","direct",
 "IL12B","Crohn Disease","effector","lowers","mAb","Ustekinumab","direct",
 "IL12B","Psoriasis","effector","lowers","mAb","Ustekinumab, briakinumab","direct",
 "IL12B","Arthritis, Psoriatic","effector","lowers","mAb","Ustekinumab","direct",
 "IL4R","Asthma","shed_receptor","lowers","mAb","Dupilumab","direct",
 "IL4R","Eosinophilic Esophagitis","shed_receptor","lowers","mAb","Dupilumab","direct",
 "KIT","Thrombocytopenia","shed_receptor","lowers","Small molecule","Imatinib, dasatinib, avapritinib","inverse",
 "PCSK9","Hypercholesterolemia","effector","lowers","mAb/siRNA","Evolocumab, inclisiran","direct",
 "PCSK9","Hyperlipoproteinemia Type II","effector","lowers","mAb","Alirocumab, evolocumab","direct",
 "PCSK9","Hyperlipidemias","effector","lowers","mAb","Alirocumab, evolocumab","direct",
 "PCSK9","Atherosclerosis","effector","lowers","mAb","Evolocumab","direct",
 "SERPINA1","alpha 1-Antitrypsin Deficiency","effector","raises","Augmentation","Alpha-1 proteinase inhibitor","direct",
 "VWF","Hemophilia A","effector","raises","Replacement","VWF/FVIII concentrate, desmopressin","direct",
 "VWF","von Willebrand Diseases","effector","raises","Replacement","Vonicog alfa","direct")

notes <- tribble(~gene,~indication_mesh_term,~notes,
 "APOB","Hyperlipoproteinemia Type II","PAV present: cis index SNPs rs1367117 (missense, CADD 22.1) and rs17240441 (inframe deletion, CADD 16.5), both in APOB, r2=0.93 (TOP-LD EUR, r2>=0.6, MAF>=0.01)",
 "CSF3","Neutropenia","cis/trans opposite",
 "CSF3R","Leukopenia","composite WBC trait; cis/trans opposite",
 "CSF3R","Neutropenia","cis/trans opposite (17_39871710_C_G)",
 "EGLN1","Anemia","trans-association only",
 "F2","Thrombosis","",
 "F2","Venous Thrombosis","",
 "FLT3","Thrombocytopenia","approved indications are myeloid leukaemia/myelofibrosis/mastocytosis; thrombocytopenia label follows MeSH mapping; GWAS trait plateletcrit",
 "FLT3","Thrombocytosis","approved indications are myeloid malignancies; thrombocytosis label follows MeSH mapping; GWAS trait plateletcrit",
 "IL12B","Colitis, Ulcerative","13 concordant cis associations across 6 traits",
 "IL12B","Crohn Disease","13 concordant cis associations across 6 traits",
 "IL12B","Psoriasis","single cis association; discrepancy vs concordant IBD associations unexplained; plasma p40 may not reflect skin tissue levels",
 "IL12B","Arthritis, Psoriatic","single cis association; skin trait for MSK indication; discrepancy vs concordant IBD associations unexplained",
 "IL4R","Asthma","shed receptor decoy",
 "IL4R","Eosinophilic Esophagitis","shed receptor decoy",
 "KIT","Thrombocytopenia","approved indications are CML/mastocytosis/myeloma; thrombocytopenia label follows MeSH mapping; GWAS trait plateletcrit; PAV absent",
 "PCSK9","Hypercholesterolemia","",
 "PCSK9","Hyperlipoproteinemia Type II","",
 "PCSK9","Hyperlipidemias","",
 "PCSK9","Atherosclerosis","",
 "SERPINA1","alpha 1-Antitrypsin Deficiency","",
 "VWF","Hemophilia A","",
 "VWF","von Willebrand Diseases","direction indication-specific")

tbl <- lau %>%
  left_join(cis, by="ti_uid") %>% left_join(trans, by="ti_uid") %>%
  left_join(ann, by=c("gene","indication_mesh_term")) %>%
  left_join(notes, by=c("gene","indication_mesh_term")) %>%
  mutate(
    mr_evidence_used = if_else(!is.na(cis_sign), "cis", "trans only"),
    beta_sign = coalesce(cis_sign, trans_sign),
    beta      = coalesce(cis_beta, trans_beta),
    n_assoc   = coalesce(n_cis, n_trans),
    matched_gwas_trait = cis_traits,
    # a priori derivation
    higher_protein_raises_activity = analyte_class %in% c("effector","enzyme"),
    more_activity_raises_disease   = drug_direction == "lowers",
    exp_disease_scale = if_else(higher_protein_raises_activity == more_activity_raises_disease, "+", "-"),
    expected_sign = if_else(trait_orientation == "inverse",
                            if_else(exp_disease_scale=="+","-","+"), exp_disease_scale),
    alignment = case_when(is.na(beta_sign) ~ "not assessable",
                          beta_sign=="mixed" ~ "mixed sign",
                          beta_sign==expected_sign ~ "aligned",
                          TRUE ~ "not aligned")) %>%
  arrange(desc(alignment=="aligned"), gene, indication_mesh_term)

out <- tbl %>% transmute(
  Target=gene, MeSH=indication_mesh_id, Indication=indication_mesh_term,
  Area=therapeutic_area, Pair=ti_uid,
  Direction=if_else(drug_direction=="lowers","Inhibitor","Agonist"),
  Example=example_drug,
  analyte_class, matched_gwas_trait, trait_orientation,
  mr_evidence=mr_evidence_used, n_assoc, mr_beta=beta, mr_beta_sign=beta_sign,
  expected_sign, alignment, trans_beta, notes)

write_tsv(out, "output/r2_5_alignment_table.tsv", na="")
wb <- createWorkbook(); addWorksheet(wb,"R2.5 directional alignment")
writeData(wb,"R2.5 directional alignment", out)
freezePane(wb,"R2.5 directional alignment", firstActiveRow=2, firstActiveCol=1)
setColWidths(wb,"R2.5 directional alignment", cols=1:ncol(out), widths="auto")
saveWorkbook(wb,"output/r2_5_alignment_table.xlsx", overwrite=TRUE)

cat("\n=== ALIGNMENT: all 23 launched pQTL-supported pairs ===\n\n")
print(as.data.frame(out %>% transmute(Target, Indication=str_trunc(Indication,24), Direction,
  class=str_replace(analyte_class,"_"," "), orient=trait_orientation,
  ev=mr_evidence, n=n_assoc, beta=mr_beta, exp=expected_sign, alignment,
  notes=str_trunc(notes,34))), row.names=FALSE, right=FALSE)
cat("\n=== TALLY ===\n"); print(as.data.frame(out %>% count(alignment)), row.names=FALSE)
cat(sprintf("\naligned: %d / %d (%.0f%%)\n", sum(out$alignment=="aligned"), nrow(out),
            100*mean(out$alignment=="aligned")))
cat("\n-> output/r2_5_alignment_table.tsv\n-> output/r2_5_alignment_table.xlsx\n\n")
