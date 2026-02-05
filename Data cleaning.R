# Load required libraries
library(dplyr)
library(lubridate)
library(stringr)
library(readr)
library(data.table)
InpatDiagCodes <- read_csv("InpatDiagCodes.csv")
InpatEpisodes <- read_csv("InpatEpisodes.csv")
Demographic <- read_csv("Demographics.csv")
CancerDiagCodes <- read_csv("/home/shared/IORD_CancerAMRSurveillance_56_20250430/CancerDiagCodes.csv")
antibiotic <- read_csv('/home/shared/IORD_CancerAMRSurveillance_56_20250430/AntibioticsEPR.csv') 
######All patients
# Merge in patient episodes with demographics
df_episode <- InpatEpisodes %>%
  left_join(Demographic, by = "ClusterID")

#  Convert date columns to Date type
df_episode <- df_episode %>%
  mutate(
    AdmissionDate = as.POSIXct(AdmissionDate), # Convert both date and time with an associated time zone
    DischargeDate = as.POSIXct(DischargeDate),
    LinkedDeathdate = as.POSIXct(LinkedDeathdate, origin="1970-01-01", format="%Y-%m-%d"), #This version is typically used when:
    #LinkedDeathdate is stored as numeric/charactor values — e.g., number of seconds or days since 1970-01-01 (UNIX epoch)
    #Tell R how to interpret those numbers/charactors as actual dates
    LinkedBirthMonth = as.POSIXct(LinkedBirthMonth)
  )

# Filter out invalid or missing admission/discharge dates
df_episode <- df_episode %>%
  filter(
    !is.na(AdmissionDate),
    !is.na(DischargeDate),
    DischargeDate > AdmissionDate
  )

# Remove records where death occurred before admission, and keep rows where there's no recorded death date
df_episode <- df_episode %>%
  filter(is.na(LinkedDeathdate) | LinkedDeathdate >= AdmissionDate)

# Remove records with missing sex. Because: A person’s sex at birth doesn’t typically change across records.
# So, multiple recorded sexes for a single patient (e.g., both "M" and "F") is usually a data quality issue, not reality.
# If a patient has two different values, we can't tell which one is true — so keeping the "most frequent" (mode) could propagate errors.
# Also by using mode assumes the majority value is correct
df_episode <- df_episode %>%
  filter(!is.na(LinkedSex))

# Keep only patients with consistent sex across all episodes
valid_ids <- df_episode %>%
  group_by(ClusterID) %>%
  filter(n_distinct(LinkedSex) == 1) %>%
  pull(ClusterID) %>%
  unique()

df_episode <- df_episode %>%
  filter(ClusterID %in% valid_ids)

# Recode EthnicGroupCode into simplified 'ethnicity'
df_episode <- df_episode %>%
  mutate(
    EthnicGroupCode = str_to_upper(str_trim(as.character(EthnicGroupCode))),
    ethnicity = case_when(
      EthnicGroupCode %in% c("A", "B", "C") ~ "White",
      EthnicGroupCode %in% c("99", "Z")     ~ NA_character_,
      TRUE                                  ~ "Non-white"
    )
  )

# Clean the ethnicity by replacing with the mode of ethnicity
### Assign modal ethnicity per patient (ClusterID)
find_mode <- function(x) {
  x <- na.omit(x) # exclude NA in ethnicity for each patient
  tab <- table(x)
  if (length(tab) == 0) {
    return(NA_character_) # if the patient only has NA ethnicity then retyrn NA in ethnicity
  }
  if (sum(tab == max(tab)) == 1) {  # If the patient (ClusterID)'s recorded ethnicities have only one maximum, then the mode is just names(which.max(tab))
    return(names(which.max(tab)))
  }
  return(NA_character_)  # return NA if there's a tie
}

df_episode <- df_episode %>%
  group_by(ClusterID) %>%
  mutate(ethnicity = find_mode(ethnicity)) %>% # This will apply find_mode() individually to the ethnicity values within each ClusterID
  ungroup()

# Remove overlapping admissions

# Drop duplicate spell records based on ID + dates to avoid unuseful information (e.g. wrong/random information)
df <- df_episode %>%
  select(ClusterID, SpellID, AdmissionDate, DischargeDate) %>%
  distinct()

# Convert to data.table for efficient shift operations
setDT(df)
setorder(df, ClusterID, AdmissionDate)

# Previous discharge date per patient
df[, Discharge_prev := shift(DischargeDate, type = "lag"), by = ClusterID]

# Overlap if admission starts before previous discharge
df[, overlap := as.integer(AdmissionDate < Discharge_prev)]

# Look ahead to see if this admission overlaps with the next one. If without fill=0, the last row in each ClusterID, where there is no next row, overlap will return NA
df[, overlap2 := shift(overlap, type = "lead", fill = 0), by = ClusterID]

# Mark overlapping SpellIDs
remove_IDs <- df[overlap == 1 | overlap2 == 1, unique(SpellID)]


# Remove overlapping admissions from original data. Remove both overlapping in case of any ambiguity (unsure which admission leads to the outcome)
df_episode <- df_episode %>%
  filter(!SpellID %in% remove_IDs)

max(df_episode$AdmissionDate, na.rm = TRUE)
# The data is to 24/04/2025, we only use data to 01/04/2025 as microbiology data is up to this date
df_episode <- df_episode %>% filter(AdmissionDate < as.POSIXct("2025-04-01")) 

######### Identify cancer episodes
cancer_diagnosis = read_csv('/home/shared/IORD_CancerAMRSurveillance_56_20250430/CancerDiagCodes.csv')
cancer_diagnosis <- cancer_diagnosis %>% 
  mutate(cancer_type_sub= case_when(
    grepl(paste(c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09",
                  "C10", "C11", "C12", "C13", "C14"), collapse = "|"),DiagCode) ~ "Lip/Oral/Pharynx",
    grepl(paste(c("C15", "C16", "C17", "C18", "C19",
                  "C20", "C21", "C22", "C23", "C24", "C25", "C26"), collapse = "|"),DiagCode) ~ "Digestive organs",
    grepl(paste(c("C30", "C31", "C32", "C33", "C34", "C35", "C36", "C37", "C38", "C39"), collapse = "|"),DiagCode) ~ "Respiratory organs",
    grepl(paste(c("C40", "C41"), collapse = "|"),DiagCode) ~ "Bone/Articular cartilage",
    grepl(paste(c("C43", "C44"), collapse = "|"),DiagCode) ~ "Skin",
    grepl(paste(c("C45", "C46", "C47", "C48", "C49"), collapse = "|"),DiagCode) ~ "Soft tissue",
    grepl(paste(c("C50"), collapse = "|"),DiagCode) ~ "Breast",
    grepl(paste(c("C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58"), collapse = "|"),DiagCode) ~ "Female genital organs",
    grepl(paste(c("C60","C61","C62","C63"), collapse = "|"),DiagCode) ~ "Male genital organs",
    grepl(paste(c("C64","C65","C66","C67","C68"), collapse = "|"),DiagCode) ~ "Urinary tract",
    grepl(paste(c("C69","C70","C71","C72"), collapse = "|"),DiagCode) ~ "Central nervous system",
    grepl(paste(c("C73","C74","C75"), collapse = "|"),DiagCode) ~ "Thyroid/Endocrine",
    grepl(paste(c("C78","C79"), collapse = "|"),DiagCode) ~ "Secondary",
    grepl(paste(c("C76","C77","C80"), collapse = "|"),DiagCode) ~ "Ill-defined/Unspecified",
    grepl(paste(c("C81", "C82", "C83", "C84", "C85", "C86", "C87", "C88", "C90", "C91","C92","C93","C94","C95","C96"),collapse = "|"),DiagCode) ~ "Lymphoid/Haematopoietic",
    grepl(paste(c("C97"), collapse = "|"),DiagCode) ~ "Multiple sites"
  )) %>% 
  mutate(cancer_group = case_when(
    grepl(paste(c("C50"), collapse = "|"),DiagCode) ~ "Breast",
    grepl(paste(c("C61"), collapse = "|"),DiagCode) ~ "Prostate",
    grepl(paste(c("C34"), collapse = "|"),DiagCode) ~ "Lung",
    grepl(paste(c("C18", "C19", "C20"), collapse = "|"),DiagCode) ~ "Colorectal",
    grepl(paste(c("C15", "C16","C17","C22","C23","C24","C25"), collapse = "|"),DiagCode) ~ "Upper GI and hepatobiliary",
    grepl(paste(c("C81", "C82", "C83", "C84", "C85", "C86", "C87", "C88", "C90", "C91","C92","C93","C94","C95","C96"),collapse = "|"),DiagCode) ~ "Lymphoid/Haematopoietic",
    grepl(paste(c("C43"), collapse = "|"),DiagCode) ~ "Melanoma",
    grepl(paste(c("C44"), collapse = "|"),DiagCode) ~ "Other skin",
    TRUE~"Others"
  )) 


table(cancer_diagnosis$cancer_group, useNA = "always")


first_cancer_code <- cancer_diagnosis %>% 
  group_by(ClusterID) %>%
  slice_min(EpisodeStart) %>%
  slice_min(DiagNumber) %>%
  slice(1) %>% 
  ungroup() %>%
  rename(FirstCancerDate = EpisodeStart) %>% 
  distinct(ClusterID, FirstCancerDate, DiagCode, cancer_type_sub, cancer_group) 

table(first_cancer_code$cancer_group)


#clean secondary codes
secondary_ID <- first_cancer_code %>% 
  filter(cancer_type_sub == "Secondary") %>% pull(ClusterID)

secondary <- cancer_diagnosis %>% 
  filter(ClusterID %in% secondary_ID) %>% 
  filter(cancer_type_sub != "Secondary") %>% 
  group_by(ClusterID) %>% 
  slice_min(EpisodeStart) %>%
  slice_min(DiagNumber) %>%
  ungroup() %>% 
  select(ClusterID, cancer_type_sub, DiagCode) %>% 
  rename(cancer_type_sub2=cancer_type_sub,
         DiagCode2=DiagCode)


first_cancer_code = first_cancer_code %>% 
  left_join(secondary) %>% 
  mutate(DiagCode3=ifelse(is.na(DiagCode2), DiagCode, DiagCode2)) %>% # if a patient only has "secondary" as first cancer diagnosis then use the secondary diagnosis code
  # as first cancer diagnosis, otherwise, use the one doesn't belong to "secondary" after the "secondary" first cancer diagnosis as first cancer
  select(ClusterID, FirstCancerDate, DiagCode3) %>% 
  rename(DiagCode=DiagCode3) %>% 
  mutate(cancer_group = case_when(
    grepl(paste(c("C50"), collapse = "|"),DiagCode) ~ "Breast",
    grepl(paste(c("C61"), collapse = "|"),DiagCode) ~ "Prostate",
    grepl(paste(c("C34"), collapse = "|"),DiagCode) ~ "Lung",
    grepl(paste(c("C18", "C19", "C20"), collapse = "|"),DiagCode) ~ "Colorectal",
    grepl(paste(c("C15", "C16","C17","C22","C23","C24","C25"), collapse = "|"),DiagCode) ~ "Upper GI and hepatobiliary",
    grepl(paste(c("C81", "C82", "C83", "C84", "C85", "C86", "C87", "C88", "C90", "C91","C92","C93","C94","C95","C96"),collapse = "|"),DiagCode) ~ "Lymphoid/Haematopoietic",
    grepl(paste(c("C43"), collapse = "|"),DiagCode) ~ "Melanoma",
    grepl(paste(c("C44"), collapse = "|"),DiagCode) ~ "Other skin",
    TRUE~"Others"
  )) 

table(first_cancer_code$cancer_group, useNA = "always")

df_episode <- df_episode %>% 
  left_join(first_cancer_code)


df_episode <- df_episode %>% 
  filter(FirstCancerDate>="2000-04-01"&FirstCancerDate<="2025-03-17") %>% 
  mutate(Age_diag = as.numeric(as_date(FirstCancerDate)-as_date(LinkedBirthMonth))/365.24) %>% 
  mutate(Age_adm = as.numeric(as_date(AdmissionDate)-as_date(LinkedBirthMonth))/365.24) %>% 
  filter(Age_diag>=16) # Remove the patients without any cancer

write_csv(df_episode, '/home/anniehu/df_episode.csv')
df= df_episode %>% 
  distinct(ClusterID, LinkedBirthMonth, LinkedSex, ethnicity, FirstCancerDate, DiagCode, cancer_group) 

table(df$cancer_group)
write_csv(df,"df_patient.csv")

### antibiotics
unique(antibiotic$Drug)
unique(micro$DrugName)
df_abx = antibiotic %>% 
  mutate(Drug = str_to_sentence(Drug)) %>% 
  mutate(Drug = str_replace(Drug, " pcc", "")) %>% 
  mutate(Drug = str_replace(Drug, " ophthalmic", "")) %>% 
  mutate(Drug = str_replace(Drug, " \\(name check\\)", "")) %>% 
  mutate(Drug = str_replace(Drug, "/", " + ")) %>% 
  mutate(Drug = str_replace(Drug, "tazobactam .*", "tazobactam")) %>% 
  mutate(Drug = str_replace(Drug, " topical", "")) %>% 
  mutate(Drug = case_when(
    str_detect(Drug, "mbisome") ~ "Ambisome",
    Drug == "Benzylpenicillin sodium" ~ "Benzylpenicillin",
    Drug == "Co-amoxiclav (penicillin base)" ~ "Co-amoxiclav",
    Drug == "Colistin (colistimethate sodium)" ~ "Colistin",
    Drug == "Colistimethate sodium" ~ "Colistin",
    Drug == "Tazocin (piperacillin + tazobactam)" ~ "Piperacillin + tazobactam",
    Drug == "Moxifloxicin" ~ "Moxifloxacin",
    Drug == "Sodium fusidate" ~ "Fusidic acid",
    str_detect(Drug, "Erythromycin") ~ "Erythromycin",
    TRUE ~ Drug
  )) %>% 
  mutate(start = FirstAdministrationDateTime,
         end = LastAdministrationDateTime) %>% 
  filter(Drug != "Dalteparin", Drug != "Amantadine", Drug != "Pabrinex: wernicke's encephalopathy") %>% 
  select(ClusterID, Drug, Route, PrescriptionType, start, end) 
unique(df_abx$Drug)
# remove records without a route or start or end date
df_use = df_abx %>% 
  filter(!is.na(Drug), !is.na(Route), !is.na(start), !is.na(end))
# only keep antibacterials
drug_remove_list = c("Abacavir + lamivudine", "Aciclovir", "Albendazole", "Ambisome", 
                     "Amphotericin b (fungizone)", "Artemether + lumefantrine", "Artesunate", 
                     "Atovaquone", "Bedaquiline", "Caspofungin", "Chloroquine phosphate", 
                     "Chloroquine sulphate", "Cidofovir", "Clofazimine", "Clotrimazole", 
                     "Cycloserine", "Dapsone", "Ethambutol", "Famciclovir", "Fluconazole", 
                     "Flucytosine", "Foscarnet", "Ganciclovir", "Griseofulvin", 
                     "Hydrocortisone + clotrimazole", "Hydrocortisone + gentamicin", 
                     "Isavuconazole", "Isoniazid", "Itraconazole", "Ivermectin", 
                     "Ketoconazole", "Mebendazole", "Mefloquine", "Micafungin", 
                     "Nitazoxanide", "Oseltamivir", "Pentamidine", "Posaconazole", 
                     "Praziquantel", "Primaquine", "Proguanil", 
                     "Proguanil + atovaquone", "Pyrazinamide", 
                     "Pyrimethamine", "Remdesivir", "Rifabutin", 
                     "Rifampicin + isoniazid", "Rifampicin + isoniazid + pyrazinamide", 
                     "Sodium stibogluconate", "Sulfadiazine", "Terbinafine", 
                     "Tinidazole", "Tobramycin", "Valaciclovir", 
                     "Valganciclovir", "Voriconazole", "Zanamivir")
df_use = df_use %>% 
  filter(!Drug %in% drug_remove_list)

#clean up drug names - do this another time
drug_list = df_use %>% count(Drug)


dfc <- df_use %>%
  mutate(
    start = na_if(start, "NULL"),  # Replace "NULL" with actual NA
    start = ymd_hms(start, truncated = 3),  # Parse datetime with fractional seconds
    end = na_if(end, "NULL"),  
    end = ymd_hms(end, truncated = 3)  
  )%>% 
  filter(!is.na(start)&!is.na(end))

dfc <- dfc %>%
  mutate(
    prev_end_24h = as.Date(c(NA, as.Date((end[1:(nrow(dfc)-1)])) ))+ 1
  )

dfc$lag_ClusterID = c(NA, dfc$ClusterID[1:(nrow(dfc)-1)])
dfc$lag_Drug = c(NA, dfc$Drug[1:(nrow(dfc)-1)])
dfc$lag_Route = c(NA, dfc$Route[1:(nrow(dfc)-1)])
dfc$same_gp = dfc$lag_ClusterID == dfc$ClusterID & dfc$lag_Drug == dfc$Drug & dfc$lag_Route == dfc$Route
dfc$same_gp[1] = FALSE
dfc$same_gp_new_episode = ifelse(dfc$same_gp & dfc$start < dfc$prev_end_24h, 0, 1) # If same gp's prescription and start is within 24hrs after the last end date,
# then we don't consider it as new episode.
dfc$episode = cumsum(dfc$same_gp_new_episode)

abx_final = dfc %>% 
  group_by(ClusterID, Drug, Route, episode) %>% 
  summarise(start = min(start), end = max(end)) %>% 
  ungroup()

write_csv(abx_final,"/home/anniehu/AMR Cancer/abx_final.csv")