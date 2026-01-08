# src/config.py

import os

BASE_DIR = '~/projects/gpasd'
RAW_DATA_PATH = os.path.join(BASE_DIR, 'data/plp_3726.xlsx')
PHENOTYPE_PATH = os.path.join(BASE_DIR, 'data/spark_clinical_infor_44962.phenotype.csv')

# WGS data
WGS_PHENO_PATH = os.path.join(BASE_DIR, 'data/SPARK_WGS_sample_meta_2237.csv')
WGS_PLP_PATH = os.path.join(BASE_DIR, 'data/wgs.83genes.sample_2237.44plp.tsv') 

OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
PROCESSED_DATA_DIR = os.path.join(OUTPUT_DIR, 'processed_data')

os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

# Traits
# Scales (continuous)
SCQ = ['SCQ_total', 'SCQ_communication', 'SCQ_interaction', 'SCQ_stereotype']
RBSR = ['RBSR_total', 'RBSR_stereotyped', 'RBSR_injurious', 'RBSR_compulsive', 'RBSR_ritualistic', 'RBSR_sameness', 'RBSR_restricted']
DCDQ = ['DCDQ_total', 'DCDQ_movement', 'DCDQ_handwriting', 'DCDQ_coordination']
ABC = ['ABC_total', 'ABC_communication', 'ABC_daily_living_skills', 'ABC_socialization', 'ABC_motor_skills']

# Developmental milestones (continuous)
DEV_MILESTONES = [
    'Smiling', 'Sitting_upright', 'Crawling', 'Walking', 'Spoon_feeding_self',
    'Speaking_first_word', 'Speaking_with_combined_word', 'Speaking_first_phrase',
    'Attaining_bladder_control', 'Attaining_bowel_control'
]

# Mental health (binary)
MH_PC = [
    'ADHD', 'Conduct_disorder', 'Intermittent_explosive_disorder', 'Oppositional_defiant_disorder',
    'Anxiety_disorder', 'Bipolar_disorder', 'Major_depressive_disorder', 'DMDD', 'Hoarding_disorder',
    'Obsessive_compulsive_disorder', 'Separation_anxiety_disorder', 'Social_anxiety_disorder',
    'Schizophrenia', 'Sleep_disorder', 'Mutism', 'Social_communication_disorder', 'Eating_disorder', 'Encopresis',
    'Enuresis'
]

# Neurodevelopmental (binary)
NDD_RELATED = [
    'Tourette_Syndrome', 'Epilepsy', 'Intellectual_disability', 'Language_disorder',
    'Learning_disability', 'DCD', 'Speech_sound_disorder',
    'Hearing_loss.deafness', 'Strabismus', 'Feeding_problems'
]

# Perinatal (binary)
PERINATAL = [
    'Perinatal_complications', 'Fetal_alcohol_syndrome', 'Intracranial_hemorrhage',
    'Birth_hypoxia', 'Prenatal_infection', 'Premature_birth'
]

# Birth defects (binary)
BIRTH_DEFECTS = [
    'Birth_defects', 'Bone.limb_defect', 'CNS_defect',
    'Craniofacial_defect', 'Gastrointestinal_defect', 'Heart.lung_defect', 'Genitourinary_defect'
]

# Growth (binary)
GROWTH = [
    'Growth_abnormalities', 'Macrocephaly', 'Microcephaly',
    'Underweight', 'Obesity', 'Short_stature'
]

CONTINUOUS_TRAITS = SCQ + RBSR + DCDQ + ABC + DEV_MILESTONES
BINARY_TRAITS = MH_PC + NDD_RELATED + PERINATAL + BIRTH_DEFECTS + GROWTH
ALL_TRAITS = CONTINUOUS_TRAITS + BINARY_TRAITS


if __name__ == "__main__":
    print(f"Continuous Traits: {len(CONTINUOUS_TRAITS)}")
    print(f"Binary Traits: {len(BINARY_TRAITS)}")
    print(f"Total Traits: {len(ALL_TRAITS)}")
