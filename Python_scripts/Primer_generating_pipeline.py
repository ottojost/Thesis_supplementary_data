#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script for RPA primer pair generation suited for the DoE screening.

Key operations:
- Applies a penalty scoring system to sample a diverse set of unique primers from modified PrimedRPA output
- Separates forward and reverse primers of each pair
- Create new primer pairs if they are similar in characteristics and form a valid amplicon (100-250bp)
- Integrates PrimedRPA logic to calculate the 'Max Dimerisation Percentage' and calculates other necessary parameters
- Returns Excel file containing primer pairs with similar forward and reverse primers
"""

from pathlib import Path
import pandas as pd
import numpy as np
import sys

# --- Configuration & Paths ---

# Path to folder in which separate folders containing PrimedRPA outputs for each primer length are situated
Input_folder = Path("")

# Create a separate output directory
Output_folder = Input_folder / "Primer_generating_pipeline"
Output_folder.mkdir(exist_ok=True)

# Template on which primers should be generated (required for amplicon features)
Template = 'GCGAAGAACTAGTTAACCGCGTTCTGTGGCTGGTTCTCTGCTCGCCTTCCTTGGCTTCCTGCCCTCTCCACCCTTCATTGTGTTGTGGCATGCTCCGTTACGTAACGTTATTTCCTTCTCCAGCTTCTTGGGCGTTTGATGATTGTTGGCGAGAAGTGTGCTGCTAGCCTGGGCTTGACCGATGGATTCCGGATGGCTGTGAGATACCCACCCTCAGTCCCTTCAGACTACCGCGCGCGGCTCTGTATTCTGGGTGGCCGTCAGTTGGGCCAGCCTCCTGGCTAAGATGTTTGCACCGCCGGTGTTGCTGCACGCGTACGGATCGCCACCGAATGGGTTTCACGTGTTGCCCGTCAGCCTAGCCACCGGTGACATGTAATTGTTTTTGGTGGGTGCCTATGGAGGGTAATGAAAAGCTTTGAGCAGCATTTGCAGAATAAAGATGGAGCATGGGGATATCATCTCTGTCTCGCGTGTCTTACTGATGGGAATGGTTCTGCAAGTTAAAACGGTGGGTCTCCGTCTCCCAAAGATAATATAGGATGTGCTGTGTGTTTGAGGCCTTGCTGTCTTCCCCCAGCTTCCAAATAAAAGGGCGTGGGTTTTATCGCGTACTAGCAGCATCGACAAGCAGCAGGCGTAGAATTAGCGAAGAATAAGGTAGGAAGAAAGCTTCCAGCTGTCGTCTGTTGTATGGCCCCGTCGCGAGGCCAAATAAACTCTACCTGCATGATTTCTGACCGGTGCGCGTCAGCCTGTTCTTAACCAAACTCGGAGTTTTGCCACTGCTTTCTTTTCAGTCGCGGATTTTGCGACGTAGAGAATTAATTTCTCCCGCGGAAACTTTCCCAAATGGAAAAGACTTTGTGATACTTGAGTCGGCGATTGTTGCGCTTCTGCTGTCTTTGAGACTCGGCTGCCCTTGTCCATTTCTCTTTTTTTTACATATGAACAAAACATAAAAAAACAATCAATCAAACAAACAAACAAACAATAAAGCAAAAAAAAACCAAAAAAACAAAAGCCGTAGTATTGTTCCATAACCAGTAGGAGGGGTAGCTACAGGCACCGTAAGGCTTTTGCATGTAATCACCTCTCGTGCTTTGTTAGGCCTCTGGCTTTACTTGAGTCTTTGTTGAATTGCGATTGCATTCAAGGCTTACTGGAGTGGGGGAGTTTCAAAATGCTAATAAATATTTGCATTAAAATCTAAATAGTTGAAACTTGCTAAATGAAATTACTGTGTGGCCGTTCCTCAAGGAATCTCGATTTGTCAGCAGTAGATCCGCTTGCGTGAGTGACGATCGCAGCTGTGAATGACTCTGAATTAAACCGTTTGAAGAAGCGGAATAAGACAACTTGTTACGTTTTTCCTCATATCCAACCTGAAGCTCCCCCCGGCGCAACTTAAGGCCGTTACCTCTTGTCCCGTTGCTAGTTACCTGGGAGAAGCGGCCAACCCTCACCTCGCTACAACCACCTTTCGGGTAGTTGTAAAGAGGCGATAAGGTCTCTCTCCCCTGAGCCTCCCCTCCTCCAGACGAAACAATCCCAGTTCCCTCAACTGCTACCCATAGGTCTTGTGCTCCAGTCCCTTCACCGGCTTCGTTGCCCTCCTCTGGCCACGCTCCAGAGCCTCAACGTCTTTCTTGTAGTGAGGGGCCCAAAGCTGCACGCAGTACTCGAGGTGCGGCCTCACCAGCGCCGAGTACAGAGGGATGATCACCTCCCTCGTCCTGCTGACTCCGTGATTTCTGATCCAAGCCAGGATGCCGTTGGCCTTCTTGGCCCCCTGGGCACGCCGCTGGCTCATGTTCAGGCGGCTGTTGACTAATGCCCCCAGATCCGTTTCTTCCGCGCAGCTTTCCAGCCACTCTGCCCCGAGCCTGTTGCGTTGCATGGGGTTGTTGTGACCAAAGTGCAGGACCCGGCGCTTGGTCTCGTAGATGCTCCTACAATTGGCCTCGGCCCAACGATCCAGCCCATCCGGATCCCTGTGTAGGGCCTTCCTACCCTCAGGCAGATCAAGCCTTCCTCCCAACATGGTGTGGCCTGCAAACTTATTGAGGGTGCGCTCAATCCCCTCGTCCAGATCATCATTCGAGATATTAAACAGGACTGGCCCAAACGCTGACCCCCTGGGGAACACCGCTAGCGACCACCGCCAGCTGGATTGAACTCCATTCAGCACCACTCTCTGGGCCCGGCCGTCCAGCCAGGTTTTGACGCAGCAGAGTACACCTGCCCAAGCCATGGGCTGCCACCTTCTCTGGGAGAATACCGTGGGAGACAGAGTCAAAGGCTTCACTGAAGTCTAGGTAGACTACATCCACAGCCTTTGCCTCATCCACTAGTCGGGTCACCTGGTCATAGAAGGAGATCGGGTCGGTCAAGCAGGACCTGTCTTTAAGGAACCCGCGCTGGCTAGGCCCGGTCCCCTGGTTGTCCCGCACACGCCCTGTGATCTCACTCCGGGTGATCTGCTCCATGACCTTCCCCGATACCGAGGTCAGGCTGACAGGCCTGTAGTTCCCCGGATCCTCCATCCAACCCTTCTTGGAGGTGGGCGTCTGTCTGGCTCCAAGTATCGGTAGCAAGACGGTTCGTTATATACGTATAGCCAACGTGAGATTAAAACTCAGAAAGGAGACAAATATTTGCCCAACCGGTTTATTGTAAAACGTTTTGGTTAGCTAGGCCTGGGGATGGGAGGTGGGATGACGAAAGCAATAGTGAAAAGACAGAAAAGAAAAGCAAGGATTTGTTACAGGAATGCAGAGATGGTTACCACCGCAGTGTGACAATGTTTCATGATTCCGGTTGGAAGC'

# Parameters & Constraints
MIN_AMPLICON_LEN, MAX_AMPLICON_LEN = 100, 250
MIN_AMPLICON_GC, MAX_AMPLICON_GC = 30, 70
MIN_AMPLICON_LC, MAX_AMPLICON_LC = 70, 100

MIN_PRIMER_LEN, MAX_PRIMER_LEN = 27, 36
MIN_PRIMER_GC, MAX_PRIMER_GC = 30, 70
MIN_PRIMER_LC, MAX_PRIMER_LC = 70, 100
MIN_PRIMER_HP_Tm, MAX_PRIMER_HP_Tm = 0, 80

# Sampling Constraints
SAMPLE_SIZE_PER_FILE = 1000 
MIN_DISTANCE_THRESHOLD = 5

# Logging configuration
log_messages = []
def log(message):
    print(message)
    log_messages.append(str(message))

# --- Helper Functions ---

def SSIdentification(SeqOne, SeqTwo, ReverseOrientation, FixedBack=False, threshold=4):
    """
    Calculates the secondary structure and dimerisation potential between two sequences.
    Ported from PrimedRPA logic (Higgins M et al. Submitted. 2018).
    """
    NucDict = {'A':'T',
               'T':'A',
               'C':'G',
               'G':'C',
                }

    MaxBindingSites = 0
    MaxBindingString = ''
    MaxBindingPercentage = 0

    MaxHardFailScore = 0
    HardFail = False
    HardFailString = ''

    MaxLength = len(SeqOne) + len(SeqTwo) - 1

    SeqOneNorm = SeqOne + ' '*(MaxLength-len(SeqOne))
    SeqTwoNorm = SeqTwo + ' '*(MaxLength-len(SeqTwo))
    SeqTwoReversed = SeqTwo[::-1] + ' '*(MaxLength-len(SeqTwo))

    if ReverseOrientation == True:
        CombosToCheck = [(SeqTwoNorm,SeqOneNorm),
                        (SeqOneNorm,SeqTwoNorm),
                        (SeqOneNorm,SeqTwoReversed),
                        (SeqTwoReversed,SeqOneNorm)]
    else:
        CombosToCheck = [(SeqTwoNorm,SeqOneNorm),
                        (SeqOneNorm,SeqTwoNorm)]

    for VarPair in CombosToCheck:
        for i in list(range(VarPair[0].count(' ')+1)):

            if i == 0:
                DynamicSeq = VarPair[0]
            else:
                DynamicSeq = ' '*i + VarPair[0][:-i]

            DyanimicSeqList = list(DynamicSeq)
            FixedSeqList = list(VarPair[1])

            PossibleBindingString = ''

            for StringPos in list(range(len(DyanimicSeqList))):
                if DyanimicSeqList[StringPos] in list(NucDict.keys()):
                    if NucDict[DyanimicSeqList[StringPos]] == FixedSeqList[StringPos]:
                        PossibleBindingString += '|'
                    else:
                        PossibleBindingString += '-'
                else:
                    PossibleBindingString += '-'

            NumberOfBindingMatches = PossibleBindingString.count('|')
            CompleteBindingString = DynamicSeq + '\n' + PossibleBindingString + '\n' + VarPair[1]

            if VarPair[1] == SeqOneNorm:
                FiveIndex = 0
                ThreeIndex = len(SeqOne)
            else:
                FiveIndex = i
                ThreeIndex = i + len(SeqOne)

            UpperLimit = FiveIndex+22
            LowerLimit = ThreeIndex-22
            if LowerLimit < 0:
                LowerLimit = 0

            FivePrimeCounts = PossibleBindingString[FiveIndex:UpperLimit]
            ThreePrimeCounts = PossibleBindingString[LowerLimit:ThreeIndex]

            weightings = [3,2,1.5]
            ThreePrimeLoc = [-1,-2,-3]
            FivePrimeLoc = [0,1,2]

            OriginalScore = [FivePrimeCounts,ThreePrimeCounts]
            IndexLocations = [FivePrimeLoc,ThreePrimeLoc]

            AdjustedWeights = []
            for ib in [0,1]:
                TempScore = OriginalScore[ib].count('|') - OriginalScore[ib].count('-')
                Indexes = IndexLocations[ib]
                for iz in list(zip(Indexes,weightings)):
                    if iz[0] < len(OriginalScore[ib]):
                         if OriginalScore[ib][iz[0]] == '|':
                             TempScore += iz[1] - 1
                AdjustedWeights.append(TempScore)

            if max(AdjustedWeights) >= 21.5:
                HardFail = True
                HardFailString = CompleteBindingString

            if MaxBindingSites < NumberOfBindingMatches:
                MaxBindingSites = NumberOfBindingMatches
                MaxBindingString = CompleteBindingString

                if FixedBack == False:
                    Denom = min([len(SeqOne),len(SeqTwo)])
                else:
                    Denom = len(SeqOne)

                MaxBindingPercentage = (MaxBindingSites/Denom)*100

    return (MaxBindingPercentage, MaxBindingString, HardFail, HardFailString)

def calculate_gc_vectorized(series):
    s = series.astype(str).str.upper()
    lengths = s.str.len()
    gc_count = s.str.count('G') + s.str.count('C')
    return np.where(lengths > 0, (gc_count / lengths) * 100, 0)

def three_prime_vectorized(series, window=5):
    s = series.astype(str).str.upper()
    tails = s.str[-window:]
    lengths = tails.str.len()
    gc_counts = tails.str.count('G') + tails.str.count('C')
    gc_perc = np.where(lengths > 0, (gc_counts / lengths) * 100, 0)
    return np.where(gc_perc >= 80, "Yes", "No")

def linguistic_complexity_simple(seq, max_k=6):
    if not isinstance(seq, str): return 0
    seq = seq.upper()
    n = len(seq)
    if n == 0: return 0
    total_obs, total_poss = 0, 0
    for k in range(1, min(max_k, n) + 1):
        observed = len(set(seq[i:i+k] for i in range(n - k + 1)))
        possible = min(4**k, n - k + 1)
        total_obs += observed
        total_poss += possible
    return (total_obs / total_poss)*100 if total_poss > 0 else 0

def select_diverse_primers_scored(df, n_samples, distance_threshold=5):
    """
    Selects primers based on a penalty scoring system to ensure spatial diversity.
    Prioritizes maximum distance between selected binding sites.
    """
    if len(df) <= n_samples:
        log(f"   -> Dataset smaller than requested samples ({len(df)} vs {n_samples}). Retaining all.")
        return df

    df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    fp_starts = df_shuffled["FP Binding Start Site"].values
    rp_starts = df_shuffled["RP Binding Start Site"].values
    
    N = len(df_shuffled)
    
    penalties = np.random.rand(N) * 0.0001 
    selected_indices = []
    
    log(f"   -> Initiating penalty sampling (Target: {n_samples}, Pool: {N})...")
    
    for i in range(n_samples):
        best_idx = np.argmin(penalties)
        selected_indices.append(best_idx)
        
        chosen_fp = fp_starts[best_idx]
        chosen_rp = rp_starts[best_idx]
        
        penalties[best_idx] = np.inf
        
        fp_dist = np.abs(fp_starts - chosen_fp)
        rp_dist = np.abs(rp_starts - chosen_rp)
        
        is_neighbor = (fp_dist < distance_threshold) & (rp_dist < distance_threshold)
        penalties[is_neighbor] += 1.0
        
        if (i + 1) % 1000 == 0:
            log(f"      Selected {i + 1}/{n_samples} primers...")

    log(f"   -> Completed. {len(selected_indices)} primers selected with maximum diversity.")
    return df_shuffled.iloc[selected_indices].copy()

# --- Step 1: Load & Process Data ---

FP_columns = [
    "FP Binding Start Site", "FP GC%", "Forward Primer (FP)", "Forward Primer Length",
    "FP Hairpin", "FP Hairpin Tm", "FP Hairpin Structure", "FP Hairpin MFE", "Set_ID"
]
RP_columns = [
    "RP Binding Start Site", "RP GC%", "Reverse Primer (RP)", "Reverse Primer Length",
    "RP Hairpin", "RP Hairpin Tm", "RP Hairpin Structure", "RP Hairpin MFE", "Set_ID" 
]

xlsx_files = []
log(f"[Step 1] Searching for files in: {Input_folder}")

if Input_folder.exists():
    count_scanned = 0
    for path in Input_folder.rglob("*.xlsx"):
        count_scanned += 1
        name = path.stem.lower()
        is_target = "3_adding_hairpins_and_tm" in name
        in_output_folder = "separated" in str(path).lower()
        if is_target and not in_output_folder:
            xlsx_files.append(path)
else:
    log(f"CRITICAL ERROR: Input directory {Input_folder} does not exist.")
    sys.exit("Script terminated.")

log(f"[Step 1] Found {len(xlsx_files)} valid input files.")
if not xlsx_files:
    sys.exit("Script terminated: No files found.")

all_FP = []
all_RP = []
global_id_counter = 1

for xlsx_path in xlsx_files:
    try:
        df = pd.read_excel(xlsx_path)
    except Exception as e:
        log(f"Error reading {xlsx_path.name}: {e}")
        continue

    if "Set_ID" not in df.columns:
        df.insert(0, "Set_ID", range(1, len(df) + 1))

    missing = [col for col in FP_columns if col not in df.columns]
    if missing:
        continue
    
    df_sample = select_diverse_primers_scored(df, n_samples=SAMPLE_SIZE_PER_FILE, distance_threshold=MIN_DISTANCE_THRESHOLD)
    
    df_FP = df_sample[FP_columns].copy()
    df_RP = df_sample[RP_columns].copy()

    ids = np.arange(global_id_counter, global_id_counter + len(df_sample))
    df_FP.insert(0, "Primer_ID", ["FP_" + str(i) for i in ids])
    df_RP.insert(0, "Primer_ID", ["RP_" + str(i) for i in ids])
    
    global_id_counter += len(df_sample)
    
    df_FP["Primer_LC%"] = df_FP["Forward Primer (FP)"].apply(linguistic_complexity_simple)
    df_FP["GC_enrichment_3p"] = three_prime_vectorized(df_FP["Forward Primer (FP)"])
    
    df_RP["Primer_LC%"] = df_RP["Reverse Primer (RP)"].apply(linguistic_complexity_simple)
    df_RP["GC_enrichment_3p"] = three_prime_vectorized(df_RP["Reverse Primer (RP)"])
    
    all_FP.append(df_FP)
    all_RP.append(df_RP)

if not all_FP:
    sys.exit("No data found.")

final_FP = pd.concat(all_FP, ignore_index=True)
final_RP = pd.concat(all_RP, ignore_index=True)

# --- Step 2: Grouping ---

def Grouping_primers(df, primer_type):
    rename_map = {}
    target = "Forward Primer" if primer_type == "FP" else "Reverse Primer"
    short = "(FP)" if primer_type == "FP" else "(RP)"
    prefix = "FP" if primer_type == "FP" else "RP"
    
    for col in df.columns:
        new = col.replace(target, "Primer").replace(short, "").replace(prefix, "Primer")
        new = new.strip().replace("  ", " ")
        rename_map[col] = new
    
    df = df.rename(columns=rename_map)
    
    if "Primer GC%" in df.columns:
        df["GC_group"] = pd.cut(df["Primer GC%"], bins=range(30, 71, 3), right=False)
    if "Primer Length" in df.columns:
        df["Primer_length_group"] = pd.cut(df["Primer Length"], bins=range(27, 37, 2), right=False)
    
    if "Primer Hairpin Tm" in df.columns:
        df["Primer Hairpin Tm"] = df["Primer Hairpin Tm"].fillna(25)
        df["Hairpin_Tm_group"] = pd.cut(df["Primer Hairpin Tm"], bins=range(25, 101, 5), right=False)
    if "Primer_LC%" in df.columns:
        df["LC_group"] = pd.cut(df['Primer_LC%'], bins=range(70,101,5), right=False)
        
    return df

FP_with_groups = Grouping_primers(final_FP, "FP")
RP_with_groups = Grouping_primers(final_RP, "RP")

# --- Step 3 & 4: Combine & Calculate Dimerisation ---

def Combine_primers(FP_df, RP_df):
    log("[Step 4] Combining primers...")
    
    base_groups = ["GC_group", "Hairpin_Tm_group", "Primer_length_group", "LC_group"]
    extra_keys = ["GC_enrichment_3p"]
    merge_keys = [c for c in base_groups + extra_keys if c in FP_df.columns and c in RP_df.columns]
    
    log(f"Merging on keys: {merge_keys}")
    
    merged = FP_df.merge(RP_df, on=merge_keys, suffixes=(" FP", " RP"))
    
    if merged.empty:
        log("No combinations found.")
        return merged

    for col in extra_keys:
        if col in merged.columns:
            merged[f"{col} FP"] = merged[col]
            merged[f"{col} RP"] = merged[col]

    starts = merged["Primer Binding Start Site FP"].values.astype(int)
    rp_starts = merged["Primer Binding Start Site RP"].values.astype(int)
    rp_lens = merged["Primer Length RP"].values.astype(int)
    ends = rp_starts + rp_lens
    
    merged["Amplicon_Seq"] = [Template[s:e] for s, e in zip(starts, ends)]
    merged["Amplicon_length"] = merged["Amplicon_Seq"].str.len()
    
    merged = merged[merged["Amplicon_length"].between(MIN_AMPLICON_LEN, MAX_AMPLICON_LEN)]
    if merged.empty: return merged
    
    merged["Amplicon_GC"] = calculate_gc_vectorized(merged["Amplicon_Seq"])
    merged = merged[merged["Amplicon_GC"].between(MIN_AMPLICON_GC, MAX_AMPLICON_GC)]
    
    merged["Amplicon_LC"] = merged["Amplicon_Seq"].apply(linguistic_complexity_simple)
    merged = merged[merged["Amplicon_LC"].between(MIN_AMPLICON_LC, MAX_AMPLICON_LC)]
    
    merged = merged[
        (merged["Primer GC% FP"].between(MIN_PRIMER_GC, MAX_PRIMER_GC)) &
        (merged["Primer Length FP"].between(MIN_PRIMER_LEN, MAX_PRIMER_LEN)) &
        (merged["Primer_LC% FP"].between(MIN_PRIMER_LC, MAX_PRIMER_LC)) &
        (merged["Primer Hairpin Tm FP"].between(MIN_PRIMER_HP_Tm, MAX_PRIMER_HP_Tm))
    ].copy()
    
    log(f"Calculating PrimedRPA Dimerisation Scores for {len(merged)} pairs...")
    dimer_scores = []
    
    for idx, row in merged.iterrows():
        fp_seq = row["Primer FP"]
        rp_seq = row["Primer RP"]
        
        perc, string, fail, fail_string = SSIdentification(fp_seq, rp_seq, ReverseOrientation=True)
        dimer_scores.append(perc)
        
    merged["Max Dimerisation Percentage"] = dimer_scores

    desired_order = [
        "Primer_ID FP", "Primer_ID RP", "Set_ID FP",
        "Amplicon_length", "Amplicon_GC", "Amplicon_LC", "Amplicon_Seq",
        "Max Dimerisation Percentage",
        "Primer FP", "Primer GC% FP", "Primer Length FP", "Primer_LC% FP", "Primer Hairpin Tm FP",
        "GC_enrichment_3p FP", "Primer Binding Start Site FP",
        "Primer RP", "Primer GC% RP", "Primer Length RP", "Primer_LC% RP", "Primer Hairpin Tm RP",
        "GC_enrichment_3p RP", "Primer Binding Start Site RP"
    ]
    
    final_cols = [c for c in desired_order if c in merged.columns]
    merged_clean = merged[final_cols].copy()
    
    log(f"[Step 4] {len(merged_clean)} valid pairs extracted.")
    
    output_path = Output_folder / "Generated_primer_pairs.xlsx"
    if not merged_clean.empty:
        merged_clean.to_excel(output_path, sheet_name="All_combinations", index=False)
        
    return merged_clean

combined = Combine_primers(FP_with_groups, RP_with_groups)


