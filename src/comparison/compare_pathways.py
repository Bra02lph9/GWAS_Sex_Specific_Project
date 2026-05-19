import pandas as pd


female_df = pd.read_csv("../GWAS_Sex_Specific_Project/3_tools_results/gprofiler/female_pathways.tsv")
male_df = pd.read_csv("../GWAS_Sex_Specific_Project/3_tools_results/gprofiler/male_pathways.tsv")

columns_to_keep = [
    "term_name",
    "source",
    "adjusted_p_value",
    "intersection_size"
]

female_df = female_df[columns_to_keep]
male_df = male_df[columns_to_keep]

female_pathways = set(female_df["term_name"])
male_pathways = set(male_df["term_name"])

shared_pathways = female_pathways.intersection(male_pathways)

female_specific = female_pathways - male_pathways
male_specific = male_pathways - female_pathways


shared_df = pd.concat([
    female_df[female_df["term_name"].isin(shared_pathways)],
    male_df[male_df["term_name"].isin(shared_pathways)]
]).drop_duplicates()

female_specific_df = female_df[
    female_df["term_name"].isin(female_specific)
]

male_specific_df = male_df[
    male_df["term_name"].isin(male_specific)
]

shared_df.to_csv("shared_pathways.csv", index=False)

female_specific_df = female_specific_df.sort_values("adjusted_p_value")
female_specific_df.to_csv(
    "female_specific_pathways.csv",
    index=False
)

male_specific_df = male_specific_df.sort_values("adjusted_p_value")
male_specific_df.to_csv(
    "male_specific_pathways.csv",
    index=False
)


print("===== PATHWAY COMPARISON =====")

print(f"Female pathways: {len(female_pathways)}")
print(f"Male pathways: {len(male_pathways)}")

print(f"Shared pathways: {len(shared_pathways)}")
print(f"Female-specific pathways: {len(female_specific)}")
print(f"Male-specific pathways: {len(male_specific)}")

print("\nFiles generated:")
print("- shared_pathways.csv")
print("- female_specific_pathways.csv")
print("- male_specific_pathways.csv")