import pandas as pd
from pathlib import Path

df_af2 = pd.read_csv("../data/af2_metrics.tsv", sep="\t")
df_af3 = pd.read_csv("../data/af3_predictions_5_seeds_sorted.tsv", sep="\t")

df_af2_dockq = pd.read_csv("../data/alphafold2multimer_25_predictions_sr.tsv", sep="\t")
df_af3_dockq = pd.read_csv("../data/alphafold3_25_predictions_sr.tsv", sep="\t")

pdbs = sorted(df_af2["pdb"].unique())
af3_pdbs = (df_af3["pdb"].unique())
print(f"Number of pdbs in AF2: {len(pdbs)}")
print(f"Number of pdbs in AF3: {len(af3_pdbs)}")
# what is not part of what?
af3_not_af2 = []
for pdb in af3_pdbs:
    if pdb not in pdbs:
        af3_not_af2.append(pdb)
print(f" pdbs in AF3 not in AF2: {af3_not_af2}")

pdbs_capri = sorted(df_af2_dockq["pdb"].unique())
if not pdbs == pdbs_capri:
    print("PDBs are not the same")
    for n in range(max(len(pdbs), len(pdbs_capri))):
        if n < len(pdbs):
            if pdbs[n] not in pdbs_capri:
                print(f"{pdbs[n]} not in pdbs_capri")
        if n < len(pdbs_capri):
            if pdbs_capri[n] not in pdbs:
                print(f"{pdbs_capri[n]} not in pdbs")

assert len(pdbs) == 40, f"Expected 40 pdbs, got {len(pdbs)}"

af3_dockqs = []
af3_iptms = []
af3_rscore = []

af2_dockqs = []
af2_iptms = []

low_iptm_af3 = []
low_iptm_af2 = []
af3_failed_pdbs = []
af2_failed_pdbs = []

for pdb in pdbs_capri:
    if pdb == "7nbbD" or pdb == "7te8B":
        continue
    if pdb.startswith("7qne") or pdb.startswith("7qbg"):
        print(f"{pdb} has multi-chain antigen")
        # alphafold failed on both multi-chain antigens
        af2_failed_pdbs.append(pdb)
        continue
    # af2_dockq = df_af2_dockq[df_af2_dockq["pdb"] == pdb]["dockq"].values[0]
    af2_t10 = df_af2_dockq[(df_af2_dockq["pdb"] == pdb) & (df_af2_dockq["rank"] < 11)]
    af2_dockq_t10 = af2_t10["dockq"].max()
    # look at the failed structures
    af2_failed = True
    for i in range(1, 11):
        af2_t10_model = af2_t10[af2_t10["rank"] == i]
        if af2_t10_model["fnat"].values[0] > 0.1 and (af2_t10_model["lrmsd"].values[0] < 10 or af2_t10_model["irmsd"].values[0] < 4):
            af2_failed = False
            break
    if af2_failed:
        af2_failed_pdbs.append(pdb)

    af2_iptm, af2_ptm = df_af2[df_af2["pdb"] == pdb]["iptm"].values[0], df_af2[df_af2["pdb"] == pdb]["ptm"].values[0]
    af2_dockqs.append(af2_dockq_t10)
    af2_iptms.append(af2_iptm)

    af3_failed = True
    af3_dockq_t10 = df_af3_dockq[(df_af3_dockq["pdb"] == pdb) & (df_af3_dockq["rank"] < 11)]["dockq"].max()
    af3_iptm, af3_ptm = df_af3[df_af3["pdb"] == pdb]["iptm"].values[0], df_af3[df_af3["pdb"] == pdb]["ptm"].values[0]
    af3_dockqs.append(af3_dockq_t10)
    # print(f"{pdb} has AF3 iptm {af3_iptm} and af3_dockq_t10 {af3_dockq_t10} while AF2 iptm is {af2_iptm} and af2_dockq_t10 is {af2_dockq_t10}")
    af3_iptms.append(af3_iptm)
    af3_rscore.append(0.8*af3_iptm + 0.2*af3_ptm)
    for i in range(1, 11):
        af3_t10_model = df_af3_dockq[(df_af3_dockq["pdb"] == pdb) & (df_af3_dockq["rank"] == i)]
        if af3_t10_model["fnat"].values[0] > 0.1 and (af3_t10_model["lrmsd"].values[0] < 10 or af3_t10_model["irmsd"].values[0] < 4):
            af3_failed = False
            break
    if af3_failed:
        af3_failed_pdbs.append(pdb)
    # else:
    #     af3_dockq_t10 = "NA"
    #     af3_iptm, af3_ptm = "NA", "NA"

    # print(f"{pdb}\t{af2_dockq_t10}\t{af2_iptm}\t{af3_dockq_t10}\t{af3_iptm}")
    
    if af2_iptm < 0.6:
        #print(f"{pdb} has low AF2 iptm score")
        low_iptm_af2.append(pdb)
    # if af3_iptm is a float
    if af3_iptm != "NA" and af3_iptm < 0.6:
        #print(f"{pdb} has low AF3 iptm score")
        low_iptm_af3.append(pdb)

# calculate correlation between AF3 and iptm
from scipy.stats import spearmanr
corr_iptm = spearmanr(af3_iptms, af3_dockqs)
corr_rscore = spearmanr(af3_rscore, af3_dockqs)
print(f"Correlation between iptm and AF3 dockq: {corr_iptm}")
print(f"Correlation between rscore and AF3 dockq: {corr_rscore}")
# scatter plot for iptm
import matplotlib.pyplot as plt
plt.scatter(af3_iptms, af3_dockqs)
plt.xlabel("IPTM")
plt.title("AF3 IPTM vs DOCKQ-T1")
plt.ylabel("DOCKQ in the top 10 models")
plt.savefig(Path("figures", "alphafold", "af3_iptm_vs_dockq.png"), dpi=300)
plt.close()
# correlations for alphafold multimer
af2_corr_iptm = spearmanr(af2_iptms, af2_dockqs)
print(f"Correlation between iptm and AF2 dockq: {af2_corr_iptm}")
plt.scatter(af2_iptms, af2_dockqs)
plt.xlabel("IPTM")
plt.ylabel("DOCKQ in the top 1 model")
plt.title("AF2 IPTM vs DOCKQ-T1")
plt.savefig(Path("figures", "alphafold", "af2_iptm_vs_dockq.png"), dpi=300)

# failed_af2: should be 73%
assert round(len(af2_failed_pdbs)/len(pdbs_capri),3) == 0.725, f"Expected 72.5% failed structures, got {len(af2_failed_pdbs)/len(pdbs)}"
print(f"Number of failed AF2 structures: {len(af2_failed_pdbs)}")
print(f"Number of low AF2 iptm scores: {len(low_iptm_af2)}")
print(f"Number of low AF3 iptm scores: {len(low_iptm_af3)}")
print(f"Number of failed AF3 structures: {len(af3_failed_pdbs)}")

# now let's iterate over the low_iptm values and check the HADDOCK values
haddock_df_loose_unb = pd.read_csv("../data/capri/loose_capri_predictions_semibound.tsv", sep="\t")
haddock_df_twohit_unb = pd.read_csv("../data/capri/twohit_capri_predictions_semibound.tsv", sep="\t")
haddock_df_real_unb = pd.read_csv("../data/capri/real_capri_predictions_semibound.tsv", sep="\t")

loose_unb_t10 = 0
twohit_unb_t10 = 0
real_unb_t10 = 0
loose_unb_t1 = 0
twohit_unb_t1 = 0
real_unb_t1 = 0

for pdb in low_iptm_af2:
    dock_loose_unb_top10 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] < 11)]
    # print(dock_loose_unb_top10)
    if "acceptable" in dock_loose_unb_top10["capri"].values or "medium" in dock_loose_unb_top10["capri"].values or "high" in dock_loose_unb_top10["capri"].values:
        max_dockq = dock_loose_unb_top10["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_loose_unb_top10}) in the loose scenario")
        loose_unb_t10 += 1
    dock_twohit_unb_top10 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] < 11)]
    if "acceptable" in dock_twohit_unb_top10["capri"].values or "medium" in dock_twohit_unb_top10["capri"].values or "high" in dock_twohit_unb_top10["capri"].values:
        max_dockq = dock_twohit_unb_top10["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_twohit_unb_top10}) in the two-hit scenario")
        twohit_unb_t10 += 1
    dock_real_unb_top10 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] < 11)]
    if "acceptable" in dock_real_unb_top10["capri"].values or "medium" in dock_real_unb_top10["capri"].values or "high" in dock_real_unb_top10["capri"].values:
        max_dockq = dock_real_unb_top10["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_real_unb_top10}) in the real scenario")
        real_unb_t10 += 1
    # now only with the best ranked model
    dock_loose_unb_top1 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] == 1)]
    if "acceptable" in dock_loose_unb_top1["capri"].values or "medium" in dock_loose_unb_top1["capri"].values or "high" in dock_loose_unb_top1["capri"].values:
        max_dockq = dock_loose_unb_top1["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_loose_unb_top1}) in the loose scenario")
        loose_unb_t1 += 1
    dock_twohit_unb_top1 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] == 1)]
    if "acceptable" in dock_twohit_unb_top1["capri"].values or "medium" in dock_twohit_unb_top1["capri"].values or "high" in dock_twohit_unb_top1["capri"].values:
        max_dockq = dock_twohit_unb_top1["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_twohit_unb_top1}) in the two-hit scenario")
        twohit_unb_t1 += 1
    dock_real_unb_top1 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] == 1)]
    if "acceptable" in dock_real_unb_top1["capri"].values or "medium" in dock_real_unb_top1["capri"].values or "high" in dock_real_unb_top1["capri"].values:
        max_dockq = dock_real_unb_top1["dockq"].max()
        # print(f"{pdb} has low AF2 iptm score but high HADDOCK dockq ({dock_real_unb_top1}) in the real scenario")
        real_unb_t1 += 1

print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the loose scenario (T10): {loose_unb_t10/len(low_iptm_af2):.2f}")
print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the two-hit scenario (T10): {twohit_unb_t10/len(low_iptm_af2):.2f}")
print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the real scenario (T10): {real_unb_t10/len(low_iptm_af2):.2f}")
print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the loose scenario with the best model (T10): {loose_unb_t1/len(low_iptm_af2):.2f}")
print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the two-hit scenario with the best model (T10): {twohit_unb_t1/len(low_iptm_af2):.2f}")
print(f"Number of low AF2 iptm scores with high HADDOCK dockq in the real scenario with the best model (T10): {real_unb_t1/len(low_iptm_af2):.2f}")

# same with failed AF2 structures
loose_unb_t10 = 0
twohit_unb_t10 = 0
real_unb_t10 = 0
loose_unb_t1 = 0
twohit_unb_t1 = 0
real_unb_t1 = 0

for pdb in af2_failed_pdbs:
    dock_loose_unb_top10 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] < 11)]
    if "acceptable" in dock_loose_unb_top10["capri"].values or "medium" in dock_loose_unb_top10["capri"].values or "high" in dock_loose_unb_top10["capri"].values:
        max_dockq = dock_loose_unb_top10["dockq"].max()
        loose_unb_t10 += 1
    dock_twohit_unb_top10 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] < 11)]
    if "acceptable" in dock_twohit_unb_top10["capri"].values or "medium" in dock_twohit_unb_top10["capri"].values or "high" in dock_twohit_unb_top10["capri"].values:
        max_dockq = dock_twohit_unb_top10["dockq"].max()
        twohit_unb_t10 += 1
    dock_real_unb_top10 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] < 11)]
    if "acceptable" in dock_real_unb_top10["capri"].values or "medium" in dock_real_unb_top10["capri"].values or "high" in dock_real_unb_top10["capri"].values:
        max_dockq = dock_real_unb_top10["dockq"].max()
        real_unb_t10 += 1
    # now only with the best ranked model
    dock_loose_unb_top1 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] == 1)]
    if "acceptable" in dock_loose_unb_top1["capri"].values or "medium" in dock_loose_unb_top1["capri"].values or "high" in dock_loose_unb_top1["capri"].values:
        max_dockq = dock_loose_unb_top1["dockq"].max()
        loose_unb_t1 += 1
    dock_twohit_unb_top1 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] == 1)]
    if "acceptable" in dock_twohit_unb_top1["capri"].values or "medium" in dock_twohit_unb_top1["capri"].values or "high" in dock_twohit_unb_top1["capri"].values:
        max_dockq = dock_twohit_unb_top1["dockq"].max()
        twohit_unb_t1 += 1
    dock_real_unb_top1 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] == 1)]
    if "acceptable" in dock_real_unb_top1["capri"].values or "medium" in dock_real_unb_top1["capri"].values or "high" in dock_real_unb_top1["capri"].values:
        max_dockq = dock_real_unb_top1["dockq"].max()
        real_unb_t1 += 1

print(f"Number of failed AF2 structures with high HADDOCK dockq in the loose scenario (T10): {loose_unb_t10/len(af2_failed_pdbs):.2f}")
print(f"Number of failed AF2 structures with high HADDOCK dockq in the two-hit scenario (T10): {twohit_unb_t10/len(af2_failed_pdbs):.2f}")
print(f"Number of failed AF2 structures with high HADDOCK dockq in the real scenario (T10): {real_unb_t10/len(af2_failed_pdbs):.2f}")
print(f"Number of failed AF2 structures with high HADDOCK dockq in the loose scenario with the best model (T10): {loose_unb_t1/len(af2_failed_pdbs):.2f}")
print(f"Number of failed AF2 structures with high HADDOCK dockq in the two-hit scenario with the best model (T10): {twohit_unb_t1/len(af2_failed_pdbs):.2f}")
print(f"Number of failed AF2 structures with high HADDOCK dockq in the real scenario with the best model (T10): {real_unb_t1/len(af2_failed_pdbs):.2f}")

# same with failed AF3 structures
loose_unb_t10 = 0
twohit_unb_t10 = 0
real_unb_t10 = 0
loose_unb_t1 = 0
twohit_unb_t1 = 0
real_unb_t1 = 0

for pdb in af3_failed_pdbs:
    dock_loose_unb_top10 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] < 11)]
    if "acceptable" in dock_loose_unb_top10["capri"].values or "medium" in dock_loose_unb_top10["capri"].values or "high" in dock_loose_unb_top10["capri"].values:
        max_dockq = dock_loose_unb_top10["dockq"].max()
        loose_unb_t10 += 1
    dock_twohit_unb_top10 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] < 11)]
    if "acceptable" in dock_twohit_unb_top10["capri"].values or "medium" in dock_twohit_unb_top10["capri"].values or "high" in dock_twohit_unb_top10["capri"].values:
        max_dockq = dock_twohit_unb_top10["dockq"].max()
        twohit_unb_t10 += 1
    dock_real_unb_top10 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] < 11)]
    if "acceptable" in dock_real_unb_top10["capri"].values or "medium" in dock_real_unb_top10["capri"].values or "high" in dock_real_unb_top10["capri"].values:
        max_dockq = dock_real_unb_top10["dockq"].max()
        real_unb_t10 += 1
    # now only with the best ranked model
    dock_loose_unb_top1 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] == 1)]
    if "acceptable" in dock_loose_unb_top1["capri"].values or "medium" in dock_loose_unb_top1["capri"].values or "high" in dock_loose_unb_top1["capri"].values:
        max_dockq = dock_loose_unb_top1
        loose_unb_t1 += 1
    dock_twohit_unb_top1 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] == 1)]
    if "acceptable" in dock_twohit_unb_top1["capri"].values or "medium" in dock_twohit_unb_top1["capri"].values or "high" in dock_twohit_unb_top1["capri"].values:
        max_dockq = dock_twohit_unb_top1
        twohit_unb_t1 += 1
    dock_real_unb_top1 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] == 1)]
    if "acceptable" in dock_real_unb_top1["capri"].values or "medium" in dock_real_unb_top1["capri"].values or "high" in dock_real_unb_top1["capri"].values:
        max_dockq = dock_real_unb_top1
        real_unb_t1 += 1

print(f"Number of failed AF3 structures with high HADDOCK dockq in the loose scenario (T10): {loose_unb_t10/len(af3_failed_pdbs):.2f}")
print(f"Number of failed AF3 structures with high HADDOCK dockq in the two-hit scenario (T10): {twohit_unb_t10/len(af3_failed_pdbs):.2f}")
print(f"Number of failed AF3 structures with high HADDOCK dockq in the real scenario (T10): {real_unb_t10/len(af3_failed_pdbs):.2f}")
print(f"Number of failed AF3 structures with high HADDOCK dockq in the loose scenario with the best model (T10): {loose_unb_t1/len(af3_failed_pdbs):.2f}")
print(f"Number of failed AF3 structures with high HADDOCK dockq in the two-hit scenario with the best model (T10): {twohit_unb_t1/len(af3_failed_pdbs):.2f}")
print(f"Number of failed AF3 structures with high HADDOCK dockq in the real scenario with the best model (T10): {real_unb_t1/len(af3_failed_pdbs):.2f}")

loose_unb_t10 = 0
twohit_unb_t10 = 0
real_unb_t10 = 0
loose_unb_t1 = 0
twohit_unb_t1 = 0
real_unb_t1 = 0

for pdb in low_iptm_af3:
    dock_loose_unb_top10 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] < 11)]
    if "acceptable" in dock_loose_unb_top10["capri"].values or "medium" in dock_loose_unb_top10["capri"].values or "high" in dock_loose_unb_top10["capri"].values:
        max_dockq = dock_loose_unb_top10["dockq"].max()
        loose_unb_t10 += 1
    dock_twohit_unb_top10 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] < 11)]
    if "acceptable" in dock_twohit_unb_top10["capri"].values or "medium" in dock_twohit_unb_top10["capri"].values or "high" in dock_twohit_unb_top10["capri"].values:
        max_dockq = dock_twohit_unb_top10["dockq"].max()
        twohit_unb_t10 += 1
    dock_real_unb_top10 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] < 11)]
    if "acceptable" in dock_real_unb_top10["capri"].values or "medium" in dock_real_unb_top10["capri"].values or "high" in dock_real_unb_top10["capri"].values:
        max_dockq = dock_real_unb_top10["dockq"].max()
        real_unb_t10 += 1
    # now only with the best ranked model
    dock_loose_unb_top1 = haddock_df_loose_unb[(haddock_df_loose_unb["pdb"] == pdb) & (haddock_df_loose_unb["rank"] == 1)]
    if "acceptable" in dock_loose_unb_top1["capri"].values or "medium" in dock_loose_unb_top1["capri"].values or "high" in dock_loose_unb_top1["capri"].values:
        max_dockq = dock_loose_unb_top1["dockq"].max()
        loose_unb_t1 += 1
    dock_twohit_unb_top1 = haddock_df_twohit_unb[(haddock_df_twohit_unb["pdb"] == pdb) & (haddock_df_twohit_unb["rank"] == 1)]
    if "acceptable" in dock_twohit_unb_top1["capri"].values or "medium" in dock_twohit_unb_top1["capri"].values or "high" in dock_twohit_unb_top1["capri"].values:
        max_dockq = dock_twohit_unb_top1["dockq"].max()
        twohit_unb_t1 += 1
    dock_real_unb_top1 = haddock_df_real_unb[(haddock_df_real_unb["pdb"] == pdb) & (haddock_df_real_unb["rank"] == 1)]
    if "acceptable" in dock_real_unb_top1["capri"].values or "medium" in dock_real_unb_top1["capri"].values or "high" in dock_real_unb_top1["capri"].values:
        max_dockq = dock_real_unb_top1["dockq"].max()
        real_unb_t1 += 1

print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the loose scenario (T10): {loose_unb_t10/len(low_iptm_af3):.2f}")
print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the two-hit scenario (T10): {twohit_unb_t10/len(low_iptm_af3):.2f}")
print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the real scenario (T10): {real_unb_t10/len(low_iptm_af3):.2f}")
print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the loose scenario with the best model: {loose_unb_t1/len(low_iptm_af3):.2f}")
print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the two-hit scenario with the best model: {twohit_unb_t1/len(low_iptm_af3):.2f}")
print(f"Number of low AF3 iptm scores with high HADDOCK dockq in the real scenario with the best model: {real_unb_t1/len(low_iptm_af3):.2f}")