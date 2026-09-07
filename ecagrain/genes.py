"""Genes blocked from driving metacell grouping (never from the aggregated output matrix):
lateral programmes (cell cycle, dissociation stress, heat shock) and technical families (mitochondrial, ribosomal,
Malat1; haemoglobin is kept because erythrocyte datasets exist). Matched case-insensitively so mouse and human symbols both work; the OSP var flags are reused when present."""

from __future__ import annotations

import numpy as np
import pandas as pd
from osp import (
    DISSOCIATION_GENES_HS,
)  # van den Brink et al. 2017 dissociation-stress panel, kept by osp

# Tirosh et al. 2016 cell-cycle markers (the Seurat cc.genes lists).
S_GENES = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
]
G2M_GENES = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]
# Technical families blocked from grouping (case-insensitive). Haemoglobin is deliberately NOT blocked: the atlas
# has erythrocyte datasets and Hba/Hbb are their identity. Heat-shock genes (HSPA/HSPB/HSPD/HSPE/HSPH/HSP90, DNAJA/B)
# are blocked as a stress programme on top of the dissociation panel.
FAMILY_REGEX = {
    "mt": r"^MT-",
    "ribo": r"^(RP[SL]\d+[A-Z]?\d*(L\d*)?|RPLP\d|RPSA)$",  # RPS/RPL proteins and paralogs; not RPS6K kinases or pseudogenes
    "hsp": r"^(HSP[ABDEH]\d|HSP90|DNAJ[AB]\d)",
}


def lateral_mask(var_names) -> np.ndarray:
    """True for cell-cycle and dissociation-stress genes."""
    lateral = {g.upper() for g in (*S_GENES, *G2M_GENES, *DISSOCIATION_GENES_HS)}
    return np.array([str(g).upper() in lateral for g in var_names])


def blocked_mask(var: pd.DataFrame) -> np.ndarray:
    """Lateral genes + mt / ribo / heat-shock / Malat1 (not haemoglobin). Uses the boolean var columns written by OSP when present,
    otherwise the same regexes."""
    names = var.index.astype(str)
    mask = lateral_mask(names)
    for fam, rx in FAMILY_REGEX.items():
        mask |= var[fam].to_numpy(bool) if fam in var else np.asarray(names.str.match(rx, case=False), bool)
    mask |= var["malat1"].to_numpy(bool) if "malat1" in var else np.asarray(names.str.upper() == "MALAT1", bool)
    return mask
