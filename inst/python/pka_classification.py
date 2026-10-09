# pka_classification.py
#
# Implements the pKa site classification method used to build the
# mp.prediction training dataset:
#   1. For each ionizable site pkasolver reports, classify it as ACID or
#      BASE by comparing formal charge at the reaction-center atom between
#      the protonated and deprotonated structures.
#   2. Identify the specific functional group at that atom via SMARTS
#      substructure matching (checked in the protonated structure for
#      acids, deprotonated structure for bases), using a hand-curated,
#      literature-sourced pattern library (see SOURCES below) -- NOT
#      Dimorphite-DL's bundled table.
#   3. Only sites with an identifiable, approved functional group are
#      kept -- the charge comparison fixes the ACID/BASE label, but does
#      not by itself qualify a site for inclusion.
#   4. Sites are further restricted to pKa in the physiologically relevant
#      range 3-11.
#   5. If multiple valid sites remain: keep the lowest pKa among acidic
#      sites and the highest pKa among basic sites.
#
# SOURCES for SMARTS patterns (see inline comments per pattern):
#   [RDKit]     - rdkit.Chem.FunctionalGroups.BuildFuncGroupHierarchy()
#   [SMARTS-RX] - Journal of Cheminformatics (2025), MolecularAI/smartsrx GitHub
#   [Brenk]     - Brenk et al., ChemMedChem 2008, Supplementary Table S1
#   [Daylight]  - Daylight SMARTS Tutorial (official SMARTS specification)

#

#

import contextlib
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from rdkit import Chem
from pkasolver.query import calculate_microstate_pka_values


def _env_python():
    """Path to the python executable of the CURRENT conda environment."""
    if os.name == "nt":
        candidate = os.path.join(sys.exec_prefix, "python.exe")
    else:
        candidate = os.path.join(sys.exec_prefix, "bin", "python")
    return candidate if os.path.exists(candidate) else sys.executable


@contextlib.contextmanager
def _pkasolver_runtime():
    """Make pkasolver's Dimorphite-DL subprocess call work everywhere.

    pkasolver.query._call_dimorphite_dl() has three portability problems,
    all worked around here WITHOUT modifying pkasolver:

    1. It calls subprocess.run(["python", ...], stderr=STDOUT) without
       setting stdout/stdin. Inside RStudio on Windows there is no console
       handle to inherit, so Windows raises
       "OSError: [WinError 6] The handle is invalid". We default stdin to
       DEVNULL and stdout to PIPE.
    2. It launches the literal command "python", which may be a different
       interpreter (or the Microsoft Store stub) on another machine. We
       substitute this environment's own interpreter.
    3. The helper script writes its results to "test.pkl" in the current
       working directory, which pkasolver then reads back. We run inside a
       private temp directory so nothing is written into the user's
       working directory, and clean up afterwards.
    """
    original_run = subprocess.run
    old_cwd = os.getcwd()
    tmp_dir = tempfile.mkdtemp(prefix="mp_predict_")
    env_python = _env_python()

    def patched_run(cmd, *args, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd and cmd[0] == "python":
            cmd = [env_python] + list(cmd[1:])
        if kwargs.get("stdin") is None:
            kwargs["stdin"] = subprocess.DEVNULL
        if kwargs.get("stdout") is None:
            kwargs["stdout"] = subprocess.PIPE
        return original_run(cmd, *args, **kwargs)

    subprocess.run = patched_run
    os.chdir(tmp_dir)
    try:
        yield
    finally:
        subprocess.run = original_run
        os.chdir(old_cwd)
        shutil.rmtree(tmp_dir, ignore_errors=True)

PKA_RANGE_MIN = 3
PKA_RANGE_MAX = 11

ACIDIC_PATTERNS = {
    "carboxylate": "C(=O)[O;H,-]", 
    "phenol": "[O;H1;$(O-!@c)]",  
    "sulfonamide": "[N;$(NS(=O)(=O))]", 
    "heterocyclic_N_acidic": "[n;D2;H1;$(n1cccc1)]",  
    "hydroxamate": "C(=O)N[OH]",  # [Brenk]
    "phosphate": "[O;$([O;D1]P(=O));!$(OP(=O)(=O))]",  
    "tetrazole": "[nH]1nnnc1",  
    "thiol": "[SH]",  # [Brenk]
    "alcohol": "[O;H1;$(O-!@[C;!$(C=!@[O,N,S])])]",  
    "acidic_amide": "[CX3](=O)[NX3H2,NX3H1]", 
    "acidic_aniline": "[NX3H2,NX3H1][cX3]:[c]",  
    "carbamate": "[N;!$(NC(=O)OC(C)(C)C);!$(NC(=O)OC1c2ccccc2-c3c1cccc3);$(NC(=O)[O;D2])]",  
    "hydrazide": "C(=O)N[NH2]",  
    "imide": "[N;$(N(C(=O)[#6])C(=O)[#6])]",  
    "sulfate": "OS(=O)(=O)[O-]",  
    "carbonic_acid": "[CX3](=[OX1])([OX2H])[OX2H]", 
}

BASIC_PATTERNS = {
    "heterocyclic_N_basic_6ring": "[n;D2;$(n1ccccc1)]", 
    "heterocyclic_N_basic_polyhet": "[n;D2;$(n1aaaaa1);!$(n1ccccc1)]", 
    "heterocyclic_N_basic_5ring_polyhet": "[n;H0;$(n1aaaa1);!$(n1cccc1)]",  
    "aliphatic_amine": "[N;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",  
    "guanidine": "[N;$(NC(~N)=N)]",  
    "amidine": "[N;!R;$(N([#6])=[C;$(C(~N)~N);!$(C(~N)(~N)~N)]),$([N;D1]=[C;$(C(~N)~N);!$(C(~N)(~N)~N)])]",  
    "aniline": "c1cc([NH2])ccc1",  
    
}

_ACIDIC_COMPILED = {name: Chem.MolFromSmarts(smarts) for name, smarts in ACIDIC_PATTERNS.items()}
_BASIC_COMPILED = {name: Chem.MolFromSmarts(smarts) for name, smarts in BASIC_PATTERNS.items()}


def identify_functional_group(mol, atom_idx, compiled_patterns):
    """Return the name of the first functional group pattern that atom_idx
    belongs to, or None if no pattern matches."""
    for fg_name, pattern in compiled_patterns.items():
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            if atom_idx in match:
                return fg_name
    return None


def classify_site(state):
    """Classify a single pkasolver microstate as ACID or BASE, identify its
    functional group, and apply the range/functional-group filter.

    Returns a dict with keys: pka, direction ("ACID"/"BASE"), functional_group
    (str or None), method ("SMARTS" or "charge_fallback"), kept (bool).
    """
    idx = state.reaction_center_idx
    pka = float(state.pka)

    prot_charge = state.protonated_mol.GetAtomWithIdx(idx).GetFormalCharge()
    deprot_charge = state.deprotonated_mol.GetAtomWithIdx(idx).GetFormalCharge()

    if prot_charge > deprot_charge:
        # losing a proton relieves a positive charge -> basic site
        charge_direction = "BASE"
    else:
        # losing a proton creates/reveals a negative charge -> acidic site
        charge_direction = "ACID"

    fg_acid = identify_functional_group(state.protonated_mol, idx, _ACIDIC_COMPILED)
    fg_base = identify_functional_group(state.deprotonated_mol, idx, _BASIC_COMPILED)

    if fg_acid is not None:
        direction, fg, method = "ACID", fg_acid, "SMARTS"
    elif fg_base is not None:
        direction, fg, method = "BASE", fg_base, "SMARTS"
    else:
        # no functional group match at all -- charge fallback fixes the
        # label but the site is NOT kept (no identifiable functional group)
        direction, fg, method = charge_direction, None, "charge_fallback"

    in_range = PKA_RANGE_MIN <= pka <= PKA_RANGE_MAX
    kept = (fg is not None) and in_range

    return {
        "pka": pka,
        "direction": direction,
        "functional_group": fg,
        "method": method,
        "kept": kept,
    }


def extract_relevant_pka(smiles):
    """Run pkasolver on a SMILES string, classify every predicted site, and
    return the final filtered acid/base pKa values plus supporting detail.

    Returns a dict with keys:
      pka_acid, pka_base   - float or None
      fg_acid, fg_base     - str or None (functional group name)
      all_sites            - list of per-site classification dicts (for
                              transparency / debugging), including rejected
                              sites
      error                - None, or a short error code string
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "pka_acid": None, "pka_base": None,
            "fg_acid": None, "fg_base": None,
            "all_sites": [],
            "error": "invalid_smiles",
        }

    states = calculate_microstate_pka_values(mol, only_dimorphite=False)

    sites = [classify_site(s) for s in states]

    acid_candidates = [s for s in sites if s["kept"] and s["direction"] == "ACID"]
    base_candidates = [s for s in sites if s["kept"] and s["direction"] == "BASE"]

    best_acid = min(acid_candidates, key=lambda s: s["pka"]) if acid_candidates else None
    best_base = max(base_candidates, key=lambda s: s["pka"]) if base_candidates else None

    return {
        "pka_acid": best_acid["pka"] if best_acid else None,
        "pka_base": best_base["pka"] if best_base else None,
        "fg_acid": best_acid["functional_group"] if best_acid else None,
        "fg_base": best_base["functional_group"] if best_base else None,
        "all_sites": sites,
        "error": None,
    }
