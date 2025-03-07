import pandas as pd
def msi(msi_file: str) -> str:
    """Parse MSI file and return the MSI score."""
    msi = pd.read_csv(msi_file, sep="\t")
    msi_score = msi["%"].values[0]
    return msi_score

def hrd(hrd_file: str) -> str:
    """Parse HRD file and return the HRD score."""
    hrd = pd.read_csv(hrd_file, sep="\t")
    hrd_score = hrd["HRD_score"].values[0]
    return hrd_score

def tmb(tmb_file: str) -> str:
    """Parse TMB file and return the TMB score."""
    with open(tmb_file, 'r') as file:
        for line in file:
            if line.startswith("TMB:"):
                tmb_score = line.split()[1]
                return tmb_score
    return None
