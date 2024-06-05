import string
from pathlib import Path

import numpy as np
import pandas as pd


def gaussian_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def create_input_data(
    dir_path, seqs=["0001", "0002"], tps=[0, 3, 6], n_repl=2, half_lives=[3, 1]
):
    n_seqs = len(seqs)
    n_tps = len(tps)

    # create EPGs
    x = np.linspace(0, 3000, 3001)
    epgs = np.zeros((n_repl, n_seqs, n_tps, len(x)))
    for repl in range(n_repl):
        for seq_idx in range(n_seqs):
            for tp_idx, tp in enumerate(tps):
                epgs[repl, seq_idx, tp_idx] = gaussian_pdf(x, 1500, 100) * 0.5 ** (
                    tp / half_lives[seq_idx]
                )
                epgs[repl, seq_idx, tp_idx] += gaussian_pdf(x, 150, 10)  # control peak

    for repl in range(n_repl):
        plate_map = np.empty((n_seqs, n_tps), dtype=object)
        epg_file_content = {"Size (nt)": x}
        for seq_idx, (row, seq_id) in enumerate(
            zip(string.ascii_uppercase[:n_seqs], seqs)
        ):
            for col, tp in enumerate(tps):
                epg_file_content[f"{row}{col + 1}"] = epgs[repl, seq_idx, col].round(4)
                plate_map[seq_idx, col] = f"{seq_id}_TP{tp}"

        path = Path(dir_path) / f"plate_{repl+1}"
        path.mkdir(parents=True, exist_ok=True)
        with open(path / "epg.csv", "w") as f:
            pd.DataFrame(epg_file_content).to_csv(f, index=False)

        with open(path / "plate_map.csv", "w") as f:
            pd.DataFrame(plate_map).to_csv(f, index=False, header=False)
