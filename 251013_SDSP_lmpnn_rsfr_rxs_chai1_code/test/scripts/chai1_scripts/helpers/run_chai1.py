from pathlib import Path
import numpy as np
import torch
from chai_lab.chai1 import run_inference
import sys

def run_chai1(fasta_path,output_dir): 
    torch.cuda.empty_cache()
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    output_cif_paths = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        # 'default' setup
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=42,
        device=torch.device("cuda:0"),
        use_esm_embeddings=True,
    )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_chai1.py <chai1_run_input_fasta_path> <chai1_run_root_output_dir_path>")
        sys.exit(1)

    chai1_run_input_fasta_path = sys.argv[1]
    chai1_run_output_path = sys.argv[2]

    run_chai1(chai1_run_input_fasta_path, chai1_run_output_path)