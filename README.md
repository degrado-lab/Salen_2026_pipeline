# SDSP: Sequence Design / Structure Prediction pipeline

An iterative enzyme active-site design pipeline that alternates fixed-backbone
sequence design with flexible-backbone relaxation, then validates the resulting
designs with two independent structure predictors.

```
LigandMPNN  <->  Rosetta FastRelax        (3 rounds)
                        |
                        v
              RaptorX-Single              (single chain, protein only)
                        |
                        v
                     Chai-1               (ligand-aware, holo + apo)
                        |
                        v
                  analysis / filtering
```

Each design round expands the sequence pool: round 1 produces `init_num_seqs`
sequences, round 2 produces `init_num_seqs^2`, round 3 `init_num_seqs^3`.
After the final relaxation, designs are scored (pLDDT, Cα RMSD, net charge,
ligand ΔSASA, Rosetta total score and packstat), filtered, and the survivors
are passed to Chai-1 for ligand-aware prediction.

## Important: this pipeline is cluster-specific

The pipeline is driven from a Jupyter notebook that writes shell scripts and
submits them to an **SGE / Grid Engine** scheduler with `qsub`. It was developed
and run on UCSF's Wynton cluster, and the job scripts contain site-specific
resource requests (`-q gpu.q`, `-l gpu_mem=8G`, `-l hostname=qb3-atgpu*`,
`module load cuda`).

To run it elsewhere you will need to adapt the job headers in the notebook
cells and in `test/scripts/chai1_scripts/prep_submit_analyze_chai1.py`.
On a Slurm cluster the scheduler calls (`qsub`, `qstat`, `$SGE_TASK_ID`,
`$SGE_GPU`) need replacing with their Slurm equivalents. The scripts under
`test/scripts/` are plain Python and are scheduler-independent.

## Requirements

### Third-party software (installed separately)

| Tool | Used for | Source |
| --- | --- | --- |
| LigandMPNN | ligand-aware sequence design | https://github.com/dauparas/LigandMPNN |
| Rosetta | FastRelax, packstat scoring, ligand `.params` generation | https://rosettacommons.org |
| RaptorX-Single | protein-only structure prediction | https://github.com/AndersJing/RaptorX-Single |
| Chai-1 (`chai_lab`) | ligand-aware structure prediction | https://github.com/chaidiscovery/chai-lab |

You also need the LigandMPNN side-chain packing checkpoint
(`model_params/ligandmpnn_sc_v_32_002_16.pt`) and the RaptorX-Single
parameter file (`params/RaptorX-Single-ESM1b-ESM1v-ProtTrans.pt`), both
obtained through their respective repositories.

`download.py` fetches the Chai-1 model components, cached RDKit conformers,
and the ESM-2 tokenizer/model into the `chai_lab` download path. Run it once,
from an environment with `chai_lab` installed, before the first Chai-1 job.

### Python packages

Scripts import `numpy`, `pandas`, `biopython` (`Bio.PDB`), `torch`,
`chai_lab` and `transformers`. Package versions were not pinned during
development; if you need an exact environment, generate one from your
working install (`conda env export > environment.yml`).

### Conda environments

The job scripts activate environments **by name**. Either create environments
with these names or edit the `conda activate` lines:

| Environment | Activated by |
| --- | --- |
| `ligandmpnn_env` | LigandMPNN and Rosetta FastRelax rounds |
| `RaptorX-Single` | RaptorX-Single prediction |
| `lmpnn_rs_ipynbs` | analysis and filtering scripts |
| `Chai1` | Chai-1 inference |
| `chai1_analysis` | Chai-1 output analysis |

Note that `conda activate` inside a non-interactive batch job usually requires
sourcing conda's shell hook first, e.g.
`source "$(conda info --base)/etc/profile.d/conda.sh"`. Add this to the job
headers if your cluster does not configure it globally.

## Directory layout

A run directory must have the structure below. `test/` is a worked example
you can copy as a template — its `input/` holds a small test case and its
`output/` subdirectories are empty placeholders that the pipeline fills.

In the example input, note that the Rosetta ligand (`sal.params`) carries a
`ZN1` atom typed as `Zn2p` and the constraints in `constraints.cst` are written
against it, while the `chai1_ligand_smiles` setting in the notebook describes
the corresponding cobalt complex. If you adapt this example, make sure your
`.params` atom names match the atom names used in your `.cst` file.

```
<run_directory>/
├── input/
│   ├── 1_lmpnn_in/
│   │   └── test.pdb                     starting structure
│   ├── constraints/
│   │   └── constraints.cst              Rosetta-format constraints (filename is fixed)
│   ├── omit_aa_per_residue/
│   │   └── omit_aa_per_residue.json     per-position LigandMPNN restrictions (filename is fixed)
│   └── params/
│       ├── sal.params                   Rosetta ligand params
│       └── sal_0001.pdb                 ligand conformer
├── output/
│   ├── 1_rs_in/                         LigandMPNN round 1 output
│   ├── 2_lmpnn_in/                      FastRelax round 1 output
│   ├── 2_rs_in/                         LigandMPNN round 2 output
│   ├── 3_lmpnn_in/                      FastRelax round 2 output
│   ├── 3_rs_in/                         LigandMPNN round 3 output
│   ├── final_out/                       FastRelax round 3 output
│   ├── rxs_predictions/
│   │   └── fastas/                      RaptorX-Single input fastas
│   ├── chai1_predictions/
│   │   ├── input_fastas/                Chai-1 input fastas (holo + apo)
│   │   └── output/                      one subdirectory per Chai-1 prediction
│   ├── analysis/                        analysis and filtered CSVs
│   ├── ipynb_job_scripts/               generated SGE job scripts
│   └── wynton/                          SGE stdout/stderr
└── scripts/
    ├── lmpnn_rsfr_rxs_scripts/
    │   ├── create_pdb_json.py           builds the multi-PDB JSON for LigandMPNN
    │   ├── fast_relax_packstat_csts.xml RosettaScripts protocol (round 3)
    │   ├── filter_pdbs.py               applies pLDDT / RMSD / charge cutoffs
    │   ├── find_pdb_filename.py         resolves a PDB path from an array-task tag
    │   ├── lmpnn_rsfr_rxs_analysis.py   per-design metrics -> CSV
    │   ├── mpnn_out_to_raptor_x_input.py  LigandMPNN .fa -> one fasta per sequence
    │   ├── mpnn_out_to_single_file_fasta.py  LigandMPNN .fa -> single combined fasta
    │   ├── pdb_tag.py                   appends numeric array tags to PDB filenames
    │   ├── remove_packed_tag.py         strips LigandMPNN "_packed_N_M" tags
    │   ├── remove_tag.py                strips a trailing "_tag" (standalone; unused by the pipeline)
    │   └── rename_raptor_x_single.py    strips the RaptorX model suffix
    └── chai1_scripts/
        ├── analyze_chai1.py             collects Chai-1 metrics -> CSV
        ├── prep_submit_analyze_chai1.py preps fastas, submits and analyses Chai-1
        └── helpers/
            ├── analyze_chai1_helpers.py
            ├── find_fasta_file.py
            ├── prep_and_submit_chai1_helpers.py
            └── run_chai1.py
```

## Running the pipeline

1. Copy `test/` to a new run directory, or recreate the layout above.
2. Place your starting structure in `input/1_lmpnn_in/`, your Rosetta
   constraints in `input/constraints/constraints.cst`, your ligand params in
   `input/params/`, and your per-residue restrictions in
   `input/omit_aa_per_residue/omit_aa_per_residue.json`.
   Generate the ligand `.params` file with Rosetta's `mol2_to_params.py`.
3. Open `SDSP_lmpnn_rsfr_raptorx_chai1.ipynb` and edit **the first code cell
   only**. It holds every run-specific setting: `run_id`, `input_pdb`,
   `init_num_seqs`, `params_file`, `pipeline_run_dir`, the LigandMPNN
   restriction flags, the Chai-1 ligand name and SMILES, the analysis cutoffs,
   and the paths to your LigandMPNN / Rosetta / RaptorX-Single installations.
   The installation paths are placeholders (`/path/to/...`) and must be set
   before the first run.
4. Run the cells in order. Each cell writes a job script into
   `output/ipynb_job_scripts/` and submits it with a dependency on the
   previous job, so the whole chain can be queued in one pass. Monitor with
   `qstat`.

Results land in `output/analysis/` as
`<run_id>_lmpnn_rsfr_rxs_analysis.csv` (all designs),
`<run_id>_lmpnn_rsfr_rxs_analysis_filtered.csv` (designs passing the cutoffs)
and `chai1_analysis_output.csv` (Chai-1 metrics, tab-separated).

## Known limitations

- The Chai-1 preparation step processes at most 100 designs per run
  (`MAX_DESIGNS` in `helpers/prep_and_submit_chai1_helpers.py`).
- `lmpnn_rsfr_rxs_analysis.py` and `filter_pdbs.py` write their CSV output to
  the current working directory; the notebook `cd`s into `output/analysis`
  before calling them.
- RMSD comparison requires the model and the prediction to have the same
  number of Cα atoms. Mismatches are logged to `mismatched_comparisons.txt`
  and that design is skipped.
- In the Chai-1 analysis, a prediction whose CIF or NPZ fails to parse
  contributes `NaN` to the corresponding mean rather than being dropped.
- The notebook's round-1 LigandMPNN cell passes
  `--omit_AA_per_residue "./inputs/omit_AA_per_residue.json"`, a path relative
  to the LigandMPNN repository, rather than the `lmpnn_omit_aa_per_resi`
  variable defined in the settings cell. Rounds 2 and 3 do not pass the flag
  at all. This is the behaviour the published results were generated with and
  has been left unchanged; set the flag explicitly if you want per-residue
  restrictions applied in every round.

## Citation

If you use this pipeline, please cite the accompanying paper.

> Citation details to be added on publication.

## License

MIT. See [LICENSE](LICENSE).

## Credits

Developed in the [DeGrado Lab](https://github.com/degrado-lab), UCSF.

Pipeline scripts written by Ian Bakanas. Analysis script and adaptation to the
Jupyter notebook by Stephanie Ouchida.

Third-party methods are the work of their respective authors; please cite
LigandMPNN, Rosetta, RaptorX-Single and Chai-1 alongside this repository.
