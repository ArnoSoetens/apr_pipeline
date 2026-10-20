# APR pipeline (config-driven)

```
project/
├── config.toml            <- the one file you edit
├── apr_pipeline/          python modules (apr_config.py reads config.toml;
│                          datapoints.py knows the state of each data point)
├── scripts/
│   ├── main.py            setup: isomers -> structure -> placement -> parametrize -> APR windows
│   │                      (serial: python main.py;  or per data point name, as array tasks do)
│   ├── analysis.py        free energies -> results.csv
│   ├── run.sh             runs all windows (pmemd choice + .in values come from config.toml)
│   ├── open_windows.sh    cpptraj + VMD:  ./open_windows.sh iso_00001 a000-a005
│   ├── job1.slurm / job1_gpu.slurm      single jobs (one node, all data points)
│   ├── submit_arrays.py   slurm arrays: status | sim | analysis | all  (--submit to really submit)
│   ├── job_array.slurm    shared by the arrays: modules, conda, amber, then array_task.sh
│   ├── array_task.sh      what one array task does (its slice of the worklist)
│   ├── collect_results.py builds results.csv from all results.json
│   ├── config_env.py      bash cannot read TOML: prints `export` lines for run.sh / open_windows.sh
│   ├── make_inputs.py     templates/*.in  ->  inputs/*.in  (filled from config.toml)
│   └── templates/         the three Amber .in files with ${placeholders}
└── data/                  BCD.pdb / .mol2 / .frcmod, data_input.xml, iso_00001/, ...
```

* Needs Python >= 3.11 (tomllib), i.e. your `paprika` environment.
* `APR_CONFIG=/path/to/other.toml` makes everything use a different config.
* A variable exported in the environment beats the config for bash settings
  (`export PMEMD=pmemd`, `export MAX_JOBS=4`).
* `scripts/` assumes `apr_pipeline/` sits next to it.

Not in the config on purpose: the `#SBATCH` lines (sbatch reads them before anything else runs),
the BCD-specific per-glucose restraint atom numbers in `apr_setup.set_restraints`,
and the output file names inside `run.sh` (production.nc, ...).

## Setup: serial or as an array

```
python scripts/main.py                   # serial: folders from the XML, then every data point that needs setup
python scripts/main.py --force           # ... also redo the complete ones
python scripts/main.py iso_00001 iso_00002   # just these (what an array task does)

python scripts/submit_arrays.py setup --submit        # creates data point folders from the XML on first run
                                                      # one array task = a few data points (slurm.setup_*)
python scripts/submit_arrays.py setup --retry-failed --submit
```

* A data point that fails gets `setup_failed.txt` with the reason (shown by `status`); success removes it.
  The `setup` array leaves failed data points out unless you pass `--retry-failed`.
* Wait for the setup array to finish before `sim`: the simulation list is made from the windows setup creates.

## Slurm arrays

```
python scripts/submit_arrays.py status            # what is finished, what is left
python scripts/submit_arrays.py all               # shows the sbatch commands, submits nothing
python scripts/submit_arrays.py all --submit      # simulation arrays -> analysis array -> results.csv
python scripts/submit_arrays.py sim --windows-per-task 8 --submit    # smaller packages for a small job
```

* Every window that still has to be simulated goes into one list (`iso_00001/a000`, ...).
  The list is cut into work packages of `slurm.sim_cpus x slurm.sim_windows_per_cpu` windows
  (6 CPUs x 45 = 270; 4 CPUs x 45 = 180). One array task runs one package on `sim_cpus` CPUs,
  one window per CPU at a time; a package may span several data points.
* More tasks than `max_array_size` are split over several arrays (each with an offset into the list).
* A finished window gets a `production2.done` file; resubmitting only puts unfinished windows in the list.
* For a small job (few tasks) the script suggests a smaller `--windows-per-task`, because then all
  tasks start at once and the job finishes much sooner.
* The analysis array starts when all simulation arrays have ended (`afterany`) and skips data points
  whose windows are not all finished. `collect_results.py` then writes `results.csv`.
* Windows simulated before the markers existed: if you know they are complete,
  `for w in data/*/windows/*/; do touch "${w}production2.done"; done` marks them
  (a plain `touch data/*/windows/*/production2.done` does NOT work: the `*` only matches files that already exist).
* `run.sh` takes data point names (all their windows) or `data_point/window` pairs, and keeps
  `MAX_JOBS` windows running over data point borders; without arguments it runs every window it finds.
