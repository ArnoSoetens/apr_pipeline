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
python scripts/submit_arrays.py all               # shows setup -> sim -> analysis -> collection commands
python scripts/submit_arrays.py all --submit      # submits setup -> simulation -> analysis -> results.csv
python scripts/submit_arrays.py all --batches 6 --submit
python scripts/submit_arrays.py all --batches 6 --clean-after-analysis --submit
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
* Later stages wait until prerequisite jobs have ended (`afterany`), whether they succeeded or failed.
  This lets unrelated data points and later batches continue after a failure. Analysis skips points
  whose windows are not all finished; failed analysis points keep their files and are omitted from
  `results.csv`.
* `all` includes data points that need setup, and submits the later stages with Slurm dependencies.
  Previously failed setup points are left out unless `--retry-failed` is passed.
* `all --batches N` splits the automatically generated data-point worklist into up to N balanced,
  non-empty batches. One setup array handles all points; each batch then runs its simulation array(s)
  followed by its analysis array. The next batch's simulations wait for the previous batch's analysis
  tasks to end, but still run if an analysis task fails.
  `N=1` keeps the single-batch workflow. All windows belonging to one data point stay in the same batch.
* `--clean-after-analysis` is opt-in. After a data point analyzes successfully, transient files are
  removed while preserving its configured guest inputs, complex PDB, SMILES, results JSON, analysis
  plots, the first five window trajectories, `a000`'s topology, and `windows/restraints.json`.
  Files are deleted only within that data point's directory; directories are kept. `open_windows.sh`
  uses the retained `a000` topology when another selected window has no local topology.
* Windows simulated before the markers existed: if you know they are complete,
  `for w in data/*/windows/*/; do touch "${w}production2.done"; done` marks them
  (a plain `touch data/*/windows/*/production2.done` does NOT work: the `*` only matches files that already exist).
* `run.sh` takes data point names (all their windows) or `data_point/window` pairs, and keeps
  `MAX_JOBS` windows running over data point borders; without arguments it runs every window it finds.

## Analysis results and diagnostics

Every analysis run computes both estimators. `results.json` always contains
`phases.attach` and `phases.pull`. Each phase contains one `ti-block` and one
`mbar-autoc` result, each with a single total free-energy estimate (`fe`) and
Paprika standard error (`sem`), in kcal/mol. No per-window free energies are
duplicated into `results.json`. TI uses Paprika's block-based TI analysis; MBAR
uses its autocorrelation-based sampling treatment.
The top-level `dG_bind` and `dG_sem` remain the values from the configured
`analysis.method`, so existing summary CSVs and comparisons retain their previous
meaning. The configured method's plot profile is still used internally to draw
`binding_free_energy.jpg`, but its per-window values are not repeated in JSON.

Every data point also gets a `convergence_plot.jpg`, with separate attach and
pull panels. Each panel shows the TI and MBAR phase free-energy estimates at
10%, 20%, ..., 100% of each window's trajectory frames, with Paprika's
fraction-specific SEM as error bars. Each fraction uses the first indicated
portion of the available frames in every window. This lets you see whether the
estimate settles toward its full-trajectory value; it is a diagnostic, not a
proof that the simulation has sampled all relevant states. Set
`analysis.convergence_plot = false` in `config.toml` to skip creating this plot;
it defaults to `true`.

`attach_distance_outside_threshold` contains the fraction and percent of
attachment distance samples satisfying
`abs(sample - initial_distance) > 7 Å`. Values exactly 7 Å from the initial
distance are not counted as outside. The attachment phase shares its fully
attached endpoint with the first pulling window (`p000`), so this calculation
uses only attachment simulation windows `a000` through `a015` and excludes the
shared `p000` samples. The JSON lists the included windows and sample count.
Because trajectory frames are saved at a fixed interval, the sample fraction
is also the sampled-time fraction. Non-finite samples cause an explicit error
instead of being silently discarded.

`analysis.enhanced` continues to control the supplementary histogram plots and
coordinate metadata in `results.json`. The convergence plot and phase estimates
are generated regardless of that setting.
