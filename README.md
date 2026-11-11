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
└── data/                  BCD.pdb / .mol2 / .frcmod, SMILES .txt, iso_00001/, ...
```

* Needs Python >= 3.11 (tomllib), i.e. your `paprika` environment.
* `APR_CONFIG=/path/to/other.toml` makes everything use a different config.
* A variable exported in the environment beats the config for bash settings
  (`export PMEMD=pmemd`, `export MAX_JOBS=4`).
* `scripts/` assumes `apr_pipeline/` sits next to it.

`scripts/open_windows.sh` also supports static structure inspection:

```
./scripts/open_windows.sh --complex-only iso_00001 - iso_00010
./scripts/open_windows.sh --complex-only
./scripts/open_windows.sh --failed-clashes
./scripts/open_windows.sh --failed-clashes iso_00001-iso_00010
```

`--complex-only` opens the saved host+guest complex PDBs for the selected data
points. With no data points supplied, it opens every configured complex PDB in
the data directory. Data-point ranges accept either `iso_00001-iso_00010` or
`iso_00001 - iso_00010`. `--failed-clashes` scans `setup_failed.txt` files for
clash reports, loads the configured host PDB once and each failed guest PDB as
a separate molecule, and highlights the nearest guest atom with a Licorice
representation.
Failed guest poses are written to the normal guest PDB (`UNL.pdb` by default).
Only the host is initially visible; use VMD's molecule list to switch to a
guest molecule. Re-run setup/placement to regenerate guest poses after changing
the placement or clash-detection algorithm.

## Guest-only implicit-solvent simulation

To simulate a guest separately from the host, use a guest-only PDB and its
parameterized MOL2:

```sh
python scripts/simulate_guest.py guest.pdb guest.mol2
python scripts/simulate_guest.py guest.pdb guest.mol2 \
  --frcmod guest.frcmod --output-dir guest_run
```

The script builds an Amber topology with LEaP and runs minimization,
equilibration, and production using the force field, GB model, temperature,
timestep, thermostat, and run lengths configured in `config.toml`. If
`--frcmod` is omitted, `parmchk2` generates one from the MOL2. The PDB must
contain one guest residue whose residue name and atom names match the MOL2;
the MOL2 must already contain the desired GAFF atom types and partial charges.

This is **GB implicit solvent**, not a periodic box of explicit water. The
simulation is unrestrained and writes its topology, input files, logs, restart
files, and trajectories to the output directory (by default
`<PDB-stem>_implicit`). AmberTools (`parmchk2`, `tleap`) and the configured
Amber `pmemd` executable must be available in the environment.

To run every immediate subdirectory of the configured data directory that has
the configured guest PDB, MOL2, and FRCMOD files, use the batch runner.
Directory names do not need to use the `iso_` prefix:

```sh
python scripts/simulate_guest_batch.py --cpus-per-task 24
python scripts/simulate_guest_batch.py --simulation apr-p000 --cpus-per-task 24
```

In a Slurm job, it defaults to `SLURM_CPUS_PER_TASK`, so the job script can run
`python "$SLURM_SUBMIT_DIR/scripts/simulate_guest_batch.py" --simulation apr-p000`
after activating the Python and Amber environments to select the p000 workflow.
Each data point writes into its own output directory (`guest_implicit/` for the
default mode or `guest_apr_p000/` for p000), uses serial `pmemd` to keep each
concurrent simulation to one CPU, and records command output in
`batch.log`. Directories missing any of the three guest files are reported and
skipped. Simulation failures are reported without stopping other data points;
after all jobs finish, the batch command exits nonzero if any simulation failed.

The `apr-p000` mode instead runs one host-free APR pull window per data point:
the guest plus DM1/DM2/DM3, aligned and restrained using pAPRika and the APR
settings in `config.toml`. It writes to `guest_apr_p000/`. Use
`--simulation apr-p000` in a Slurm batch job to select this mode. Host-anchor
and host-conformational restraints are omitted. The third guest restraint
matches the current APR setup, where its three masks define an angle in
pAPRika despite the config field being named `torsion_target`.

Not in the config on purpose: the `#SBATCH` lines (sbatch reads them before anything else runs),
the BCD-specific per-glucose restraint atom numbers in `apr_setup.set_restraints`,
and the output file names inside `run.sh` (production.nc, ...).

An optional attach-phase distance wall is configured under
`[apr.restraints.walls.distance]`. Set `enabled = true` to turn it on and
configure `r1`–`r4` in Angstroms. `r2`–`r3` is the flat, zero-energy region;
the regions outside it are parabolic between the bounds and linear beyond
`r1`/`r4`. It is disabled by default.

## Setup: serial or as an array

```
python scripts/main.py                   # serial: folders from the configured SMILES list, then every data point that needs setup
python scripts/main.py --force           # ... also redo the complete ones
python scripts/main.py iso_00001 iso_00002   # just these (what an array task does)

python scripts/submit_arrays.py setup --submit        # creates data point folders from the configured SMILES list
                                                      # one array task = a few data points (slurm.setup_*)
python scripts/submit_arrays.py setup --retry-failed --submit
```

For custom guest structures, create a folder under `data/` for each input and
put exactly one input PDB (with any filename) in each folder. Then run:

```
python apr_pipeline/placement.py             # process all data-point folders
python apr_pipeline/placement.py my_guest    # process only data/my_guest/
```

Placement writes `UNL.pdb` and the host-guest complex PDB to the same folder as
the input. The configured host (`data/BCD.pdb` by default) and placement
settings in `config.toml` are used. Existing generated `UNL.pdb` and
`BCD-UNL.pdb` files are ignored when finding the input; folders must contain
exactly one other PDB.

To run APR setup on a datapoint whose guest files are already prepared, place
`UNL.pdb`, `UNL.mol2`, `UNL.frcmod`, and `BCD-UNL.pdb` (using the configured
guest and complex prefixes) in its data folder, then run:

```
python scripts/main.py <data-point>
```

`main.py` uses those files as-is and runs setup without regenerating or
parametrizing them. Such folders are also discovered when `main.py` is run
without names. Folders without a `.smi` file and without all four prepared
files are skipped.

To parametrize guest PDBs separately, run:

```
python apr_pipeline/parametrize.py                       # all folders containing UNL.pdb
python apr_pipeline/parametrize.py my_guest --net-charge 0
```

When `parametrize.net_charge = "auto"` in `config.toml`, standalone
parametrization requires `--net-charge` because a PDB file alone does not
reliably provide the guest's formal charge. `scripts/main.py` continues to
derive the charge from each SMILES molecule and pass it to the same
parametrization function.

Create or edit the SMILES text file locally before submitting. By default, the array script reads
`data/Seven_PFAS.txt` (`directories.smiles_list` in `config.toml`), which must contain one SMILES
per line. Array setup creates data-point folders from that file only when no `.smi` or fully
prepared data-point folders exist. Serial `main.py` can use existing `.smi` folders without the
list file, and it does not regenerate guest inputs when fully prepared data points are present.
APR setup and analysis discover both SMILES-based and fully prepared data points. Analysis only
runs when every configured window contains the configured production trajectory (normally
`production.nc`). Result collection includes analyzed prepared data points without a `.smi` file,
leaving the SMILES field empty. XML-to-text conversion must be done separately before running the
pipeline. The `setup` command then submits the usual per-data-point setup jobs.

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
  By default, the script calculates windows per CPU from `slurm.sim_task_hours` and
  `slurm.sim_minutes_per_window`, reserving 25% of the task time. For example, 4 hours and
  5 minutes per window gives 36 windows per CPU. One array task runs one package on `sim_cpus`
  CPUs, one window per CPU at a time; a package may span several data points. Use
  `--windows-per-task` to override the calculated package size.
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
  plots, trajectories from `a000`-`a004` and the last attach window (for example `a015`), the
  topology in `a000`, and `windows/restraints.json`. Other files in window folders are removed,
  but existing window directories are kept. Log files within a data point are removed, along with
  other directories left empty by cleanup. Files are deleted only within that data point; the shared
  configured Slurm logs directory is otherwise untouched. With `--clean-after-analysis` in `all`,
  a finalizer removes that batch's simulation and analysis array output logs only after those arrays
  succeed and every batch point is confirmed analyzed. If any array fails or a point remains
  incomplete, those logs are retained. Worklists and logs from setup, other batches, or non-batch
  jobs are not removed. `open_windows.sh` uses the retained `a000` topology when another selected
  window has no local topology.
  To clean already-analyzed data points manually, run
  `python scripts/cleanup.py --data-dir ../data`. It applies the same retention rules and skips
  points without `results.json`; add data-point names to limit cleanup to selected folders.
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

When a two-atom distance wall is present in the attachment windows'
`disang.rest` files, `results.json` also contains `distance_wall`. The analysis
reconstructs the restraint potential from each sampled wall distance and its
Amber `r1`–`r4` and `rk2`/`rk3` parameters. In the parabolic regions the Amber
energy is `rk * (distance - boundary)²` (without a one-half factor); beyond
`r1`/`r4`, the potential continues linearly with the boundary slope. The output
reports the mean and maximum instantaneous potential energy, the fraction of
frames outside the flat region, and signed net wall work from the first to last
saved frame **for each independent window**. Positive work means the wall did
net work on the simulated system; negative work means the system did net work
against the wall. Because the restraint parameters are fixed during a window,
this net work is `U(first frame) - U(last frame)`. Per-window work values are
not summed across separate simulations. Mean potential energy is a restraint
energy statistic, not work or a binding free-energy contribution. No
`distance_wall` section is added when no such restraint is present.

`analysis.enhanced` continues to control the supplementary histogram plots and
coordinate metadata in `results.json`. The convergence plot and phase estimates
are generated regardless of that setting.

## Restraint-angle histograms

`paprika.analysis.fe_calc.collect_data()` reads trajectories to report the
coordinates of restraints in `windows/restraints.json`. This helper plots the
first angle restraint's sampled coordinate by attachment window, using the
atom indices and target stored in that JSON rather than assuming a particular
guest atom name. It uses MDTraj to read the trajectories.

```sh
python scripts/angle_histogram.py
```

By default this reads `data/PFBA/windows/{a000,a001,a002,a003,a004,a015}/production.nc`
and plots the first angle restraint, `DM2:DUM - DM1:DUM - UNL:C1`, with its
configured target shown as a dashed line. It writes the per-window histogram
and frame-by-frame CSV to `data/PFBA/analysis/`. Use `--windows`, `--data-dir`,
`--output-dir`, or `--bins` to override those defaults.

## Dihedral distributions

Compare rotatable guest dihedrals between the unrestrained `a000` and restrained
`a015` windows by giving just the data-point name:

```sh
python other/dihedral_check.py iso_00001
```

If no dihedrals are supplied, the script reads the guest mol2 from the
datapoint folder (using `naming.guest_prefix` in `config.toml`) and
automatically finds its rotatable torsions. It calls RDKit's
`TorsionFingerprints.CalculateTorsionLists`, takes its non-ring torsions
(using RDKit's rotatable-bond rules), then maps each torsion's mol2 atom names
to the matching guest-residue atoms in the Amber topology. GAFF atom-type
labels are normalized for RDKit parsing; atom names and bond records are
preserved. The topology is read from `a000`; trajectories and the data
directory follow `config.toml`. By default, `a000` is compared with `a015`;
pass `--windows a003 a010` to choose another pair:

```sh
python other/dihedral_check.py iso_00013 --windows a003 a010 \
  -first 149,150,151,152
```

You can still provide atom quadruplets explicitly with `-first`, `-second`,
`-third`, positional arguments, `--dihedral`, or `--dihedral-file`. Automatic
detection runs by default only when none are specified; `--guest-mol2` can be
used to select another mol2 file for detection.

The first selected window is treated as unrestrained and the second as
restrained. The topology continues to come from `a000`. Results go to
`data/<datapoint>/dihedral_out/`, including
`histogram_counts.csv` (counts per angle bin), count-based histogram plots,
summary statistics, and the timeseries CSV. A single shared timeseries figure
contains a subplot for each dihedral; its filename includes all selected atom
sets so different selections do not overwrite each other's plots. The existing
`--dihedral`, `--dihedral-file`, and explicit trajectory options remain
available.

## Acknowledgments and references

Isomer enumeration uses `enu` from [CombiFF](https://github.com/csms-ethz/CombiFF).
The APR setup, simulation, and analysis code uses
[pAPRika](https://github.com/GilsonLabUCSD/pAPRika).

Please cite these projects as follows:

- Oliveira, M. P., Andrey, M., Rieder, S. R., Kern, L., Hahn, D. F., Riniker,
  S., Horta, B. A. C., and Hünenberger, P. H. “Systematic Optimization of a
  Fragment-Based Force Field against Experimental Pure-Liquid Properties
  Considering Large Compound Families: Application to Saturated Haloalkanes.”
  *Journal of Chemical Theory and Computation* 16(12), 7525–7555 (2020).
  [https://doi.org/10.1021/acs.jctc.0c00683](https://doi.org/10.1021/acs.jctc.0c00683)
- Slochower, D., Setiadi, J., Henriksen, N., Boothroyd, S., Thompson, M., and
  Wang, L. *pAPRika*, version 1.1.0 (2021).
  [https://github.com/GilsonLabUCSD/pAPRika](https://github.com/GilsonLabUCSD/pAPRika).
  This follows the citation metadata in pAPRika's
  [CITATION.cff](https://github.com/GilsonLabUCSD/pAPRika/blob/master/CITATION.cff).
