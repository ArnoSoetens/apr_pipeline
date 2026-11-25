# APR pipeline (config-driven)

```
project/
├── apr_pipeline/          python modules (apr_config.py reads data/inputs/config.toml;
│                          datapoints.py knows the state of each data point; the setup stages
│                          structure.py, placement.py, parametrize.py, apr_setup.py also run on
│                          their own: python apr_pipeline/placement.py [names])
├── scripts/
│   ├── main.py            setup: structure -> placement -> parametrize -> APR windows, each point from
│   │                      the stage it needs (--from STAGE; --force redoes simulated points)
│   ├── analysis.py        free energies -> results.json, plots/, data/inputs/results.csv
│   ├── run.py             runs APR windows: picks them, writes data/inputs/*.in + run.env, starts run.sh
│   ├── run.sh             the worker: pmemd on the given window folders, N at a time (plain bash)
│   ├── job_gpu.slurm      one GPU job: main.py then run.py --cpus 1 (sbatch scripts/job_gpu.slurm [names])
│   ├── run_implicit.py    guest-only GB runs: --implicit-only / --p000 -> <dp>/windows_implicit/
│   ├── open_windows.sh    cpptraj + VMD:  ./open_windows.sh iso_00001 a000-a005
│   ├── submit_arrays.py   slurm arrays: status | setup | sim | analysis | all | implicit  (--submit to submit)
│   ├── job_array.slurm    shared by the arrays: modules, conda, amber, then array_task.sh
│   ├── array_task.sh      what one array task does (its slice of the worklist)
│   ├── collect_results.py builds results.csv from all results.json
│   ├── cleanup.py         removes transient files of analyzed points (--dry-run first)
│   ├── cleanup_batch_logs.py   removes a finished batch's slurm logs (started by submit_arrays)
│   ├── config_env.py      bash cannot read TOML: prints `export` lines for array_task.sh / open_windows.sh
│   └── templates/         the Amber .in templates with ${placeholders} (APR windows and guest-only)
└── data/                  the project: everything a campaign needs and produces
    ├── inputs/
    │   ├── config.toml    <- the one file you edit (tracked by git)
    │   ├── BCD.pdb / .mol2 / .frcmod   host files (tracked by git)
    │   ├── <list>.txt     SMILES list (input.smiles_list)
    │   └── *.in, results.csv          written by the scripts
    ├── logs/              slurm output and worklists of submitted arrays
    └── iso_00001/ ...     data points
```

* Needs Python >= 3.11 (tomllib), i.e. your `paprika` environment.
* Every script reads `data/inputs/config.toml`; the other folders follow from the
  data folder (there is no `[directories]` section any more). To work on another
  data folder, move it to `data/`, or set `APR_DATA=/path/to/that/folder` (it needs its
  own `inputs/config.toml`).
* Archiving or copying `data/` keeps every setting and input of the simulations.
* `scripts/` assumes `apr_pipeline/` sits next to it.

Starting a new campaign: put the SMILES list in `data/inputs/` and set
`smiles_list` under `[input]` in `data/inputs/config.toml`, then
`python scripts/submit_arrays.py all --submit` (or `python scripts/main.py` locally) creates
the data point folders and sets them up. For a GPU run set `run.pmemd = "cuda"`.

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
Placement writes a rejected guest pose to `UNL_clash.pdb` (`<guest_prefix>_clash.pdb`)
and leaves `UNL.pdb` and the complex PDB alone, so the input is still found next
time; a successful placement removes `UNL_clash.pdb`. For folders from before
this, `--failed-clashes` falls back to `UNL.pdb` when there is no complex PDB.
Only the host is initially visible; use VMD's molecule list to switch to a
guest molecule. Re-run setup/placement to regenerate guest poses after changing
the placement or clash-detection algorithm.

## Guest-only implicit-solvent simulation

Two guest-only GB runs per data point, next to the APR windows:

```sh
python scripts/run_implicit.py --implicit-only              # the free guest, every parametrized point
python scripts/run_implicit.py --p000 iso_00001 iso_00002   # the guest alone in the p000 window
python scripts/run_implicit.py --implicit-only --p000 ...   # both
python scripts/run_implicit.py ... --list todo.txt          # points from a file (also ranges)
python scripts/run_implicit.py ... --cpus 4                 # 4 data points at a time (default 1)
python scripts/run_implicit.py ... --clean                  # keep only .prmtop, production.nc, done marker
python scripts/run_implicit.py ... --dry-run                # only list what would run
```

| run | folder | what |
|---|---|---|
| `--implicit-only` | `<dp>/windows_implicit/implicit_only/` | the free, unrestrained guest: tleap (`guest.prmtop`), minimize, equilibrate, production with `data/inputs/implicit_only/*.in` (production continues equilibration with its velocities) |
| `--p000` | `<dp>/windows_implicit/implicit_p000/` | one host-free APR p000 window: the guest plus DM1/DM2/DM3, aligned and restrained with the guest's APR restraints at p000 (`disang.rest`, `restraints.json`), minimize, equilibrate, production with the APR window `.in` files in `data/inputs/` |

A point needs its `UNL.pdb`, `UNL.mol2` and `UNL.frcmod` (after parametrize); without them it
is skipped with the reason. The p000 restraints use the G1/G2 anchors from `guest_input.json`,
the same atoms as the complex (or the chain ends, lower one first, while that is empty).
The pmemd command is `run.pmemd` from the config. A run that finishes gets the done marker
and is never redone; an unfinished one starts again from scratch. A failed run gets
`run_failed.txt` (step, exit code, last line of its output) and the others go on; the exit
code is 1 if a run failed. Host-anchor and host-conformational restraints are omitted in
p000. The third guest restraint matches the current APR setup, where its three masks
define an angle in pAPRika despite the config field being named `torsion_target`.
This is **GB implicit solvent**, not a periodic box of explicit water.

`python scripts/cleanup.py` also cleans finished guest-only runs the way `--clean` does,
whether or not the point has been analyzed.

Not in the config on purpose: the `#SBATCH` lines (sbatch reads them before anything else runs),
the BCD-specific per-glucose restraint atom numbers in `apr_setup.set_restraints`,
and the output file names inside `run.sh` (production.nc, ...).

The optional attach-phase distance wall is configured under
`[apr.restraints.walls.distance]` and is disabled by default. It is one
DM1-G1 flat-well restraint: `r2` and `r3` set the zero-energy interval, while
`r1`/`r4` set the parabolic-to-linear transition points and `rk2`/`rk3` set the
two force constants. With `mode = "manual"`, the configured `r1`–`r4` values
are used directly. For example, `r2 = 1` and `r3 = 12` mean the flat region
spans DM1-G1 distances from 1 to 12 Å. Manual mode does not inspect the host
or require rim atom numbers. With `mode = "rim"`, setup fits planes to the
seven configured oxygen atoms on each host rim, verifies that they align with
the APR axis and bracket G1/G2, then uses the DM1-to-rim distances for `r2`
and `r3`; `parabolic_width` sets the shoulder widths (`r1` is clamped to zero
if needed). `rim_offset` expands the flat region: setting it to 3 Å makes `r2`
3 Å smaller and `r3` 3 Å larger than their fitted rim distances. Rim atom
groups must be supplied in the config when rim mode is selected.

`guest = "centroid"` (default `"g1"`) puts the wall on the distance from DM1 to the
centroid of all guest atoms instead of G1, as an Amber group restraint (`iat=…,-1` with
`igr2`). Amber uses the unweighted, coordinate-averaged position of a group, not the centre
of mass, so the force is shared equally over the guest atoms and does not turn the guest. In
rim mode the flat region is then near/far rim ∓ (guest radius + `radius_offset`), where the
guest radius is the largest centroid-to-heavy-atom distance in the setup structure:
at `radius_offset = 0` the wall starts where the whole guest has left the cavity along the
axis, and `-1` starts it 1 Å earlier. `rim_offset` is not used then. Setup refuses a guest
whose centroid starts outside the flat region. The same `guest` key under
`[apr.restraints.walls.cone]` puts the cone on the D2-D1-guest centroid angle (and the
start check on the centroid). The analysis measures a group wall from the trajectory with the
same unweighted centroid.

## Setup: serial or as an array

```
python scripts/main.py                       # every data point, each from the stage it still needs
python scripts/main.py iso_00001 iso_00002   # just these (what an array task does)
python scripts/main.py --from placement      # redo placement and everything after it (all, or the named ones)

python scripts/submit_arrays.py setup --submit        # creates data point folders from the configured SMILES list
                                                      # one array task = a few data points (slurm.setup_*)
python scripts/submit_arrays.py setup --retry-failed --submit
```

Setup has four stages, run in this order for each data point:

| stage | needs | writes |
|---|---|---|
| `structure` | `<name>.smi` | the 3D guest, `structure.output_pdb` (default `UNL.pdb`) |
| `placement` | one guest PDB | `UNL.pdb` + `BCD-UNL.pdb` (guest in the host cavity) |
| `parametrize` | `UNL.pdb` | `UNL.mol2` + `UNL.frcmod` |
| `apr_setup` | the four files above | `complex/` and `windows/` |

`main.py` looks at the files in each folder and starts at the first stage whose
output is missing (existing files are trusted), so a folder can start from a
SMILES, from any guest PDB, from a placed guest, or from all four prepared
files. `--from STAGE` starts there regardless of what exists. A data point that
fails gets `setup_failed.txt` (`<stage>: <reason>`) and the others go on.
A data point with simulation output (any `windows/*/production.nc`) is never
touched by setup, not even with `--from` or by name; only `main.py --force` redoes it
(see "Changing settings and redoing setup").

Choosing data points works the same in `main.py` and every stage script:

```
python scripts/main.py                          # every data point that needs it
python scripts/main.py iso_00001 PFBA_m_down    # these
python scripts/main.py iso_00001-iso_00010      # a range (only numbers that have a folder)
python scripts/main.py --list todo.txt          # names/ranges from a file, one or more per line
```

In a `--list` file, blank lines and `# comments` are ignored, and entries can be
mixed with names on the command line. A range needs the same prefix and the
same number of digits on both ends; numbers without a folder are listed in a
note.

Each stage can also be run on its own, with the same logging and
`setup_failed.txt` as `main.py`:

```
python apr_pipeline/placement.py              # every folder whose next stage is placement
python apr_pipeline/placement.py my_guest     # placement for data/my_guest, even if it was placed
```

The same works for `structure.py`, `parametrize.py` and `apr_setup.py`. By
name, only that stage is redone: if the folder already had later stages'
output, a note says so and gives the `main.py NAME --from STAGE` command that
redoes them too. A stage script exits with code 1 when a data point failed.
APR setup always removes an earlier `complex/` and `windows/` before building,
so no window from an older window list (e.g. a longer pull) is left behind for
`run.sh` to simulate. That is safe because setup never runs on a folder with a
trajectory.

The SMILES list (`input.smiles_list`, a file in `data/inputs/`) is only read when the data
directory has no data point folders at all.

For custom guest structures, create a folder under `data/` for each input and
put exactly one input PDB (with any filename) in each folder. Then run:

```
python apr_pipeline/placement.py             # every folder that still needs placement
python apr_pipeline/placement.py my_guest    # only data/my_guest/
```

Placement writes `UNL.pdb` and the host-guest complex PDB to the same folder as
the input. The configured host (`data/inputs/BCD.pdb`) and placement
settings in `data/inputs/config.toml` are used. Any `.pdb` in the folder is an input
except `BCD-UNL.pdb` and antechamber's `sqm.pdb`; `UNL.pdb` only counts when it
is the only one (next to another PDB it is an earlier placement's output).
Folders must contain exactly one input. The guest's residue name in the input
does not matter.

Instead of running the stage scripts one by one, `python scripts/main.py
my_guest` runs placement, parametrize and APR setup in one go for such a folder
(set its charge first, see `guest_input.json` below).

To parametrize guest PDBs separately, run:

```
python apr_pipeline/parametrize.py                       # every placed folder without UNL.mol2/.frcmod
python apr_pipeline/parametrize.py my_guest
```

Each data point has a `guest_input.json` that the setup stages hand to each other:

```json
{ "smiles": "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F", "net_charge": -1, "anchors": [1, 5] }
```

* `smiles` is written when the folder has a `.smi`; for a guest that starts as
  a PDB it stays `null`.
* `net_charge` is the charge parametrize uses, always taken from this file.
  While it is `null`, `parametrize.net_charge` in `data/inputs/config.toml` fills it in:
  * `"auto"`: the formal charge of the SMILES. A PDB-only guest stays `null`
    (a PDB alone does not reliably give the charge), and parametrize asks you
    to set it.
  * an integer, e.g. `-1`: that charge, for every data point. Handy for a list
    of PDB guests that all have the same charge.

  A charge already in the file is never replaced, so you can always set one
  data point by hand (e.g. the one neutral guest in a list of anions, or a
  SMILES charge that is wrong; then redo `--from parametrize`). To apply a
  changed config value to a point that already has a charge, set its
  `net_charge` back to `null`.
* `anchors` are G1 and G2, the guest atoms the APR restraints hold, as atom
  serial numbers in `UNL.pdb` (its second column). G1 is the end nearest the
  dummy atoms: G1→G2 points along `placement.cavity_axis`, the way the guest is
  pulled out. Folders prepared before this file existed get them recomputed from
  `UNL.pdb` (lower end first) and saved.
  * `null` (the default): placement takes the two ends of the longest carbon
    chain and turns the guest as little as possible, so a guest given upside
    down stays upside down. Whichever end then sits lower becomes G1, and both
    are written to the file.
  * Set by hand: placement uses them as given and turns the guest until G1→G2
    points along `cavity_axis`, however far that is. They are never
    overwritten. Either write them before the first run (the input PDB must then
    number its atoms 1 to N in order, otherwise placement refuses it), or place
    the guest once, read the numbers from `UNL.pdb`, edit the file and rerun with
    `main.py my_guest --from placement`. Set them back to `null` to return to the
    automatic choice, e.g. after replacing the input PDB with one in a different
    orientation.
  * APR setup never swaps them. If they don't match the placed guest (G1 above
    G2, e.g. after editing them and running only `--from apr_setup`), it stops
    and asks you to redo `--from placement`.

`structure.output_pdb` in `data/inputs/config.toml` names the file the 3D guest built from
the SMILES is written to. With the default `UNL.pdb`, placement overwrites it.

Create or edit the SMILES text file locally before submitting. By default, the array script reads
`data/inputs/<input.smiles_list>` (set in `data/inputs/config.toml`), which must contain one SMILES
per line. `main.py` and the setup array create data-point folders from that file only when the
data directory has no data point folders yet; otherwise the list file is not needed. Analysis only
runs when every configured window contains the configured production trajectory (normally
`production.nc`). Result collection includes analyzed prepared data points without a `.smi` file,
leaving the SMILES field empty. XML-to-text conversion must be done separately before running the
pipeline. The `setup` command then submits the usual per-data-point setup jobs.

* A data point that fails gets `setup_failed.txt` with the stage and reason (shown by `status`); success removes it.
  The `setup` array leaves failed data points out unless you pass `--retry-failed`.
* Wait for the setup array to finish before `sim`: the simulation list is made from the windows setup creates.

## Changing settings and redoing setup

Setup never touches a data point with simulation output, so after changing a setting, redo the
affected points with `--force`, starting at the earliest stage the setting changes:

| you change | redo with |
|---|---|
| walls, window list (`attach_string`, pull distances), restraints, dummy atoms, host anchors | `main.py --force --from apr_setup NAMES` |
| `net_charge` (in `guest_input.json` or the config), charge method, atom types | `main.py --force --from parametrize NAMES` |
| placement settings, anchors in `guest_input.json`, the host | `main.py --force --from placement NAMES` |
| structure settings | `main.py --force NAMES` (from the start: structure with a `.smi`, otherwise placement) |

`--force` first shows what it will delete and asks:

```
--force deletes the simulation output of 3 data point(s): 135 trajectories (2.1 GB), results of 2, guest-only runs of 0:
  iso_00001, iso_00002, iso_00003
Cancel running or queued Slurm jobs for these points first.
Continue? [y/N]
```

Only `y`/`yes` continues. Without a terminal (Slurm, a pipe) it refuses unless you add `--yes`.
For each point, just before its setup, it deletes `windows/`, `complex/`, `results.json`,
`plots/`, `analysis_failed.txt` and old plots in the folder itself, and from `parametrize` or
earlier also `windows_implicit/` (built from the old guest parameters). Inputs (`.smi`, input
PDBs, `guest_input.json`) and your own files are never deleted; the stages overwrite their own
outputs. `results.csv` is rebuilt afterwards. Then run the windows again (`run.py`, or
`submit_arrays.py sim --submit`). The array tasks never use `--force`.

## Running the APR windows

```
python scripts/run.py                          # every unfinished window of every set-up data point
python scripts/run.py iso_00001 iso_00002      # the unfinished windows of these
python scripts/run.py iso_00001/a000-a005      # these windows (also iso_00001/p010)
python scripts/run.py --list todo.txt          # names, ranges or windows from a file
python scripts/run.py --cpus 4 ...             # 4 CPUs; default run.max_jobs (0 = all CPUs); --cpus 1 = serial
python scripts/run.py --dry-run ...            # only list what would run
```

`run.py` picks the windows, writes the three Amber `.in` files and `run.env` (the pmemd
command and names `run.sh` needs) into `data/inputs/` from the config, and starts
`scripts/run.sh` on the window folders. `run.sh` is plain bash (no Python): it runs
minimize, production part 1 and part 2 per window, with as many windows at the same time as
there are CPUs (CPUs / `run.mpi_ranks` with `run.pmemd = "mpi"`); a free CPU takes the next
window, also from the next data point. A window that is finished (`production2.done` and
`production.nc`) is never run again, so no trajectory is ever overwritten; a window without
topology is skipped; a data point whose setup is not finished is skipped with the reason, and so is
an analyzed one (its windows no longer look finished after cleanup) unless you name its windows
(`run.py iso_00001/a000` runs that window again).
A window whose pmemd step fails gets `run_failed.txt` (date, host, step, exit code, last line
of that step's output); the other windows go on, and a later successful run removes the file.
The summary line counts done, failed and skipped windows.

GPU: set `pmemd = "cuda"` in `data/inputs/config.toml` and submit from the repository root

```
sbatch scripts/job_gpu.slurm                       # setup + every unfinished window
sbatch scripts/job_gpu.slurm iso_00001/a000-a005   # setup of iso_00001 + these windows
```

It runs `main.py` and then `run.py --cpus 1` with the same arguments (`main.py` accepts window
entries and sets up their data point), on one GPU; its log is `data/logs/gpu_<job id>.out`.

## Slurm arrays

```
python scripts/submit_arrays.py status            # what is finished, what is left, per stage
python scripts/submit_arrays.py all               # shows setup -> sim -> analysis -> collection commands
python scripts/submit_arrays.py all --submit      # submits setup -> simulation -> analysis -> results.csv
python scripts/submit_arrays.py all --batches 6 --submit
python scripts/submit_arrays.py all --batches 6 --clean-after-analysis --submit
python scripts/submit_arrays.py sim --windows-per-task 8 --submit    # smaller packages for a small job
python scripts/submit_arrays.py implicit --implicit-only --p000 --submit   # guest-only runs
```

Each array task takes its lines of the worklist (task k: lines k*N+1 ... k*N+N) and hands them
to the same script you would run locally: setup -> `main.py <names>`, analysis ->
`analysis.py --no-collect <names>`, implicit -> `run_implicit.py --cpus <slurm CPUs> <flags> <names>`,
simulation -> `run.sh` on the window folders. The scripts skip what is not ready and say why.
A simulation task runs **no Python at all**: `submit_arrays.py` writes the Amber `.in` files and
`data/inputs/run.env` (data folder, pmemd command, ...) when submitting, and the task only uses
bash and pmemd, so a slow shared file system can't break Python start-up in the most numerous tasks.

* `status` sorts the data points by the setup stage they still need (the same rule
  `main.py` uses), lists the ones whose setup failed with the stage and reason, then
  counts the set-up ones by simulation/analysis state, for example:

  ```
  72 data points in .../data
    setup (python scripts/main.py, or 'setup --submit')
         2  need placement    (a guest PDB, not placed yet)   e.g. PFBA_m_up, PFBS_mt_up
        24  failed last time  (setup_failed.txt; 'setup --retry-failed')   e.g. iso_00003, ...
                iso_00003: placement: clash: 2 host-guest atom pair(s) closer than 1 A ...
    simulation and analysis
        30  windows ready, simulation not finished   e.g. iso_00001, ...
        16  analysis done   e.g. iso_00010, ...
    windows still to simulate: 912
  ```

  `setup` submits exactly the points listed under "setup" (failed ones only with
  `--retry-failed`). A `setup_failed.txt` in a point whose setup is complete is ignored.
  Below that, only when there are any: windows that failed last time (`run_failed.txt`, with
  step and exit code; `sim` simply runs them again), points whose analysis failed
  (`analysis_failed.txt`, with the reason) and the guest-only runs done/failed per kind.
* Every window that still has to be simulated goes into one list (`iso_00001/a000`, ...).
  By default, the script calculates windows per CPU from `slurm.sim_task_hours` and
  `slurm.sim_minutes_per_window`, reserving 25% of the task time. For example, 4 hours and
  5 minutes per window gives 36 windows per CPU. One array task runs one package on `sim_cpus`
  CPUs, one window per CPU at a time; a package may span several data points. Use
  `--windows-per-task` to override the calculated package size.
* More tasks than `max_array_size` are split over several arrays (each with an offset into the list).
* A finished window gets a `production2.done` file; resubmitting only puts unfinished windows in the list.
  A window counts as finished only with that marker **and** `production.nc`: an interrupted window
  (partial trajectory, no marker) is simulated again and is never analyzed.
* Analysis tasks make the extra plots when `analysis.enhanced = true` in `data/inputs/config.toml`.
* Analysis defaults to one data point per task and `analysis_mem` of `4G` in `data/inputs/config.toml`.
  This avoids retaining PyMBAR/JAX memory across several analyses in one Python process.
  Override the task size with `--points-per-task` or adjust `analysis_mem` for the cluster.
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
* `implicit` submits `run_implicit.py` for parametrized points without the requested finished runs
  (`--implicit-only`, `--p000` or both; `--clean` keeps only what is needed). Task size and resources:
  `slurm.implicit_data_points_per_task`, `implicit_cpus` (points at a time), `implicit_minutes_per_data_point`,
  `implicit_mem`.
* `--clean-after-analysis` is opt-in. After a data point analyzes successfully, its transient files
  are removed, with the same rules as `scripts/cleanup.py` (below). With `--clean-after-analysis` in
  `all`, a finalizer removes that batch's simulation and analysis array output logs only after those
  arrays succeed and every batch point is confirmed analyzed. If any array fails or a point remains
  incomplete, those logs are retained. Worklists and logs from setup, other batches, or non-batch
  jobs are not removed.
* To clean analyzed data points by hand:

  ```
  python scripts/cleanup.py --dry-run            # what would be removed, per point, with sizes
  python scripts/cleanup.py --dry-run -v iso_00001   # the same, listing every file
  python scripts/cleanup.py                      # every point with results.json
  python scripts/cleanup.py iso_00001-iso_00010  # these;  --list todo.txt  points from a file
  ```

  Only points with `results.json` are touched. In the data point folder only setup leftovers are
  removed (`complex/`, `sqm.*`, antechamber scratch files); everything else stays: guest inputs,
  `guest_input.json`, input PDBs, `plots/`, `results.json` and any file of your own. In `windows/`
  everything is removed except `restraints.json`, the trajectories of `a000`-`a004` and the last
  attach window (for example `a015`) and the topology in `a000`; the window folders themselves
  stay. `open_windows.sh` uses the retained `a000` topology when another selected window has no
  local topology.
* Windows simulated before the markers existed: if you know they are complete,
  `for w in data/*/windows/*/; do touch "${w}production2.done"; done` marks them
  (a plain `touch data/*/windows/*/production2.done` does NOT work: the `*` only matches files that already exist).

## Analysis results and diagnostics

```
python scripts/analysis.py                      # every fully simulated point without results.json
python scripts/analysis.py iso_00001 PFBA_m_up  # these (again, if they were analyzed before)
python scripts/analysis.py iso_00001-iso_00010  # a range;  --list todo.txt reads names from a file
python scripts/analysis.py --enhanced ...       # with the convergence plot, attach grid, pull histograms
python scripts/analysis.py --no-enhanced ...    # without them
```

The default for both is `analysis.enhanced` in `data/inputs/config.toml` (`true`: with the extras).

Each data point gets `results.json` and `plots/binding_free_energy.jpg`. Enhanced adds `plots/convergence_plot.jpg`, `plots/attach_windows_grid.jpg` (DM1-G1 distance per
attach window as a mirrored histogram: 1st half of the frames up, 2nd half down, rims dotted with
the region between them shaded), `plots/attach_angle_grid.jpg` (DM2-DM1-G1 angle, same style, cone
wall dotted with the inside of the cone shaded; the legend says when the cone acts on the guest
centroid rather than the plotted G1), `plots/dihedral_overview.jpg` (guest torsions, see below),
`plots/pull_window_histograms.jpg` and some extra metadata in `results.json`
(`method`, `attach_fractions`, `pull_distances`, `overlap`, `monotonicity_attach`). Without enhanced the free
energy is only computed for the full trajectories, not for the ten fractions
the convergence plot needs, so the analysis is faster. The calculation is the
same: MBAR, the reference-state work and the wall/threshold diagnostics come out
identical. The TI-block values (and so `dG_bind`/`dG_sem` with
`method = "ti"`) are bootstrap estimates, which differ by a few hundredths of a
kcal/mol between any two runs, with or without `--enhanced`.

A named point that is not fully simulated, or that was cleaned up after an
earlier analysis (trajectories removed), is skipped with the reason. A point
whose analysis fails gets `analysis_failed.txt` (date and reason) and the loop
goes on with the next one; a later successful analysis removes the file. The
exit code stays 0, so read the summary line. At the end `results.csv` is
rebuilt from every `results.json` in the data folder, including cleaned-up
points: `data/inputs/results.csv` by default (`analysis.results_csv`, relative
to the data folder). `--no-collect` skips that (the Slurm analysis tasks use it;
a separate job collects at the end). `--clean-after-analysis` removes transient
files after each successful analysis.

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

With `--enhanced`, a data point also gets `plots/convergence_plot.jpg`: one row for attach and one
for pull. The left panel shows the TI-block free energy from the first ("forward") and from the
last ("reverse") 10%, 20%, ..., 100% of every window's frames, with ±1 SEM. The right panel shows
that SEM against the number of frames used on log-log axes, with an N^-1/2 line through the
full-data point. Every SEM is a new blocking analysis of the frames used (pAPRika itself only
rescales the full-data SEM, which would follow N^-1/2 by construction); dU, the blocking and the
bootstrap/Akima integration are pAPRika's, so at 100% the value equals the `ti-block` result up to
bootstrap noise. Forward and reverse should meet well before 100%; a gap at small fractions points
to a not yet equilibrated start (forward) or a late drift (reverse). MBAR is left out: pAPRika's
`mbar-autoc` SEM does not correct for autocorrelation (it passes the frame count instead of the
time series to `detect_equilibration`, so g = 1), so its error bars are too small. It is a
diagnostic, not a proof that the simulation has sampled all relevant states.

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

When distance walls are present in the attachment windows' `disang.rest`
files, `results.json` also contains `distance_wall`. For the configured
single-wall setup, it contains that wall's summary directly; for multiple
walls, it contains a `walls` list with one summary per restraint pair. A group wall (guest
centroid) has `-1` in `amber_atom_indices_1based` and its atoms in `amber_groups_1based`. The
analysis reconstructs the restraint potential from each sampled wall distance and its Amber `r1`–`r4` and
`rk2`/`rk3` parameters. In the parabolic regions the Amber
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

With `--enhanced`, `overlap` holds one number per phase (`attach`, `pull`): the smallest
MBAR overlap-matrix element between neighbouring windows (pymbar `compute_overlap`, all frames,
the same reduced potentials as pAPRika's MBAR). Identical neighbours give 0.5; below about 0.03 a
pair overlaps poorly and the windows there are too far apart.

`monotonicity_attach` checks that the mean restraint energy <R> (full force constants, the TI
integrand) does not rise with lambda: d<R>/dlambda = -beta Var(R) <= 0, i.e. the attach work curve
is concave. For each neighbouring pair z = (rise) / sqrt(SEM_i^2 + SEM_j^2) with pAPRika's blocking
SEM; a pair is a violation when z > z_crit = Phi^-1(1 - 0.05/(K-1)) (one-sided, Bonferroni over the
K-1 pairs; 2.61 for 12 windows). `monotonic` is false when there is any violation; `violations`
lists their lambda pair, rise (kcal/mol) and z.

`plots/dihedral_overview.jpg` (with `--enhanced`) has two heatmaps with every window folder (a000 ...
p0NN) as a row and every rotatable guest torsion as a column. The torsions come from RDKit's
non-ring torsion list of the guest, built from the bonds in the window topology (no atom types or
bond orders needed; a symmetric rotor such as CH3-CH3 has none, and then no plot is made). Left:
the number of transitions between the g+/t/g- rotamer states (cores of ±30° around 60/180/300°; a
frame between cores keeps its last state, so jitter at a boundary does not count), coloured from 0
to 10 (≥ 10 shown as 10) with the count in each cell. Right: the Jensen-Shannon divergence (bits,
10° bins) between the first and second half of the window: 0 = same distribution, 1 = no overlap.
A torsion with few transitions and a high JS was in one state for the first half and another for
the second. Note that left- and right-handed perfluoroalkyl helices (about ±165°) both fall in the
t core, so a helix reversal is not a transition; a carboxylate O-C-C-C torsion (two equivalent O)
hops often without anything physically changing.


## Dihedral distributions

Compare rotatable guest dihedrals between the unrestrained `a000` and restrained
`a015` windows by giving just the data-point name:

```sh
python other/dihedral_check.py iso_00001
```

If no dihedrals are supplied, the script reads the guest mol2 from the
datapoint folder (using `naming.guest_prefix` in `data/inputs/config.toml`) and
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
