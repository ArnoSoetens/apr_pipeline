"""Remove completed simulation and analysis array logs for a successful batch."""
import argparse
import re
import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent / "apr_pipeline"))

from apr_config import CFG
from datapoints import DATA, scan


def cleanup_batch_logs(log_dir, batch_number, simulation_ids, analysis_id, point_names):
    log_dir = Path(log_dir)
    if not log_dir.is_dir():
        raise FileNotFoundError(f"Configured Slurm log directory does not exist: {log_dir}")

    incomplete = [
        name for name in point_names if scan(DATA / name)[0] != "analyzed"
    ]
    if incomplete:
        print(
            f"preserving batch {batch_number:02d} logs; "
            f"{len(incomplete)} data point(s) are not analyzed: {', '.join(incomplete[:5])}"
        )
        return

    job_ids = [("analysis", analysis_id)]
    job_ids.extend(("sim", job_id) for job_id in simulation_ids)
    removed = 0
    for stage, job_id in job_ids:
        if not re.fullmatch(r"\d+", job_id):
            raise ValueError(f"Invalid Slurm job ID: {job_id}")
        prefix = f"{stage}_b{batch_number:02d}_{job_id}_"
        for path in log_dir.glob(f"{prefix}*.out"):
            if path.is_file() and not path.is_symlink():
                path.unlink()
                removed += 1

    print(f"removed {removed} Slurm output log(s) for batch {batch_number:02d}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_number", type=int)
    parser.add_argument("simulation_ids", help="comma-separated Slurm IDs, or '-'")
    parser.add_argument("analysis_id")
    parser.add_argument("point_names", nargs="+")
    args = parser.parse_args()
    simulation_ids = [] if args.simulation_ids == "-" else args.simulation_ids.split(",")
    cleanup_batch_logs(
        CFG["directories"]["logs"],
        args.batch_number,
        simulation_ids,
        args.analysis_id,
        args.point_names,
    )


if __name__ == "__main__":
    main()
