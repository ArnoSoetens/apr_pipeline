#!/usr/bin/env python
# Setup APR windows for BCD-FOA, orientation 1, implicit solvent (GB).
# Produces complex/ and windows/ which are renamed by run.sh.

import os
import subprocess
from pathlib import Path

import numpy as np
import parmed as pmd

from paprika import restraints
from paprika.build import align, dummy
from paprika.build.system import TLeap
from paprika.io import save_restraints
from paprika.restraints.amber import amber_restraint_line
from paprika.restraints.utils import create_window_list

import shutil

from apr_config import CFG, get_attach_fractions, get_pull_distances, get_radii_set

# --- values from config.toml ---
HOST = CFG["naming"]["host_prefix"]
GUEST = CFG["naming"]["guest_prefix"]
COMPLEX = CFG["naming"]["complex_prefix"]
ATOM_TYPE = CFG["parametrize"]["atom_type"]
RADII = get_radii_set()

INITIAL_DISTANCE = CFG["apr"]["initial_distance"]
FINAL_DISTANCE = CFG["apr"]["final_distance"]
DUMMY_Z = CFG["apr"]["dummy_z"]
DUMMY3_Y = CFG["apr"]["dummy3_y"]
HOST_ANCHOR_ATOMS = CFG["apr"]["host_anchor_atoms"]

DIST_FC = CFG["apr"]["restraints"]["distance_fc"]
ANGLE_FC = CFG["apr"]["restraints"]["angle_fc"]
TORSION_FC = CFG["apr"]["restraints"]["torsion_fc"]
ANGLE_TARGET = CFG["apr"]["restraints"]["angle_target"]
TORSION_TARGET = CFG["apr"]["restraints"]["torsion_target"]
STATIC_FC = CFG["apr"]["restraints"]["static_fc"]
CONF_FC = CFG["apr"]["restraints"]["host_conformational_fc"]
WALLS = CFG["apr"]["restraints"]["walls"]


##host = "BCD"
##guest = "UNL"


# --- Build vacuum complex with tleap ---

def build_system( data_dir_host, data_dir_guest, host_prefix=HOST, guest_prefix=GUEST, radii_set=RADII):
    os.makedirs("complex", exist_ok=True)
    system = TLeap()
    system.output_path = "complex"
    system.pbc_type = None
    system.neutralize = False
    system.template_lines = [
        f"source leaprc.{ATOM_TYPE}",
        "loadamberparams frcmod.ions1lm_126_tip3p",
        f"set default PBRadii {radii_set}",
        f"loadamberparams {data_dir_host}/{host_prefix}.frcmod",
        f"{host_prefix} = loadmol2 {data_dir_host}/{host_prefix}.mol2",
        f"loadamberparams {data_dir_guest}/{guest_prefix}.frcmod",
        f"{guest_prefix} = loadmol2 {data_dir_guest}/{guest_prefix}.mol2",
        f"model = loadpdb {data_dir_guest}/{COMPLEX}.pdb",
        "check model",
        "savepdb model vac.pdb",
        "saveamberparm model vac.prmtop vac.rst7",
    ]
    system.build(clean_files=False)
    return

# --- Align the pulling axis and add dummy atoms ---
# Orientation 1: C3 leads (faces away from host opening during pull).

def guest_masks(b1_idx, b2_idx, guest_prefix=GUEST):
    """Amber masks of the G1/G2 anchors (0-based guest atom indices), in the given order.

    Placement turns G1->G2 onto +cavity_axis, so G1 is the end nearest the dummy atoms. Raises when
    the placed guest has them the other way round (e.g. anchors edited after placement), instead
    of swapping them: the anchors in guest_input.json are used exactly as written.
    """
    from placement import order_anchors

    s = pmd.load_file("complex/vac.prmtop", "complex/vac.rst7", structure=True)
    guest_atoms = [a for a in s.atoms if a.residue.name == guest_prefix]
    if not (0 <= b1_idx < len(guest_atoms) and 0 <= b2_idx < len(guest_atoms)):
        raise IndexError(
            f"backbone atom indices ({b1_idx}, {b2_idx}) exceed the "
            f"{len(guest_atoms)} atoms in residue {guest_prefix}"
        )
    coords = [(a.xx, a.xy, a.xz) for a in guest_atoms]
    if order_anchors(coords, b1_idx, b2_idx) != (b1_idx, b2_idx):
        raise ValueError(
            f"anchors [{b1_idx + 1}, {b2_idx + 1}] don't match the placed guest: G1 (the first) must "
            "be below G2 along placement.cavity_axis; redo from placement (--from placement) so the "
            "guest is placed along these anchors"
        )
    G1 = f":{guest_prefix}@{guest_atoms[b1_idx].idx + 1}"
    G2 = f":{guest_prefix}@{guest_atoms[b2_idx].idx + 1}"
    return G1, G2


def validate_dummy_anchor_geometry(structure, G1, G2):
    dm1_atoms = [atom for atom in structure.atoms if atom.residue.name == "DM1"]
    if len(dm1_atoms) != 1:
        raise ValueError(f"expected exactly one DM1 atom, found {len(dm1_atoms)}")

    g1 = structure[G1].coordinates[0]
    g2 = structure[G2].coordinates[0]
    dm1 = np.array([dm1_atoms[0].xx, dm1_atoms[0].xy, dm1_atoms[0].xz])
    g1_distance = float(np.linalg.norm(g1 - dm1))
    g2_distance = float(np.linalg.norm(g2 - dm1))
    if g1_distance >= g2_distance:
        raise ValueError(
            "DM1 is not closer to G1 than G2 after alignment "
            f"(G1: {g1_distance:.3f} A, G2: {g2_distance:.3f} A)"
        )


def add_dummy(G1, G2, data_dir_host, data_dir_guest, host_prefix=HOST, guest_prefix=GUEST, radii_set=RADII):
    structure = pmd.load_file("complex/vac.prmtop", "complex/vac.rst7", structure=True)
    aligned = align.zalign(structure, G1, G2)
    aligned.save("complex/aligned.prmtop", overwrite=True)
    aligned.save("complex/aligned.rst7", overwrite=True)

    structure = pmd.load_file("complex/aligned.prmtop", "complex/aligned.rst7", structure=True)
    structure = dummy.add_dummy(structure, residue_name="DM1", z=DUMMY_Z[0])
    structure = dummy.add_dummy(structure, residue_name="DM2", z=DUMMY_Z[1])
    structure = dummy.add_dummy(structure, residue_name="DM3", z=DUMMY_Z[2], y=DUMMY3_Y)
    validate_dummy_anchor_geometry(structure, G1, G2)
    structure.save("complex/aligned.prmtop", overwrite=True)
    structure.save("complex/aligned.rst7", overwrite=True)
    structure.save("complex/aligned_with_dummy.prmtop", overwrite=True)
    structure.save("complex/aligned_with_dummy.rst7", overwrite=True)
    structure.save("complex/aligned_with_dummy.pdb", overwrite=True)

    dummy.write_dummy_frcmod(filepath="complex/dummy.frcmod")
    dummy.write_dummy_mol2(residue_name="DM1", filepath="complex/dm1.mol2")
    dummy.write_dummy_mol2(residue_name="DM2", filepath="complex/dm2.mol2")
    dummy.write_dummy_mol2(residue_name="DM3", filepath="complex/dm3.mol2")

# --- Final vacuum system with dummy atoms ---
    system = TLeap()
    system.output_path = "complex"
    system.pbc_type = None
    system.neutralize = False
    system.template_lines = [
        f"source leaprc.{ATOM_TYPE}",
        "loadamberparams frcmod.ions1lm_126_tip3p",
        f"set default PBRadii {radii_set}",
        f"loadamberparams {data_dir_host}/{host_prefix}.frcmod",
        f"{host_prefix} = loadmol2 {data_dir_host}/{host_prefix}.mol2",
        f"loadamberparams {data_dir_guest}/{guest_prefix}.frcmod",
        f"{guest_prefix} = loadmol2 {data_dir_guest}/{guest_prefix}.mol2",
        "loadamberparams dummy.frcmod",
        "DM1 = loadmol2 dm1.mol2",
        "DM2 = loadmol2 dm2.mol2",
        "DM3 = loadmol2 dm3.mol2",
        "model = loadpdb aligned_with_dummy.pdb",
        "check model",
        f"savepdb model {COMPLEX}-dum.pdb",
        f"saveamberparm model {COMPLEX}-dum.prmtop {COMPLEX}-dum.rst7",
    ]
    system.build()
    return

def window_setup():
    attach_fractions = get_attach_fractions()
    pull_distances = get_pull_distances()
    release_fractions = []
    windows = [len(attach_fractions), len(pull_distances), len(release_fractions)]
    return windows, attach_fractions


def _rim_oxygen_coordinates(structure, atom_numbers, host_prefix):
    """Coordinates of the configured rim oxygens (1-based host atom numbers)."""
    host_atoms = [atom for atom in structure.atoms if atom.residue.name == host_prefix]
    if any(number < 1 or number > len(host_atoms) for number in atom_numbers):
        raise ValueError(
            f"rim atom numbers {atom_numbers} must be within 1-{len(host_atoms)} "
            f"for host residue {host_prefix}"
        )

    atoms = [host_atoms[number - 1] for number in atom_numbers]
    if any(atom.atomic_number != 8 for atom in atoms):
        raise ValueError(
            f"configured rim atoms {atom_numbers} must all be oxygen atoms "
            f"in host residue {host_prefix}"
        )

    return np.array([[atom.xx, atom.xy, atom.xz] for atom in atoms])


def _rim_plane_distance_from_dm1(structure, atom_numbers, dm1_position, host_prefix):
    coordinates = _rim_oxygen_coordinates(structure, atom_numbers, host_prefix)
    center = coordinates.mean(axis=0)
    _, singular_values, right_vectors = np.linalg.svd(coordinates - center)
    if len(singular_values) < 2 or singular_values[1] < 1e-6:
        raise ValueError(f"configured rim atoms {atom_numbers} do not define a plane")

    normal = right_vectors[-1]
    axis = np.array([0.0, 0.0, 1.0])
    if np.dot(normal, axis) < 0:
        normal = -normal
    alignment = float(np.dot(normal, axis))
    rms_deviation = float(np.sqrt(np.mean(((coordinates - center) @ normal) ** 2)))
    if alignment < np.cos(np.deg2rad(15.0)) or rms_deviation > 0.25:
        raise ValueError(
            f"configured rim atoms {atom_numbers} do not form a rim plane "
            f"perpendicular to the APR axis (normal alignment {alignment:.3f}, "
            f"plane RMS {rms_deviation:.3f} A)"
        )

    axial_normal = float(np.dot(normal, axis))
    return float(np.dot(normal, center - dm1_position) / axial_normal)


def _create_distance_wall(
    structure, D1, guest_mask, attach_windows, bounds, rk2, rk3
):
    wall = restraints.DAT_restraint()
    wall.auto_apr = False
    wall.amber_index = True
    wall.topology = structure
    wall.mask1 = D1
    wall.mask2 = guest_mask

    wall.attach["fc_initial"] = 100.0
    wall.attach["fc_final"] = 100.0
    wall.attach["target"] = INITIAL_DISTANCE
    wall.attach["num_windows"] = attach_windows
    wall.custom_restraint_values.update(
        {
            "r1": bounds[0],
            "r2": bounds[1],
            "r3": bounds[2],
            "r4": bounds[3],
            "rk2": rk2,
            "rk3": rk3,
        }
    )
    wall.initialize()
    return wall


def _guest_centroid(structure, guest_prefix=GUEST):
    """Unweighted centroid of all guest atoms (what an Amber group restraint uses: the
    coordinate-averaged position, not the centre of mass) and the guest radius: the largest
    centroid-to-heavy-atom distance."""
    atoms = [atom for atom in structure.atoms if atom.residue.name == guest_prefix]
    if not atoms:
        raise ValueError(f"no atoms in guest residue {guest_prefix}")
    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in atoms])
    centroid = xyz.mean(axis=0)
    heavy = xyz[[atom.atomic_number > 1 for atom in atoms]]
    radius = float(np.max(np.linalg.norm((heavy if len(heavy) else xyz) - centroid, axis=1)))
    return centroid, radius


def create_wall_restraints(structure, D1, G1, G2, attach_windows, host_prefix=HOST,
                           guest_prefix=GUEST):
    """Flat-bottom distance wall from DM1 to G1 (guest = "g1") or to the centroid of all guest
    atoms (guest = "centroid", an Amber group restraint).

    Rim mode: r2/r3 = near/far rim plane distance from DM1 -/+ rim_offset (g1), or
    -/+ (guest radius + radius_offset) (centroid; the radius is the largest centroid-to-heavy-atom
    distance, so at radius_offset = 0 the wall starts where the whole guest has left the cavity
    along the axis). r1/r4 lie parabolic_width beyond r2/r3.
    """
    wall_config = WALLS["distance"]
    if not wall_config["enabled"]:
        return []
    on_centroid = wall_config.get("guest", "g1") == "centroid"

    dm1_position = structure[D1].coordinates[0]
    g1_position = structure[G1].coordinates[0]
    g2_position = structure[G2].coordinates[0]
    axis = np.array([0.0, 0.0, 1.0])
    g1_distance = float(np.dot(g1_position - dm1_position, axis))
    g2_distance = float(np.dot(g2_position - dm1_position, axis))
    if not 0 < g1_distance < g2_distance:
        raise ValueError(
            f"restraint atoms must lie in +z order from DM1 "
            f"(G1={g1_distance:.3f} A, G2={g2_distance:.3f} A)"
        )
    if on_centroid:
        centroid, guest_radius = _guest_centroid(structure, guest_prefix)
        start_distance = float(np.linalg.norm(centroid - dm1_position))

    if wall_config["mode"] == "rim":
        measured_rim_distances = sorted(
            _rim_plane_distance_from_dm1(
                structure, atom_group, dm1_position, host_prefix
            )
            for atom_group in wall_config["rim_atom_groups"]
        )
        near_rim, far_rim = measured_rim_distances
        if on_centroid:
            radius_offset = wall_config["radius_offset"]
            reach = guest_radius + radius_offset
            r2, r3 = near_rim - reach, far_rim + reach
            if r2 <= 0:
                raise ValueError(
                    f"guest radius {guest_radius:.3f} A + radius_offset {radius_offset:.3f} A "
                    f"must be smaller than the DM1-to-near-rim distance {near_rim:.3f} A"
                )
            if not r2 < start_distance < r3:
                raise ValueError(
                    "guest centroid starts outside the distance wall "
                    f"(DM1-centroid = {start_distance:.3f} A, flat region "
                    f"{r2:.3f}-{r3:.3f} A: near rim {near_rim:.3f} A, far rim {far_rim:.3f} A, "
                    f"guest radius {guest_radius:.3f} A, radius_offset {radius_offset:.3f} A)"
                )
        else:
            rim_offset = wall_config["rim_offset"]
            lower_bound = near_rim - rim_offset
            upper_bound = far_rim + rim_offset
            if not lower_bound < g1_distance < g2_distance < upper_bound:
                raise ValueError(
                    "configured beta-CD rim planes do not bracket the restraint atoms "
                    "from DM1 toward +z "
                    f"(near rim={near_rim:.3f} A, G1={g1_distance:.3f} A, "
                    f"G2={g2_distance:.3f} A, far rim={far_rim:.3f} A, "
                    f"rim_offset={rim_offset:.3f} A; allowed range "
                    f"{lower_bound:.3f}-{upper_bound:.3f} A); "
                    "check apr.restraints.walls.distance.rim_atom_groups"
                )
            r2, r3 = lower_bound, upper_bound
            if r2 <= 0:
                raise ValueError(
                    f"rim_offset {rim_offset:.3f} A must be smaller than the "
                    f"DM1-to-near-rim distance {near_rim:.3f} A"
                )
        parabolic_width = wall_config["parabolic_width"]
        bounds = (
            max(0.0, r2 - parabolic_width),
            r2,
            r3,
            r3 + parabolic_width,
        )
    else:
        bounds = tuple(wall_config[key] for key in ("r1", "r2", "r3", "r4"))
    if not bounds[0] < bounds[1] < bounds[2] < bounds[3]:
        raise ValueError(
            "wall bounds must satisfy r1 < r2 < r3 < r4 "
            f"(received {bounds})"
        )

    if on_centroid:
        label = "DM1-guest centroid"
        detail = (f" Guest radius {guest_radius:.3f} A; centroid starts "
                  f"{start_distance:.3f} A from DM1.")
    else:
        label, detail = "DM1-G1", ""
    print(
        f"{label} flat-well wall: "
        f"r1={bounds[0]:.3f}, r2={bounds[1]:.3f}, "
        f"r3={bounds[2]:.3f}, r4={bounds[3]:.3f} A.{detail}"
    )
    return [
        _create_distance_wall(
            structure,
            D1,
            f":{guest_prefix}" if on_centroid else G1,
            attach_windows,
            bounds,
            wall_config["rk2"],
            wall_config["rk3"],
        )
    ]


def _angle_degrees(a, b, c):
    """Angle a-b-c at b, in degrees."""
    u, v = np.asarray(a) - np.asarray(b), np.asarray(c) - np.asarray(b)
    cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def create_cone_wall(structure, D1, D2, G1, attach_windows, host_prefix=HOST, guest_prefix=GUEST):
    """Flat-bottom angle wall D2-D1-G1 (guest = "g1") or D2-D1-guest centroid (guest =
    "centroid", an Amber group restraint on all guest atoms): a cone with its tip at DM1 through
    the nearest rim.

    No force from 180 degrees down to r2 = (mean D2-D1-O angle of the oxygens of the rim plane
    nearest DM1) - angle_offset; harmonic (rk2, kcal/mol/rad^2) from r2 down to
    r1 = r2 - parabolic_width, linear below. Also prints how wide the cone is at the far rim.
    """
    cone = WALLS["cone"]
    if not cone["enabled"]:
        return []
    on_centroid = cone.get("guest", "g1") == "centroid"

    dm1 = structure[D1].coordinates[0]
    dm2 = structure[D2].coordinates[0]
    if on_centroid:
        guest_point, guest_label = _guest_centroid(structure, guest_prefix)[0], "guest centroid"
    else:
        guest_point, guest_label = structure[G1].coordinates[0], "G1"
    rims = []
    for group in WALLS["distance"]["rim_atom_groups"]:
        oxygens = _rim_oxygen_coordinates(structure, group, host_prefix)
        height = _rim_plane_distance_from_dm1(structure, group, dm1, host_prefix)
        angle = float(np.mean([_angle_degrees(dm2, dm1, o) for o in oxygens]))
        offsets = oxygens - oxygens.mean(axis=0)
        radius = float(np.mean(np.linalg.norm(offsets[:, :2], axis=1)))  # perpendicular to z
        rims.append((height, angle, radius))
    (near_h, near_angle, near_r), (far_h, far_angle, far_r) = sorted(rims)

    r2 = near_angle - cone["angle_offset"]
    r1 = r2 - cone["parabolic_width"]
    if not 0.0 < r1 < r2 < 180.0:
        raise ValueError(
            f"cone wall bounds out of range (r1={r1:.2f}, r2={r2:.2f} degrees; near-rim angle "
            f"{near_angle:.2f}, angle_offset {cone['angle_offset']}, parabolic_width "
            f"{cone['parabolic_width']}); check apr.restraints.walls.cone"
        )
    start_angle = _angle_degrees(dm2, dm1, guest_point)
    if start_angle < r2:
        raise ValueError(
            f"{guest_label} starts outside the cone wall "
            f"(D2-D1-{guest_label} = {start_angle:.2f} < r2 = {r2:.2f} degrees)"
        )

    half_angle = 180.0 - r2
    cone_r_far = far_h * np.tan(np.radians(half_angle))
    print(
        f"D2-D1-{guest_label} cone wall: r1={r1:.2f}, r2={r2:.2f}, r3=r4=180 deg "
        f"(half-angle {half_angle:.2f}; {guest_label} starts at {start_angle:.2f}); "
        f"near rim {near_h:.2f} A above DM1, radius {near_r:.2f} A, angle {near_angle:.2f}; "
        f"far rim {far_h:.2f} A, radius {far_r:.2f} A, angle {far_angle:.2f} "
        f"(half-angle {180.0 - far_angle:.2f}); cone radius at far rim {cone_r_far:.2f} A "
        f"= {cone_r_far / far_r:.2f}x the far ring."
    )

    wall = restraints.DAT_restraint()
    wall.auto_apr = False
    wall.amber_index = True
    wall.topology = structure
    wall.mask1 = D2
    wall.mask2 = D1
    wall.mask3 = f":{guest_prefix}" if on_centroid else G1
    wall.attach["fc_initial"] = cone["rk2"]
    wall.attach["fc_final"] = cone["rk2"]
    wall.attach["target"] = 180.0
    wall.attach["num_windows"] = attach_windows
    wall.custom_restraint_values.update(
        {"r1": r1, "r2": r2, "r3": 180.0, "r4": 180.0, "rk2": cone["rk2"], "rk3": cone["rk2"]}
    )
    wall.initialize()
    return [wall]


def set_restraints(G1, G2, host_prefix=HOST, guest_prefix=GUEST):

    windows, attach_fractions= window_setup()

    H = [f":{host_prefix}@{n}" for n in HOST_ANCHOR_ATOMS]
    H1, H2 = [], []
    for i in range(7):
        if i == 6:
            H1.append([f":{host_prefix}@142", f":{host_prefix}@127", f":{host_prefix}@129", f":{host_prefix}@12"])
            H2.append([f":{host_prefix}@127", f":{host_prefix}@129", f":{host_prefix}@12", f":{host_prefix}@14"])
        else:
            H1.append([f":{host_prefix}@{16+i*21}", f":{host_prefix}@{1+i*21}", f":{host_prefix}@{3+i*21}", f":{host_prefix}@{33+i*21}"])
            H2.append([f":{host_prefix}@{1+i*21}", f":{host_prefix}@{3+i*21}", f":{host_prefix}@{33+i*21}", f":{host_prefix}@{35+i*21}"])

    D1, D2, D3 = ":DM1", ":DM2", ":DM3"

    structure = pmd.load_file(f"complex/{COMPLEX}-dum.prmtop", f"complex/{COMPLEX}-dum.rst7")

    # Static restraints keep host in a fixed position/orientation throughout all windows.
    static_restraints = []
    for masks, fc in [
        ([D1, H[0]], STATIC_FC[0]),
        ([D2, D1, H[0]], STATIC_FC[1]),
        ([D3, D2, D1, H[0]], STATIC_FC[2]),
        ([D1, H[0], H[1]], STATIC_FC[3]),
        ([D2, D1, H[0], H[1]], STATIC_FC[4]),
        ([D1, H[0], H[1], H[2]], STATIC_FC[5]),
    ]:
        static_restraints.append(
            restraints.static_DAT_restraint(
                restraint_mask_list=masks,
                num_window_list=windows,
                ref_structure=structure,
                force_constant=fc,
                amber_index=True,
            )
        )
    for masks in H1 + H2:
        static_restraints.append(
            restraints.static_DAT_restraint(
                restraint_mask_list=masks,
                num_window_list=windows,
                ref_structure=structure,
                force_constant=CONF_FC,
                amber_index=True,
            )
        )

    # Guest restraints drive the APR coordinate: attach raises force constants,
    # pull translates the guest, release lowers force constants analytically.
    guest_restraints = []
    for mask1, mask2, mask3, attach_tgt, attach_fc, pull_tgt in [
        (D1,  G1,   None,  INITIAL_DISTANCE, DIST_FC, FINAL_DISTANCE),   # distance (Å, kcal/mol/Å²)
        (D2,  D1,   G1,    ANGLE_TARGET, ANGLE_FC, ANGLE_TARGET),   # angle (deg, kcal/mol/rad²)
        (D1,  G1,   G2,    TORSION_TARGET, TORSION_FC, TORSION_TARGET),   # torsion
    ]:
        r = restraints.DAT_restraint()
        r.mask1 = mask1
        r.mask2 = mask2
        if mask3:
            r.mask3 = mask3
        r.topology = structure
        r.auto_apr = True
        r.continuous_apr = True
        r.amber_index = True
        r.attach["target"] = attach_tgt
        r.attach["fraction_list"] = attach_fractions
        r.attach["fc_final"] = attach_fc
        r.pull["target_final"] = pull_tgt
        r.pull["num_windows"] = windows[1]
        r.initialize()
        guest_restraints.append(r)

    wall_restraints = create_wall_restraints(
        structure, D1, G1, G2, windows[0], host_prefix, guest_prefix
    ) + create_cone_wall(structure, D1, D2, G1, windows[0], host_prefix, guest_prefix)

    # Create directories and write disang.rest for each window.
    window_list = create_window_list(guest_restraints)
    for window in window_list:
        if not os.path.isdir(f"windows/{window}"):
            os.makedirs(f"windows/{window}")

    all_restraints = static_restraints + guest_restraints
    for window in window_list:
        with open(f"windows/{window}/disang.rest", "w") as f:
            for r in all_restraints:
                line = amber_restraint_line(r, window)
                if line is not None:
                    f.write(line)
            if window[0] == "a":
                for wall_restraint in wall_restraints:
                    line = amber_restraint_line(wall_restraint, window)
                    if line is not None:
                        f.write(line)

    save_restraints(guest_restraints, filepath="windows/restraints.json")

    return guest_restraints

def create_windows(G1, G2, host_prefix=HOST, guest_prefix=GUEST):

    guest_restraints = set_restraints(G1,G2, host_prefix, guest_prefix)

    window_list = create_window_list(guest_restraints)

    for window in window_list:
        if window[0] == "a":
            shutil.copy(f"complex/{COMPLEX}-dum.prmtop", f"windows/{window}/{COMPLEX}-dum.prmtop")
            shutil.copy(f"complex/{COMPLEX}-dum.rst7", f"windows/{window}/{COMPLEX}-dum.rst7")
        elif window[0] == "p":
            s = pmd.load_file(f"complex/{COMPLEX}-dum.prmtop", f"complex/{COMPLEX}-dum.rst7", structure=True)
            target_diff = (
                guest_restraints[0].phase["pull"]["targets"][int(window[1:])]
                - guest_restraints[0].pull["target_initial"]
            )
            print(f"Window {window}: translating guest {target_diff.magnitude:0.1f} Å.")
            for atom in s.atoms:
                if atom.residue.name == guest_prefix.upper():
                    atom.xz += target_diff.magnitude
            s.save(f"windows/{window}/{COMPLEX}-dum.prmtop", overwrite=True)
            s.save(f"windows/{window}/{COMPLEX}-dum.rst7", overwrite=True)

    return

def setup_datapoint(dp, data_dir_host, b1_idx, b2_idx, host_prefix=HOST, guest_prefix=GUEST):
    old = os.getcwd()
    os.chdir(dp)
    try:
        build_system(data_dir_host, dp, host_prefix, guest_prefix)
        G1, G2 = guest_masks(b1_idx, b2_idx, guest_prefix)   # needs complex/vac.prmtop
        print("masks:", G1, G2)
        add_dummy(G1, G2, data_dir_host, dp, host_prefix, guest_prefix)
        create_windows(G1, G2, host_prefix, guest_prefix)
    finally:
        os.chdir(old)


def run(dp):
    """APR setup stage: guest/complex PDBs + parameters -> complex/ and windows/.

    The G1/G2 anchors come from guest_input.json, used as written (guest_masks checks they match
    the placed guest). Folders prepared before that file existed get them from the placed guest
    PDB's backbone chain ends, lower end first (as placement picks them), and they are saved.
    An earlier setup's complex/ and windows/ are removed first, so no window from an older
    window list is left behind (there is no simulation output in them: refuse_simulated).
    """
    from rdkit import Chem
    from datapoints import refuse_simulated
    from guest_input import anchors, ensure_guest_input, save_anchors
    from placement import order_anchors
    from structure import backbone_anchors

    refuse_simulated(dp)
    dp = Path(dp).resolve()          # setup_datapoint() changes into dp, so no relative paths
    ensure_guest_input(dp)
    for name in (f"{GUEST}.pdb", f"{GUEST}.mol2", f"{GUEST}.frcmod", f"{COMPLEX}.pdb"):
        if not (dp / name).is_file():
            raise FileNotFoundError(f"APR setup needs {name}: {dp / name}")
    b1_b2 = anchors(dp)
    if b1_b2 is None:
        guest_pdb = dp / f"{GUEST}.pdb"
        mol = Chem.MolFromPDBFile(str(guest_pdb), removeHs=False)
        if mol is None:
            raise ValueError(f"could not read guest PDB: {guest_pdb}")
        b1_b2 = order_anchors(mol.GetConformer().GetPositions(),
                              *backbone_anchors(mol, source=str(guest_pdb)))
        save_anchors(dp, *b1_b2)
    for old in ("complex", "windows"):
        if (dp / old).is_dir():
            shutil.rmtree(dp / old)
    setup_datapoint(str(dp), CFG["directories"]["host"], *b1_b2)


def main(argv=None):
    """python apr_setup.py [NAME ...]: build complex/ and windows/ for each data point."""
    from setup_runner import stage_main
    return stage_main("apr_setup", argv)


if __name__ == "__main__":
    raise SystemExit(main())
