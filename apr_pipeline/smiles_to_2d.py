#!/usr/bin/env python
"""Draw 2D SMILES structures with the pipeline's restraint anchors highlighted.

Input lines can be any of
    SMILES
    SMILES <whitespace> name
    name <whitespace> SMILES
    anything<TAB>SMILES        (e.g. "iso_0022 FAILED: clash ...<TAB>CC(F)...")
Blank lines and lines beginning with '#' are skipped. Windows line endings and a
UTF-8 BOM are handled. Lines that cannot be read as SMILES are reported/skipped.

Usage:
    python smiles_to_2d.py input.txt [outdir] [--per-page 24] [--cols 4]
                                      [--size 300] [--single]

Output in outdir (default ./structures_2d):
    grid_001.png, grid_002.png, ...   paginated grids
    single/<name>.png                 one image per molecule (only with --single)
"""
import argparse
from io import BytesIO
import math
from pathlib import Path
import re
import sys

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D, Point3D

PIPELINE_MODULE_DIR = Path(__file__).resolve().parents[1] / "apr_pipeline"
sys.path.insert(0, str(PIPELINE_MODULE_DIR))

from structure import find_carbon_backbone_chain, pick_backbone_alignment_atoms

RDLogger.DisableLog("rdApp.*")

GUIDE_EXTENTS = ((-4.0, "solid"), (-3.5, "dotted"), (3.5, "dotted"), (4.0, "solid"))
CELL_FOOTER_HEIGHT = 28


def parse(line, number):
    """Return (name, molecule) or (None, None)."""
    if "\t" in line:
        label, smiles = line.rsplit("\t", 1)
        candidates = [
            (smiles.strip(), label.split()[0] if label.split() else f"mol_{number}")
        ]
    else:
        tokens = line.split()
        candidates = []
        for index, token in enumerate(tokens):
            other_tokens = [value for i, value in enumerate(tokens) if i != index]
            candidates.append(
                (token, other_tokens[0] if other_tokens else f"mol_{number}")
            )

    for smiles, name in candidates:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            return name, molecule
    return None, None


def prepare_molecule(molecule):
    """Compute 2D coordinates and put the pipeline's restraint atoms on x=0."""
    chain = find_carbon_backbone_chain(molecule)
    if len(chain) < 2:
        raise ValueError("the pipeline found fewer than two backbone atoms")

    anchors = pick_backbone_alignment_atoms(chain)
    if len(anchors) != 2:
        raise ValueError("the pipeline restraint selector must return two atom indices")

    rdDepictor.Compute2DCoords(molecule)
    conformer = molecule.GetConformer()
    first = conformer.GetAtomPosition(anchors[0])
    second = conformer.GetAtomPosition(anchors[1])
    dx, dy = second.x - first.x, second.y - first.y
    if dx == 0 and dy == 0:
        raise ValueError("the pipeline restraint atoms have identical 2D coordinates")

    angle = math.pi / 2 - math.atan2(dy, dx)
    cosine, sine = math.cos(angle), math.sin(angle)
    midpoint_x = (first.x + second.x) / 2
    midpoint_y = (first.y + second.y) / 2
    for atom_index in range(molecule.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        x, y = point.x - midpoint_x, point.y - midpoint_y
        conformer.SetAtomPosition(
            atom_index,
            Point3D(x * cosine - y * sine, x * sine + y * cosine, 0.0),
        )
    return molecule, anchors


def _draw_guides(image_draw, drawer, y_min, y_max):
    top = drawer.GetDrawCoords(Point2D(0.0, y_max))
    bottom = drawer.GetDrawCoords(Point2D(0.0, y_min))
    y_top, y_bottom = round(top.y), round(bottom.y)
    if y_top > y_bottom:
        y_top, y_bottom = y_bottom, y_top

    for x_value, style in GUIDE_EXTENTS:
        point = drawer.GetDrawCoords(Point2D(x_value, 0.0))
        x = round(point.x)
        color = (190, 35, 35) if style == "solid" else (225, 105, 105)
        if style == "solid":
            image_draw.line((x, y_top, x, y_bottom), fill=color, width=2)
            continue

        dash_length, gap_length = 5, 5
        y = y_top
        while y < y_bottom:
            image_draw.line(
                (x, y, x, min(y + dash_length - 1, y_bottom)),
                fill=color,
                width=2,
            )
            y += dash_length + gap_length


def draw_molecule(molecule, name, size, anchors):
    """Render a molecule, schematic cavity guides, and red restraint crosses."""
    conformer = molecule.GetConformer()
    points = [conformer.GetAtomPosition(index) for index in range(molecule.GetNumAtoms())]
    x_values = [point.x for point in points]
    y_values = [point.y for point in points]
    x_min = min(min(x_values), -4.0) - 0.8
    x_max = max(max(x_values), 4.0) + 0.8
    y_min = min(y_values) - 1.0
    y_max = max(y_values) + 1.0

    drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
    drawer.SetScale(
        size,
        size,
        Point2D(x_min, y_min),
        Point2D(x_max, y_max),
        molecule,
    )
    drawer.DrawMolecule(molecule, legend=name)
    drawer.FinishDrawing()
    image = Image.open(BytesIO(drawer.GetDrawingText())).convert("RGB")
    image_draw = ImageDraw.Draw(image)
    _draw_guides(image_draw, drawer, y_min, y_max)
    font = ImageFont.load_default()
    for index, anchor_name in zip(anchors, ("b1", "b2")):
        point = drawer.GetDrawCoords(index)
        x, y = int(point.x), int(point.y)
        image_draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill="white")
        image_draw.line((x - 6, y - 6, x + 6, y + 6), fill=(210, 35, 35), width=3)
        image_draw.line((x - 6, y + 6, x + 6, y - 6), fill=(210, 35, 35), width=3)
        image_draw.text((x + 8, y - 16), anchor_name, fill=(175, 20, 20), font=font)

    cell = Image.new("RGB", (size, size + CELL_FOOTER_HEIGHT), "white")
    cell.paste(image, (0, 0))
    footer = ImageDraw.Draw(cell)
    footer_font = ImageFont.load_default()
    footer.text((4, size + 1), "Red X: b1/b2 at x=0 | dotted: +/-7 A", fill="black", font=footer_font)
    footer.text((4, size + 14), "Solid: +/-8 A (schematic guides; not a measured clash)", fill="black", font=footer_font)
    return cell


def _safe_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._") or "molecule"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile")
    parser.add_argument("outdir", nargs="?", default="structures_2d")
    parser.add_argument("--per-page", type=int, default=24, help="molecules per grid image")
    parser.add_argument("--cols", type=int, default=4, help="columns per grid")
    parser.add_argument("--size", type=int, default=300, help="pixel size of each molecule cell")
    parser.add_argument("--single", action="store_true", help="also write one PNG per molecule")
    args = parser.parse_args(argv)
    if args.per_page < 1 or args.cols < 1 or args.size < 1:
        parser.error("--per-page, --cols, and --size must be positive")

    input_path = Path(args.infile)
    if not input_path.is_file():
        sys.exit(f"ERROR: input file not found: {input_path.resolve()}")

    molecules, names, anchors_by_molecule, skipped = [], [], [], []
    with input_path.open(encoding="utf-8-sig", errors="replace") as input_file:
        for number, raw_line in enumerate(input_file, 1):
            line = raw_line.strip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            name, molecule = parse(line, number)
            if molecule is None:
                skipped.append((number, line[:60]))
                continue
            try:
                molecule, anchors = prepare_molecule(molecule)
            except ValueError as error:
                skipped.append((number, f"{line[:45]} ({error})"))
                continue
            molecules.append(molecule)
            names.append(name)
            anchors_by_molecule.append(anchors)

    print(f"Read {input_path}: {len(molecules)} molecules parsed, {len(skipped)} lines skipped")
    for number, line in skipped[:10]:
        print(f"  skipped line {number}: {line!r}")
    if len(skipped) > 10:
        print(f"  ... and {len(skipped) - 10} more")
    if not molecules:
        sys.exit("ERROR: no valid SMILES with two pipeline restraint atoms; no images were made.")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    cells = [
        draw_molecule(molecule, name, args.size, anchors)
        for molecule, name, anchors in zip(molecules, names, anchors_by_molecule)
    ]
    page_count = (len(cells) + args.per_page - 1) // args.per_page
    cell_height = args.size + CELL_FOOTER_HEIGHT
    for page in range(page_count):
        start = page * args.per_page
        page_cells = cells[start : start + args.per_page]
        rows = (len(page_cells) + args.cols - 1) // args.cols
        grid = Image.new("RGB", (args.cols * args.size, rows * cell_height), "white")
        for index, cell in enumerate(page_cells):
            x = (index % args.cols) * args.size
            y = (index // args.cols) * cell_height
            grid.paste(cell, (x, y))
        grid.save(outdir / f"grid_{page + 1:03d}.png")
    print(f"Wrote {page_count} grid image(s) to {outdir}")

    if args.single:
        single_dir = outdir / "single"
        single_dir.mkdir(parents=True, exist_ok=True)
        for index, (cell, name) in enumerate(zip(cells, names), 1):
            cell.save(single_dir / f"{index:05d}_{_safe_filename(name)}.png")
        print(f"Wrote {len(cells)} single images to {single_dir}")


if __name__ == "__main__":
    main()
