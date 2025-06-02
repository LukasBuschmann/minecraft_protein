import itertools
import os
import colorsys
import shutil
from collections import defaultdict

import trimesh

from pathlib import Path

from pymol import cmd
from tqdm import tqdm

# this will not work on its own, since we don't specify a target .mtl
def add_material(obj_path, material_name):
    obj = open(obj_path).read()
    with open(obj_path, 'w') as f:
        color = f"usemtl {material_name}\n"
        obj = color + obj # use material for complete object
        f.write(obj)

def get_chain_colors(chains, shift=0.05):
    colors = {}
    for i, chain in enumerate(chains):
        hue = i / len(chains)
        colors[chain] = {
            # todo: missing edge case for hue +/- shift exceeding (-)1
            "H": colorsys.hsv_to_rgb(hue + shift, 0.7, 0.8),
            "S": colorsys.hsv_to_rgb(hue - shift, 0.7, 0.8),
            "L": colorsys.hsv_to_rgb(hue, 0.8, 0.5),
        }
    return colors

def mtl_from_colors(colors, mtl_path):
    mtl = ""
    material_names = defaultdict(dict)
    for chain, color_variants in colors.items():
        for ss, color in color_variants.items():
            material_name = f"{chain}_{ss}"
            material_names[chain][ss] = material_name
            mtl += f"newmtl {material_name}\n"
            mtl += f"Kd {' '.join(map(str, color))}\n"

    with open(mtl_path, 'w') as f:
        f.write(mtl)

    return material_names


def combine_obj(obj_paths, mtl_filename, output_path):
    all_vertices = []
    all_face_blocks = []
    vertex_offset = 0
    mtllib_written = False

    for path in obj_paths:
        with open(path, 'r') as f:
            lines = f.readlines()

        name = path.split('/')[-1].split('.')[0]  # e.g., "model1" from "model1.obj"
        vertices = [line for line in lines if line.startswith('v ')]
        faces = [line for line in lines if line.startswith('f ')]
        usemtl_lines = [line for line in lines if line.startswith('usemtl ')]

        # Use first material line in the file or fallback to default
        usemtl = usemtl_lines[0].strip() if usemtl_lines else f"usemtl {name}"

        # Clean face definitions, offset vertex indices
        adjusted_faces = []
        for face in faces:
            parts = face.strip().split()
            new_face = ['f']
            for part in parts[1:]:
                indices = part.split('/')
                v_idx = int(indices[0]) + vertex_offset
                new_face.append(str(v_idx))  # Only use vertex index
            adjusted_faces.append(' '.join(new_face) + '\n')

        all_vertices.extend(vertices)

        # Group this object's data: name + material + faces
        face_block = [f"o {name}\n", f"{usemtl}\n"]
        face_block.extend(adjusted_faces)
        all_face_blocks.append(face_block)

        vertex_offset += len(vertices)

    # Write combined file
    with open(output_path, 'w') as out:
        if not mtllib_written:
            out.write(f"mtllib {mtl_filename}\n")
            mtllib_written = True
        out.writelines(all_vertices)
        for block in all_face_blocks:
            out.writelines(block)


def show(path):
    mesh = trimesh.load(path)
    mesh.show()

def exists_or_create_dir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

if __name__ == '__main__':

    model_path = "models"
    existing_models = list(map(lambda p: p.stem, Path(model_path).iterdir()))

    protein_path = "proteins"

    proteins = Path(protein_path).iterdir()

    overwrite = False
    if not overwrite:
        proteins = list(filter(lambda p: p.stem not in existing_models, proteins))

    for protein_path in proteins:
        print(protein_path.stem)
        # if protein_path.stem != "1nu2":
        #     continue

        cmd.delete("all")

        protein = "protein_object"
        cmd.load(protein_path, protein)
        chains = cmd.get_chains(protein)
        protein_name = protein_path.stem

        exists_or_create_dir(model_path)
        tmp_path = "tmp"
        exists_or_create_dir(tmp_path)

        colors = get_chain_colors(chains)
        mtl_path = f"{model_path}/{protein_name}.mtl"
        material_names = mtl_from_colors(colors, mtl_path) # also saves .mtl file

        chain_objs = []

        for chain in tqdm(chains):
            for ss in ["H", "S", "L"]:
                material = material_names[chain][ss]
                chain_path = f"{tmp_path}/{protein_name}_{chain}_{ss}.obj"
                selection_name = f"chain_{chain}_{ss}" # some letters are reserved keywords in pymol

                selection_str = f"{protein} and chain {chain} and ss {ss}" if ss != 'L' else f"{protein} and chain {chain} and not ss H+S"
                cmd.select(selection_name, selection_str)
                # widen selection to avoid discontets in chain (due to pymol cartoon rendering)
                cmd.select(selection_name, f"{selection_name} extend 2")
                cmd.hide("all")
                cmd.show("cartoon", selection_name)
                cmd.save(chain_path)
                add_material(chain_path, material)
                chain_objs.append(chain_path)


        full_obj_path = f"{model_path}/{protein_name}.obj"
        combine_obj(chain_objs, f"{protein_name}.mtl", full_obj_path)
        shutil.rmtree(tmp_path)


