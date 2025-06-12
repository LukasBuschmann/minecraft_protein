import colorsys
import os

import trimesh
from pymol import cmd


def face_tuple(face, offset):
    vertices = face.split(' ')[1:] # f v1 v2 v3
    return [str(int(vertex.split('/')[0]) + offset) for vertex in vertices] # v1 := vId/tId/nId

def add_obj(scene, obj, color) -> tuple[int, int, str, str]:
    num_obj, offset, objs, mtls = scene

    elements = obj.split('\n')
    v = [element for element in elements if element.startswith('v ')]
    f = [f"f {' '.join(face_tuple(element, offset))}" for element in elements if element.startswith('f ')]
    obj = f"o obj{num_obj+1}\n"
    obj += "\n".join(v)
    obj += f"\nusemtl mtl{num_obj+1}\n"
    obj += "\n".join(f)
    obj += "\n"
    mtls += f"newmtl mtl{num_obj+1}\n"
    mtls += f"Kd {' '.join(map(str, color))}\n"
    mtls += "\n"

    return num_obj + 1, offset + len(v), objs + obj, mtls

def show(path):
    mesh = trimesh.load(path)
    mesh.show()

def add_colored_selection(scene, selection, rep, color) -> tuple:
    cmd.hide("all")
    cmd.show(rep, selection)
    cmd.save("tmp.obj")
    obj = open("tmp.obj").read()
    os.remove("tmp.obj")
    return add_obj(scene, obj, color)

def add_separated_colored_selection(scene, selection, rep, color) -> tuple:
    cmd.hide("all")
    cmd.create("tmp_selection", selection)
    cmd.show(rep, "tmp_selection")
    cmd.save("tmp.obj")
    obj = open("tmp.obj").read()
    os.remove("tmp.obj")
    return add_obj(scene, obj, color)

def new_scene() -> tuple:
    return 0, 0, "", ""

def rainbow(n, saturation, value): # 0.7, 0.8 looks pretty
    colors = []
    for i in range(n):
        hue = i / n
        colors.append(colorsys.hsv_to_rgb(hue, saturation, value))
    return colors

def chainbows(scene, name, rep, separate=False) -> tuple:
    chains = cmd.get_chains(name)
    colors = rainbow(len(chains), 0.7, 0.8)

    for chain, color in zip(chains, colors):
        if separate:
            scene = add_separated_colored_selection(scene, f"chain {chain}", rep, color)
        else:
            scene = add_colored_selection(scene, f"chain {chain}", rep, color)

    return scene

def save(scene, path, name):
    _, _, objs, mtls = scene
    open(path + f"/{name}.obj", 'w').write(f"mtllib {name}.mtl\n" + objs)
    open(path + f"/{name}.mtl", 'w').write(mtls)

if __name__ == '__main__':
    pass


