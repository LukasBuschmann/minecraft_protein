import pymol
from pymol import cmd

from pdb_2_obj import new_scene, chainbows, save, add_separated_colored_selection


def hsp27():
    # Hsp27 F 29 and 33 additioally tyr 66 trp 53 glu 48
    # TK 48,66,90,93,139

    def set_color(id, col):
        cmd.set("cartoon_color", col, id)
        cmd.set("surface_color", col, id)
        cmd.set("stick_color", col, id)

    cmd.delete("all")
    cmd.load("proteins/1osn.cif", "tk")
    cmd.remove("not chain A")
    cmd.remove("resn ADP")
    cmd.select("bvdu", "tk and resn BVP")
    cmd.select("bs", "tk and resi 93 or resi 139")

    set_color("tk", "palegreen")
    set_color("bs", "red")
    set_color("bvdu", "orange")

    cmd.show("sticks", "bvdu")
    cmd.show("sticks", "bs")
    cmd.show("cartoon", "tk")

    cmd.load("proteins/SAM-T08-model1-view_ConSurf.pdb", "hsp27")
    cmd.select("bs2", "hsp27 and resi 29 or (hsp27 and resi 33)")
    set_color("bs2", "red")
    cmd.show("sticks", "bs2")
    set_color("hsp27", "skyblue")
    cmd.show("cartoon", "hsp27")
    cmd.deselect()
    cmd.save("hsp27.wrl")



def gyra():
    cmd.delete("all")

    try:
        cmd.load("2y3p.pdb", "gyra")
    except:
        cmd.fetch("2y3p", "gyra", type="pdb")

    cmd.remove("chain B")

    cmd.hide("everything")

    cmd.set("cartoon_color", "palegreen", "gyra")
    cmd.set("surface_color", "palegreen", "gyra")

    cmd.select("resis", "gyra and resi 83+87")
    cmd.set("cartoon_color", "red", "resis")
    cmd.set("surface_color", "red", "resis")

    cmd.select("sm8", "resn SM8")
    cmd.set("stick_color", "orange", "sm8")

    cmd.select("pocket", "gyra and chain A within 5.0 of sm8")
    cmd.set("cartoon_color", "salmon", "pocket")
    cmd.set("surface_color", "salmon", "pocket")

    cmd.orient("gyra")
    cmd.zoom("gyra", -30)
    # cmd.show("cartoon", "gyra")
    cmd.show("surface", "gyra")
    cmd.show("sticks", "sm8")
    cmd.save("gyra.wrl")


def ribosome():
    scene = new_scene()
    cmd.load("proteins/4v6x.cif")
    scene = add_separated_colored_selection(scene, "chain A5+B2", "cartoon", (0.9, 0.8, 0.9))
    cmd.remove("chain A5+B2")
    scene = chainbows(scene, "4v6x", "surface", separate=True)
    save(scene, "models/ribosome", "ribosome")


if __name__ == '__main__':
    pymol.finish_launching()
    ribosome()