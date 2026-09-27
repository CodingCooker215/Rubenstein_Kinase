curr_pdb=[]
data=[]
#filely=open("AukB_t232a_94_c0_5121024_alignment_DFG_active_test.txt", "w")

#DFG=142-144
#G_Loop=8-13
#aC=46-53

def annoatate():
    cmd.remove("Chain B")
    cmd.remove("Chain C")
    cmd.remove("Chain D")
    cmd.remove("Chain E")
    cmd.set("cartoon_transparency", 0.8, )

    cmd.select("DFG", "resi 142-144")
    cmd.color_deep("Green", "DFG")
    cmd.set("cartoon_transparency", 0, "DFG")

    cmd.select("aC_helix", "resi 43-53")
    cmd.color_deep("Wheat", "aC_helix")
    cmd.set("cartoon_transparency", 0, "aC_helix")


    cmd.select("F219", "resi 143")
    cmd.show("sticks", "F219")
    cmd.color_deep("pink", "F219")
    cmd.set("cartoon_transparency", 0, "F219")

    cmd.select("D218", "resi 142")
    cmd.show("sticks", "D218")
    cmd.color_deep("red", "D218")
    cmd.set("cartoon_transparency", 0, "D218")

    cmd.intra_fit("aC_helix")
def MD_run():
    # Improve rendering quality
    # set ray_trace_mode, 1
    # set antialias, 2
    # set ray_shadows, on
    # set ambient, 0.2
    # set specular, 0.3


    n = cmd.count_states()
    obj_name = cmd.get_object_list()[0]
    output_dir = f"/Volumes/Jason/Simulations/{obj_name}_frames"
    os.makedirs(output_dir, exist_ok=True)


    for i in range(1, n+1, 100):
        cmd.frame(i)
        cmd.center("aC_helix")
        cmd.zoom("aC_helix", 5)

        filename = os.path.join(output_dir, f"{obj_name}_frame_{i:04d}.png")
        
        cmd.png(
            filename,
            width=1920,
            height=1080,
            dpi=300,
            ray=1
        )

    print(f"Saved {filename}")
    
annoatate()
MD_run()