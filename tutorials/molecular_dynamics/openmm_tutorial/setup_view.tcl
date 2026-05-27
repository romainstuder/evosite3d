# setup_view.tcl - Setup view for ubiquitin trajectory.
# Usage:  vmd -e setup_view.tcl
#         (or)  source setup_view.tcl   inside an already-running VMD session

# ----- Configurable paths --------------------------------------------------
set TOP_FILE  "output/step3_minimised.pdb"
set TRAJ_FILE "output/step6_trajectory.dcd"

# ----- Load files ----------------------------------------------------------
if {![file exists $TOP_FILE]} {
    puts stderr "ERROR: topology not found: $TOP_FILE (cwd: [pwd])"
    exit 1
}
if {![file exists $TRAJ_FILE]} {
    puts stderr "ERROR: trajectory not found: $TRAJ_FILE (cwd: [pwd])"
    exit 1
}

mol new     $TOP_FILE  type pdb waitfor all
mol addfile $TRAJ_FILE type dcd waitfor all

# ----- Representations -----------------------------------------------------
mol delrep 0 top

# Backbone cartoon coloured by secondary structure
mol representation NewCartoon 0.3 10 4.1 0
mol color Structure
mol selection "protein"
mol addrep top

# Side chains (licorice)
mol representation Licorice 0.2 12 12
mol color Name
mol selection "protein and sidechain"
mol addrep top

# ----- Display -------------------------------------------------------------
color Display Background white
display projection Orthographic
display depthcue off
axes location Off

# ----- Align all frames to frame 0 on backbone ----------------------------
set ref [atomselect top "protein and backbone" frame 0]
set sel [atomselect top "protein and backbone"]
set all [atomselect top "all"]

set nframes [molinfo top get numframes]
for {set i 0} {$i < $nframes} {incr i} {
    $sel frame $i
    $all frame $i
    $all move [measure fit $sel $ref]
}

# Free atomselect handles
$ref delete
$sel delete
$all delete

# ----- Final view ----------------------------------------------------------
animate goto 0
display resetview
