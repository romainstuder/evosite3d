# setup_view.tcl - Setup view for ubiquitin trajectory.
# Usage:  vmd -e setup_view.tcl
#         (or)  source setup_view.tcl   inside an already-running VMD session

# ----- Configurable paths --------------------------------------------------
set TOP_FILE  "output/step3_minimised.pdb"
set TRAJ_FILE "output/step6_trajectory.dcd"
set RMSF_FILE "output/step7_rmsf.txt"     ;# leave empty "" to disable spectrum coloring

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

# Side chains (hidden — NewCartoon on sidechain selection draws nothing)
mol representation NewCartoon 0.2 12 12
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

# ----- Optional: colour cartoon by RMSF spectrum --------------------------
if {$RMSF_FILE ne "" && [file exists $RMSF_FILE]} {
    set fp [open $RMSF_FILE r]
    set rmsf_list [split [string trim [read $fp]] "\n"]
    close $fp

    set num_res [llength $rmsf_list]
    for {set i 0} {$i < $num_res} {incr i} {
        set v [lindex $rmsf_list $i]
        set sel [atomselect top "protein and residue $i"]
        $sel set beta $v
        $sel delete
    }

    set rmsf_min [lindex $rmsf_list 0]
    set rmsf_max [lindex $rmsf_list 0]
    foreach v $rmsf_list {
        if {$v < $rmsf_min} {set rmsf_min $v}
        if {$v > $rmsf_max} {set rmsf_max $v}
    }

    mol modcolor 0 top Beta
    mol scaleminmax top 0 $rmsf_min [expr {$rmsf_max / 2.0}]
    color scale method BGR

    puts "✓ Cartoon coloured by RMSF (BGYR: blue=rigid → red=flexible)"
    puts "  Range: [format %.3f $rmsf_min] - [format %.3f $rmsf_max] nm  ($num_res residues)"
} else {
    puts "ℹ Spectrum coloring skipped (RMSF_FILE not found or disabled)"
}

# ----- Final view ----------------------------------------------------------
animate goto 0
display resetview
