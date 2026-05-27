# make_movie.tcl - Load trajectory (via setup_view.tcl), then render PNG frames.
# Usage:
#   vmd -e make_movie.tcl
#   (or)  source make_movie.tcl   inside an already-running VMD session
#
# After rendering:
#   ffmpeg -framerate 25 -i output/movie_frames/frame_%04d.tga \
#          -c:v libx264 -pix_fmt yuv420p output/movie.mp4

# ----- Configurable params -------------------------------------------------
set FRAMES_DIR "output/movie_frames"
set STRIDE     1     ;# render every Nth frame (1 = all)

# ----- Load + set up view if nothing is loaded yet -------------------------
if {[molinfo num] == 0} {
    if {![file exists "setup_view.tcl"]} {
        puts stderr "ERROR: setup_view.tcl not found in [pwd]"
        exit 1
    }
    source setup_view.tcl
}

if {[molinfo num] == 0} {
    puts stderr "ERROR: no molecule loaded after setup_view.tcl. Check paths inside that script."
    exit 1
}

# ----- Render --------------------------------------------------------------
file mkdir $FRAMES_DIR
set num_frames [molinfo top get numframes]
set total_out  [expr {($num_frames + $STRIDE - 1) / $STRIDE}]
puts "Rendering $total_out frames (stride $STRIDE) to $FRAMES_DIR"

set out_idx 0
for {set i 0} {$i < $num_frames} {incr i $STRIDE} {
    animate goto $i
    display update
    display update ui

    set filename [format "$FRAMES_DIR/frame_%04d.tga" $out_idx]
    render TachyonInternal $filename

    puts "  [expr {$out_idx + 1}]/$total_out: frame $i -> $filename"
    incr out_idx
}

puts ""
puts "Done. $out_idx frames in $FRAMES_DIR/"
puts "Combine into a movie:"
puts "  ffmpeg -framerate 25 -i $FRAMES_DIR/frame_%04d.tga \\"
puts "         -c:v libx264 -pix_fmt yuv420p output/movie.mp4"
