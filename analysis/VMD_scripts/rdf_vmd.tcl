# Define a procedure to compute g(r) and write to file
proc compute_gofr {sel1 sel2 output delta rmax usepbc selupdate first last step} {
    # Calculate g(r)
    set gr [measure gofr $sel1 $sel2 delta $delta rmax $rmax usepbc $usepbc selupdate $selupdate first $first last $last step $step]

    # Write g(r) to file
    set outfile [open $output w]
    set r [lindex $gr 0]
    set gr2 [lindex $gr 1]
    set igr [lindex $gr 2]
    foreach j $r k $gr2 l $igr {
       puts $outfile "$j $k $l"
    }
    close $outfile
}

set first_step 15000
set last_step -1
set delta 0.05
set rmax 7
set step 1

# Ge-Ge
compute_gofr [atomselect top "name Ge"] [atomselect top "name Ge"] "GeGe_rdf.dat" $delta $rmax 1 1 $first_step $last_step $step

# Ge-Te
compute_gofr [atomselect top "name Ge"] [atomselect top "name Te"] "GeTe_rdf.dat" $delta $rmax 1 1 $first_step $last_step $step

# Te-Te
compute_gofr [atomselect top "name Te"] [atomselect top "name Te"] "TeTe_rdf.dat" $delta $rmax 1 1 $first_step $last_step $step

# total rdf
compute_gofr [atomselect top "name Ge Te"] [atomselect top "name Ge Te"] "total_rdf.dat" $delta $rmax 1 1 $first_step $last_step $step

