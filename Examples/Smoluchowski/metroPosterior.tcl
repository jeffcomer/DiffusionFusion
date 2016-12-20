#!/usr/bin/tclsh
# Author: Jeff Comer <jeffcomer at gmail>

#////////////////////////////////////////////////////////////////////
# Copyright 2014-2016 Jeffrey Comer
#
# This file is part of DiffusionFusion.
#
# DiffusionFusion is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# DiffusionFusion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with DiffusionFusion. If not, see http://www.gnu.org/licenses/.
#/////////////////////////////////////////////////////////////////////

if {[llength $argv] != 5} {
    puts "Usage: tclsh $argv0 trajFile templateFile firstCycle stride outName"
    exit
}

foreach {trajFile templateFile firstCycle stride outName} $argv { break }

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	set line [string trim $line]
	if {[string match "#*" $line]} { continue }
	if {[string length $line] < 1} { continue }
	
	set tok [concat $line]
	set l {}
	foreach t $tok {
	    lappend l [string trim $t]
	}

	lappend r $l
    }

    close $in
    return $r
}

# Just read a space delimited data file.
proc readTraj {fileName firstCycle stride arrayVar} {
    upvar $arrayVar arr

    set in [open $fileName r]

    set readField {}
    set fieldList {}
    set cycleCount 0
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} { continue }
	if {[string length [string trim $line]] < 1} { continue }
	set tok [concat $line]

	# Store the previous field, if we were reading one.
	if {[string match "CYCLE*" $line] || [string match "FIELD*" $line]} {
	    # Ignore cycles before first cycle.
	    if {[llength $readField] > 0} {
		# Store.
		lappend arr($readField) $numList
		lappend arr(${readField}_cost) $cost
		puts "CYCLE $cycle $cost"
		puts "field $readField length [llength $numList]"
		set readField {}
		set numList {}
	    }
	}

	# Read the cycle info.
	if {[string match "CYCLE*" $line]} {
	    set cycle [lindex $tok 1]
	    set cost [lindex $tok 2]
	    incr cycleCount
	}

	# Skip cycles that are less than firstCycle or don't fit the stride.
	if {$cycle >= $firstCycle && $cycleCount % $stride == 0} {

	    if {[string match "FIELD*" $line]} {
		set readField [lindex $tok 1]
		# Keep a list of the fields that have been read.
		if {[lsearch $fieldList $readField] < 0} {
		    puts "Found field $readField"
		    lappend fieldList $readField
		}
	    } else {
		# Store the numeric data.
		if {[llength $readField] > 0 && [string is double [lindex $tok 0]]} {
		    lappend numList [lindex $tok 0]
		}
	    }
	}
	
    }

    # Get the last field that we were reading.
    # Store the previous field, if we were reading one.
    if {[llength $readField] > 0} {
	# Store.
	lappend arr($readField) $numList
	lappend arr(${readField}_cost) $cost
	puts "CYCLE $cycle $cost"
	puts "field $readField length [llength $numList]"
	set readField {}
	set numList {}
    }


    return $fieldList
}

set templateData [readData $templateFile]
set fieldList [readTraj $trajFile $firstCycle $stride arr]

set fieldNum [llength $fieldList]
puts "$fieldNum fields in `$trajFile'"

foreach field $fieldList {
    set samples [llength $arr($field)] 
    puts "$samples samples in field `$field'"
    set nodeNum [llength [lindex $arr($field) 0]]
    puts "$nodeNum nodes in field `$field'"
    if {$nodeNum != [llength $templateData]} {
	puts "Warning. `$templateFile' has only [llength $templateData] nodes (should have $nodeNum). Skipping field `$field'."
    }

    # Find the minimum cost.
    set minCost [lindex $arr(${field}_cost) 0]
    foreach cost $arr(${field}_cost) {	
	if {$cost < $minCost} { set minCost $cost }
    }

    # Zero the field sums.
    for {set i 0} {$i < $nodeNum} {incr i} {
	set sum($i) 0.0
	set sumSq($i) 0.0
	set count($i) 0
    }

    # Calculate the field sums.
    foreach frame $arr($field) {
	for {set i 0} {$i < $nodeNum} {incr i} {
	    set val [lindex $frame $i]
	    set sum($i) [expr {$sum($i) + $val}]
	    set sumSq($i) [expr {$sumSq($i) + $val*$val}]
	    incr count($i)
	}
    }

    # Write the output.
    set out [open $outName.$field.dat w]
    for {set i 0} {$i < $nodeNum} {incr i} {
	set mean [expr {$sum($i)/$count($i)}]
	set meanSq [expr {$sumSq($i)/$count($i)}]
	set s [expr {$meanSq - $mean*$mean}]
	set z [lindex $templateData $i 0]

	if {$s < 0} {
	    puts "Warning! <${field}^2> - <${field}>^2 < 0, presumably due to numerical error."
	    puts "We take the standard deviation to be zero."
	    puts "index $i <${field}^2> $sumSq($i) <${field}>^2 [expr {$mean*$mean}] variance $s"
	    set std 0.0
	} else {
	    set std [expr {sqrt($s)}]
	}
	puts $out "$z $mean $std"
    }

    close $out
}
