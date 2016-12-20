#!/bin/bash
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

if (( $# != 1 && $# != 2 && $# != 3)); then
    echo "Usage: $0 logField [fieldName] [secondName]"
    exit
fi

log=$1
second=$3
bas=${log%.*}

if (( $# == 1 )); then
    field=costUpdate
else
    field=$2
fi

if (( $# == 3 )); then
    paste -d " " <(awk 'NF >= 2 && /^cycle / {print $2}' $log)  <(awk "/^  $field $second / {print \$3}" $log) > $bas.$field.$second.dat
else
    paste -d " " <(awk 'NF >= 2 && /^cycle / {print $2}' $log)  <(awk "/^  $field/ {print \$2}" $log) > $bas.$field.dat
fi
