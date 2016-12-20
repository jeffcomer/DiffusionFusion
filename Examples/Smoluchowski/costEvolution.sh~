#!/bin/bash

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
