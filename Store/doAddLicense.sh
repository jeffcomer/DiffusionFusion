for f in *.H *.C; do
    cat license_header.txt $f > ../$f
done
