#!/usr/bin/zsh
mrn=$(echo $1 | tr -d ".tar.gz")
tar -xzf $1 data/h_vmac/koran/FDG_PET/$mrn --strip-components=4
find $mrn/*/*/*  ! \( -name "*.json" -o -name "*.png" -o -name "*.dcm" \) -type "f" -exec mv {} {}.dcm \;
zip -rTm $mrn.zip $mrn
