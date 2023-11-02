# sub-dddd/
#         - sub-dddd_sdit-*_pet.json

dest=$(echo $1 | rg -o --pcre2 'sub-\d{6,8}.*(?=\/data)')
echo "moving to $dest"
mv $1 $dest/