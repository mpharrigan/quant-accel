for file in *.pickl; do
    new_path=$(perl -pe 's/(result-j-runcopy-)([0-9]+)(-[0-9]+\.pickl)/j-run-$2\/$1$2$3/' <<< "$file")
    mkdir "${new_path%/*}"
    mv -vn "$file" "$new_path"
done
