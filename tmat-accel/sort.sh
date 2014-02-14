for file in *.pickl; do
    new_path=$(perl -pe 's/(result-h-runcopy-)([0-9]+)(-[0-9]+\.pickl)/h-run-$2\/$1$2$3/' <<< "$file")
    mkdir "${new_path%/*}"
    mv -vn "$file" "$new_path"
done
