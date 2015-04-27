redo-ifchange $2.c
gcc -O3 -fomit-frame-pointer -funroll-loops -Ilog -g -MD -MF $2.d -c -o $3 $2.c
read DEPS <$2.d
redo-ifchange ${DEPS#*:}
