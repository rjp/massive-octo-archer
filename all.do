DEPS="log/log.o plasma.o"
redo-ifchange $DEPS
gcc -Ofast -Ilog -o $3 $DEPS -lm
