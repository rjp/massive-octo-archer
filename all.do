DEPS="log/log.o plasma.o"
redo-ifchange $DEPS
gcc -o $3 $DEPS -lm
