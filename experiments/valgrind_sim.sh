valgrind --tool=callgrind gmx mdrun -v -deffnm em
callgrind_annotate --auto=yes callgrind.out.pid
