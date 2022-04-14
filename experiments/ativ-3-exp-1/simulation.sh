echo "Creating environment"
echo "15" | gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce > dummy
gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic > dummy
gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top > dummy
gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr > dummy
echo "13" | gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral > dummy
gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr > dummy
