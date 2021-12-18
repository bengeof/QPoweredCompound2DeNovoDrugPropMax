

#---------  FILE SETUP  -----------
#PDB FILE, will be replace in the code later by $PDB
FILE="7vh8" #<<<<<<<<<<<<<<<<<<<<<<<<< PUT THE PDB NAME HERE (without the extension)
# LIGAND NAME, if you have a ligand, it will be parametrize with acpype and the ligand name will be replace by "LIG".
LIGNAME="7DI" #<<<<<<<<<<<<<<<<<<<<<<  #PUT LIGAND NAME HERE, leave it blank if no ligand.

#---------  SIMU SETUP  -----------
BOXSIZE=1.2 #cubic simulation boxsiwe
BOXTYPE=cubic #Box type
NT=8 #Number of core.
WATER=tip3p #Water type
NUMBEROFREPLICAS=1 #Number of replica
FF=amber99sb-ildn #Force field
SIMULATIONTIME=100 #Simulation time in nanosec. Will be converted in fs and modified in the mdp file.

#---------  HPC SETUP  -----------
MPI="" #If you have to submit jobs with MPI softwares like "mpirun -np 10". Add the command here
GMX=gmx #GMX command (can be "$GMX_mpi" sometimes. Just change it here

if [ ! -z "$LIGNAME" ]
then
      #Create parameters directories
      mkdir param
      cp $FILE".pdb" param/
      cd param

      grep 'ATOM  ' $FILE".pdb" --color=none > receptor.pdb

      #Extract ligand and connect
      python_command=$(python <<EOF
ligand_atom = []
keepLine = []
with open("$FILE.pdb","r") as file:
    lines = file.readlines()
    for line in lines:
        if '$LIGNAME' in line[17:20]:
            line = line[:17]+"LIG"+line[20:]
            keepLine.append(line)
            ligand_atom.append(int(line[6:11]))
        elif "CONECT" in line[0:6]:
            idx = [int(x) for x in line.split()[1:]]
            if any(id in idx for id in ligand_atom):
                keepLine.append(line)
with open("ligand.pdb","w") as file:
    for line in keepLine:
        file.write(line)
EOF
)


    #Convert in mol2 while adding hydrogens.
    obabel -ipdb ligand.pdb -omol2 -h > ligand.mol2

    #use ACPYPE to create the ligand topology.
    #DISCLAMER! This is a "quick and dirty method", it has to be optimised with ACPYPE parameters of course and adapted to ligands
    #if you see strange MD behaviours.
    # You may also consider Automated Topology Builder (ATB) (webserver) Or LibParGen (webserver & standalone tools)
    acpype -i ligand.mol2
    mkdir ligand
    mv ligand* ligand/
    mkdir receptor
    mv receptor.pdb receptor/
    cd receptor

    #Preparing topology
    $GMX pdb2gmx -f receptor.pdb -o receptor_GMX.pdb -water $WATER -ignh -ff $FF

    #Copy files from receptors/ligand folders.
    cd ../../
    cp param/receptor/*.itp param/receptor/topol.top .
    cp param/ligand/ligand.acpype/ligand_GMX.itp ligand.itp
    grep -h ATOM param/receptor/receptor_GMX.pdb param/ligand/ligand.acpype/ligand_NEW.pdb > complex.pdb

    #Add the ligand topology inside "topol.top"
    cp topol.top topol.bak
    cat topol.top | sed '/forcefield\.itp\"/a\
  #include "ligand.itp"
  ' > topol2.top
    mv topol2.top topol.top
    echo "ligand   1" >> topol.top

#Generate idx group for ligand without hydrogens (for restraints)
ndx=$($GMX make_ndx -f param/ligand/ligand.acpype/ligand_NEW.pdb -o lig_noh.ndx <<EOF
r LIG & !a H*
name 3 LIG-H
q
EOF
)	
    echo "LIG-H" | $GMX genrestr -f param/ligand/ligand.acpype/ligand_NEW.pdb -o posre_ligand.itp -n lig_noh.ndx -fc 1000 1000 1000
	
	#Include posre_ligand.itp AT THE END!!!!!!! of  ligand.itp
	echo '
	
	 ; Include Position restraint file
#ifdef POSRES
#include "posre_ligand.itp"
#endif'  >> ligand.itp
    	
    $GMX editconf -f  complex.pdb -o intial_newbox.gro -bt cubic -d 3.3
    echo "13" | $GMX editconf -f intial_newbox.gro -o complex_newbox.gro -princ 1
	PDB=complex
else

	#set "PDB" name (all the simulation filenames are based on this variable).
	PDB=$FILE
	########################
	##   TOPOLOGIE Creation
	#######################
	$GMX pdb2gmx -f $PDB".pdb" -o $PDB"_processed.gro" -water $WATER -ignh -ff $FF
	########################
	##   Solvatation
	#######################
	#Default editconf, changer if you want...
	$GMX editconf -f  $PDB"_processed.gro" -o $PDB"_newbox.gro" -bt cubic -d 3.2
fi



#SOLVATATION
$GMX solvate -cp $PDB"_newbox.gro" -cs spc216.gro -o $PDB"_solv.gro" -p topol.top

#######################
## ADDING IONS
#######################
$GMX grompp -f ions.mdp -c $PDB"_solv.gro" -p topol.top -o ions.tpr --maxwarn 2

echo "SOL" | $GMX genion -s ions.tpr -o $PDB"_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral


$GMX grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -v -gpu_id 0 -deffnm em

$GMX grompp -f npt.mdp -c em.gro -p topol.top -r em.gro -o npt.tpr
$GMX mdrun -gpu_id 0  -deffnm npt



$GMX grompp -f md_prod.mdp -c npt.gro -p topol.top -o md_out.tpr -maxwarn 2
$GMX mdrun -gpu_id 0 -deffnm md_out 

restr=$($GMX genrestr -f md_out.gro -o Protein_chain_A.itp <<EOF
2
EOF
)

sed '/restraint file/i #ifdef POSRES_B \n#include "Protein_chain_A.itp" \n#endif' topol.top > topol_new.top


cp topol.top topol1.bak
mv topol_new.top topol.top

restrnw=$($GMX make_ndx -f md_out.gro <<EOF
name 2 Protein_chain_A
name 19 ligand
q
EOF
)


$GMX grompp -f md_pull.mdp -c md_out.gro -p topol.top -r md_out.gro -n index.ndx -t md_out.cpt -o pull.tpr -maxwarn 2
$GMX mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg -gpu_id 0 



gmx=gmx
echo 0 | $gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep
#Measuring the COM distance between the two groups
for ps in `seq 1 156`;do 
$gmx distance -s pull.tpr -f conf${ps}.gro -n index.ndx -select 'com of group "Protein_chain_A" plus com of group "ligand"' -oall dist${ps}.xvg 
d=`tail -n 1 dist${ps}.xvg | awk '{print $2}'`
echo "${ps} ${d}" >> distances.dat
rm dist${ps}.xvg
done
space=0.2 #Define the spacing of the windows
range=`cat distances.dat| awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2}; total+=$2; count+=1} END {print min" "sp" "max}' sp="$space"`
#Selecting frames for umbrella sampling input
for xx in `seq -f "%f" ${range}`;do 
aa=`awk -v c=2 -v t=${xx} '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
}' distances.dat`
grep $aa distances.dat |head -n1 >> list
done


#renaming input frames and deleting unnecessary frames
aa=1
awk '!a[$0]++' list|while read i;do num=`echo $i|awk '{print $1}'`;mv conf${num}.gro umb${aa}.gro;aa=`expr $aa + 1`;done

#Umbrella Sampling Simulations
for ii in $(seq 1 `awk '!a[$0]++' list |wc -l`);do
$GMX grompp -f npt_umbrella.mdp -c umb${ii}.gro -p topol.top -r umb${ii}.gro -n index.ndx -o npt${ii}.tpr -maxwarn 2
$GMX mdrun -deffnm npt${ii} -gpu_id 0 
$GMX grompp -f md_umbrella1.mdp -c npt${ii}.gro -t npt${ii}.cpt -p topol.top -r npt${ii}.gro -n index.ndx -o umbrella${ii}.tpr -maxwarn 2
$GMX mdrun -deffnm umbrella${ii} -gpu_id 0 
echo "umbrella${ii}.tpr" >> tpr-files.dat
echo "umbrella${ii}_pullf.xvg" >> pullf-files.dat
done
wait
$GMX wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
exit;
