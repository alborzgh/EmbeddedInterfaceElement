wipe 

model BasicBuilder -ndm 3 -ndf 3

# set the flag for using SSP elements ( 0 : standard brick - 1 : SSP brick)
set useSSP 1

# define the domain
set x_min -5.0
set x_max  5.0
set y_min -5.0
set y_max  5.0
set z_min  0.0
set z_max  20.0

# define meshing parameters
set numEle_x 10
set numEle_y 10
set numEle_z 20

# define beam location
set b_x1  0.0
set b_y1  0.0
set b_z1  10.0
set b_x2  0.0
set b_y2  0.0
set b_z2  20.0

# define beam meshing parameter
set b_numEle 10

# define number of descritization points
set nP 8
set nL 8

# create the soil geometry
set dx [expr ($x_max-$x_min)/$numEle_x]
set dy [expr ($y_max-$y_min)/$numEle_y]
set dz [expr ($z_max-$z_min)/$numEle_z]

set nodeInfo [open nodeInfo.dat w]
for {set ii 0} {$ii <= $numEle_x} {incr ii} {
	for {set jj 0} {$jj <= $numEle_y} {incr jj} {
		for {set kk 0} {$kk <= $numEle_z} {incr kk}	{
			node [expr ($numEle_y+1)*($numEle_z+1)*$ii + ($numEle_z+1)*$jj + $kk + 1] [expr $x_min + $ii * $dx] [expr $y_min + $jj  * $dy] [expr $z_min + $kk * $dz]
			puts $nodeInfo "[expr ($numEle_y+1)*($numEle_z+1)*$ii + ($numEle_z+1)*$jj + $kk + 1] [expr $x_min + $ii * $dx] [expr $y_min + $jj  * $dy] [expr $z_min + $kk * $dz]"
		}
	}
}
close $nodeInfo

# fix the boundaries of the mesh
fixX $x_min 1 1 1
fixX $x_max 1 1 1
fixY $y_min 1 1 1
fixY $y_max 1 1 1
fixZ $z_min 1 1 1

# create the soil material and elements
#          ElasticIsotropic $tag $E $poisson
nDMaterial ElasticIsotropic 1 1e8 0.33

set elemInfo [open elementInfo.dat w]
for {set ii 0} {$ii < $numEle_x} {incr ii} {
	for {set jj 0} {$jj < $numEle_y} {incr jj} {
		for {set kk 0} {$kk < $numEle_z} {incr kk}	{
			set node1 [expr ($numEle_y+1)*($numEle_z+1)*$ii + ($numEle_z+1)*$jj + $kk + 1]
			set node2 [expr ($numEle_y+1)*($numEle_z+1)*($ii+1) + ($numEle_z+1)*$jj + $kk + 1]
			set node3 [expr ($numEle_y+1)*($numEle_z+1)*($ii+1) + ($numEle_z+1)*($jj+1) + $kk + 1]
			set node4 [expr ($numEle_y+1)*($numEle_z+1)*$ii + ($numEle_z+1)*($jj+1) + $kk + 1]
			set node5 [expr $node1 + 1]
			set node6 [expr $node2 + 1]
			set node7 [expr $node3 + 1]
			set node8 [expr $node4 + 1]
		
			if {$useSSP == 1} {
				element SSPbrick [expr $numEle_y*$numEle_z*$ii + $numEle_z*$jj + $kk + 1] $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 1 
			} else {
				element stdBrick [expr $numEle_y*$numEle_z*$ii + $numEle_z*$jj + $kk + 1] $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 1 
			}
			puts $elemInfo "[expr $numEle_y*$numEle_z*$ii + $numEle_z*$jj + $kk + 1] $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 1"
		}
	}
}
close $elemInfo


# create the beam
model BasicBuilder -ndm 3 -ndf 6

set b_dx [expr ($b_x2 - $b_x1)/($b_numEle + 0.0)]
set b_dy [expr ($b_y2 - $b_y1)/($b_numEle + 0.0)]
set b_dz [expr ($b_z2 - $b_z1)/($b_numEle + 0.0)]

set nodeInfo [open nodeInfoB.dat w]
for {set ii 0} {$ii <= $b_numEle} {incr ii} {
			node [expr 100000 + $ii + 1] [expr $b_x1 + $ii * $b_dx] [expr $b_y1 + $ii * $b_dy] [expr $b_z1 + $ii * $b_dz]
			puts $nodeInfo "[expr 100000 + $ii + 1] [expr $b_x1 + $ii * $b_dx] [expr $b_y1 + $ii * $b_dy] [expr $b_z1 + $ii * $b_dz]"
}
close $nodeInfo

# geomTransf Linear $tag $xz1 $xz2 $xz3
geomTransf Linear 1   0.0 1.0 0.0

set radius 0.25
set area   [expr 3.14159 * $radius * $radius]
set E      3.0e10
set G      [expr $E/ (2.0 * (1.0 + 0.3))]
set I      [expr $area * $area / 3.14159 / 4.0]
set J      [expr 2.0 * $I]

set elemInfo [open elementInfoB.dat w]
for {set ii 0} {$ii < $b_numEle} {incr ii}	{
	set node1 [expr 100000 + $ii + 1]
	set node2 [expr $node1 + 1]

	element elasticBeamColumn $node1 $node1 $node2 $area $E $G $J $I $I 1
	puts $elemInfo "$node1 $node1 $node2 1"
}
close $elemInfo

set interfaceElemsFile "interfaceInfo.dat"
if {[file exists $interfaceElemsFile] == 1} { file delete $interfaceElemsFile }

# fix 100001   1 1 1 1 1 1

set interfaceElems {}
for {set ii 0} {$ii < $b_numEle} {incr ii}	{
	set elem [expr 100000 + $ii + 1]
	set interfaceElems [concat $interfaceElems [generateInterfacePoints $elem 1  -gPenalty -shape circle -nP $nP -nL $nL -radius $radius -file $interfaceElemsFile ]]
}
set interfaceElems [concat $interfaceElems [generateToeInterfacePoints 100001 1 100001 -gPenalty -shape circle -nP 4 -nR 2 -radius $radius -file $interfaceElemsFile ]]

pattern Plain 1 {Series -time {0 1 1e10} -values {0 1 1} -factor 1} {
    load [expr 100000 + $b_numEle + 1]  0.0  0.0    -100000.0   0.0  0.0   0.0
}

eval "recorder Element -time -file contactBeam.out -ele $interfaceElems locBeam"
eval "recorder Element -time -file contactDisp.out -ele $interfaceElems disp"
eval "recorder Element -time -file contactBeamCL.out -ele $interfaceElems beamCL"
eval "recorder Element -time -file c1.out -ele $interfaceElems c1"
eval "recorder Element -time -file c2.out -ele $interfaceElems c2"
eval "recorder Element -time -file c3.out -ele $interfaceElems c3"
recorder Node    -time -file beamDisp.out -nodeRange 100001 [expr 100001 + $b_numEle] -dof 1 2 3 4 5 6  disp
recorder Node    -time -file displacement.out -nodeRange 1 [expr ($numEle_y+1)*($numEle_z+1)*($numEle_z+1)]  -dof 1 2 3 disp

if {$useSSP == 1} {
	recorder Element -time -file stress.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] stress
	recorder Element -time -file strain.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] strain
} else {
	recorder Element -time -file stress1.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 1 stress
	recorder Element -time -file strain1.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 1 strain
	recorder Element -time -file stress2.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 2 stress
	recorder Element -time -file strain2.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 2 strain
	recorder Element -time -file stress3.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 3 stress
	recorder Element -time -file strain3.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 3 strain
	recorder Element -time -file stress4.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 4 stress
	recorder Element -time -file strain4.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 4 strain
	recorder Element -time -file stress5.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 5 stress
	recorder Element -time -file strain5.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 5 strain
	recorder Element -time -file stress6.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 6 stress
	recorder Element -time -file strain6.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 6 strain
	recorder Element -time -file stress7.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 7 stress
	recorder Element -time -file strain7.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 7 strain
	recorder Element -time -file stress8.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 8 stress
	recorder Element -time -file strain8.out -eleRange 1 [expr $numEle_x*$numEle_y*$numEle_z] material 8 strain
}

set dt 0.1
set numSteps 10

# Create analysis
constraints Transformation
# constraints Penalty 1.0e15 1.0e15
# test        NormDispIncr 1.0e-10 20 1
test EnergyIncr 1.0e-15 20 1
# test NormUnbalance 1.0e-11 20 1
# test FixedNumIter 2 1
algorithm   Newton
numberer    RCM
system      Mumps
# system      BandGeneral
integrator  LoadControl $dt
analysis    Static

analyze $numSteps

wipe

set batchFile [open writeGiD.bat w]
# puts $batchFile "GiDFlaviaWriter -solidNodeInfo nodeInfo.dat -solidElemInfo elementInfo.dat -beamNodeInfo nodeInfoB.dat -beamElemInfo elementInfoB.dat -nSolidNodes [expr ($numEle_x + 1)*($numEle_y + 1)*($numEle_z + 1)] -nSolidElems [expr $numEle_x*$numEle_y*$numEle_z] -nBeamNodes [expr $b_numEle+1] -nBeamElems $b_numEle -output res -binary -nTimeSteps $numSteps -solidDispInfo displacement.out -beamDispInfo beamDisp.out -interfacePtInfo interfaceInfo.dat -nInterfacePts [expr $b_numEle * $nP * $nL] -intPtDispInfo contactDisp.out"
puts $batchFile "GiDFlaviaWriter -solidNodeInfo nodeInfo.dat -solidElemInfo elementInfo.dat -beamNodeInfo nodeInfoB.dat -beamElemInfo elementInfoB.dat -nSolidNodes [expr ($numEle_x + 1)*($numEle_y + 1)*($numEle_z + 1)] -nSolidElems [expr $numEle_x*$numEle_y*$numEle_z] -nBeamNodes [expr $b_numEle+1] -nBeamElems $b_numEle -output res -binary -nTimeSteps $numSteps -solidDispInfo displacement.out -beamDispInfo beamDisp.out -interfacePtInfo interfaceInfo.dat -nInterfacePts 6 -intPtDispInfo contactDisp.out"
close $batchFile
