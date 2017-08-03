wipe 

model BasicBuilder -ndm 3 -ndf 3

# define the domain
set x_min -5.0
set x_max  5.0
set y_min -5.0
set y_max  5.0
set z_min  0.0
set z_max  5.0

# define meshing parameters
set numEle_x 10
set numEle_y 10
set numEle_z 10

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
		
			element stdBrick [expr $numEle_y*$numEle_z*$ii + $numEle_z*$jj + $kk + 1] $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 1 
			puts $elemInfo "[expr $numEle_y*$numEle_z*$ii + $numEle_z*$jj + $kk + 1] $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 1"
		}
	}
}
close $elemInfo


# create the beam
model BasicBuilder -ndm 3 -ndf 6

node 100001 0 0 0.0
node 100002 0 0 6.0


fix 100001    1 1 1 1 1 1
# fix 100002    1 1 1 0 0 0 

# geomTransf Linear $tag $xz1 $xz2 $xz3
geomTransf Linear 1   0.0 1.0 0.0

set radius 0.5
set area   [expr 3.14159 * $radius * $radius]
set E      20.0e9
set G      [expr $E/ (2.0 * (1.0 + 0.3))]
set I      [expr $area * $area / 3.14159 / 4.0]
set J      [expr 2.0 * $I]

element elasticBeamColumn 100001 100001 100002 $area $E $G $J $I $I 1


set interfaceElems [generateInterfacePoints 100001 1 -shape circle -nP 10 -nL 10 -radius 0.5 ]


pattern Plain 1 {Series -time {0 1 1e10} -values {0 1 1} -factor 1} {
    load 100002  500000000.0  0.0    0.0   0.0 0.0 0.0
}

# eval "recorder Element -time -file contactBeam.out -ele $interfaceElems locBeam"
# eval "recorder Element -time -file contactSolid.out -ele $interfaceElems locSolid"
# eval "recorder Element -time -file contactBeamCL.out -ele $interfaceElems beamCL"
# eval "recorder Element -time -file c1.out -ele $interfaceElems c1"
# eval "recorder Element -time -file c2.out -ele $interfaceElems c2"
# eval "recorder Element -time -file c3.out -ele $interfaceElems c3"
recorder Node    -time -file beamDisp.out -node 100001 100002 -dof 1 2 3 4 5 6  disp
recorder Node    -time -file displacement.out -nodeRange 1 [expr ($numEle_y+1)*($numEle_z+1)*($numEle_z+1)]  -dof 1 2 3 disp


# Create analysis
constraints Transformation
test        NormDispIncr 1.0e-5 20 1
algorithm   Newton
numberer    RCM
system      BandGeneral
integrator  LoadControl 0.01
analysis    Static

analyze 100

wipe