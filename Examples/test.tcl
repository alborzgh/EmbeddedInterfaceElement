wipe 

model BasicBuilder -ndm 3 -ndf 3

node  1 -1 -1 -1
node  2  1 -1 -1
node  3  1  1 -1
node  4 -1  1 -1
node  5 -1 -1  1
node  6  1 -1  1
node  7  1  1  1
node  8 -1  1  1
node  9 -1 -1  3
node 10  1 -1  3
node 11  1  1  3
node 12 -1  1  3

model BasicBuilder -ndm 3 -ndf 6

node 13 0 0 -1.0
node 14 0 0 5.0


nDMaterial ElasticIsotropic 1 1.0 0.25
element stdBrick 1 1 2 3 4 5 6 7 8 1 
element stdBrick 2 5 6 7 8  9 10 11 12 1 

geomTransf Linear 1   1.0 0.0 0.0

element elasticBeamColumn 3 13 14 10.0 1.0e5 1.0e6 1.0e-2 1.0e-2 1.0e-2 1


generateInterfacePoints 3 1 -shape circle -nP 4 -nL 7 -radius 0.5 