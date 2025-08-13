[[5,5,5],[5,-5,5],[5,5,-5],[5,-5,-5],[-5,5,5],[-5,-5,5],[-5,5,-5],[-5,-5,-5]]


$HA = new HyperplaneArrangement(HYPERPLANES=>[[1,1,1],[3,-7,3],[3,3,-7],[3,-7,-7],[-7,3,3],[-7,-7,3],[-7,3,-7],[-7,-7,-7]]);

$CD=$HA->CHAMBER_DECOMPOSITION;

print $HA->CHAMBER_SIGNATURES;



//polytope

$p3 = new Polytope<Rational>(VERTICES=>[[1,5,5,5],[1,5,-5,5],[1,5,5,-5],[1,5,-5,-5],[1,-5,5,5],[1,-5,-5,5],[1,-5,5,-5],[1,-5,-5,-5]], LINEALITY_SPACE=>[]);

print $p3->VERTICES_IN_FACETS;

print $p3->GRAPH->ADJACENCY;

