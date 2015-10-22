## Movie callback for transparent changed during simulation.

proc moviecallback { args } {
 if {$::MovieMaker::userframe == 250} {
     animate speed 0.70
 }
 if {$::MovieMaker::userframe == 350} {
     animate speed 1.0
 }
 if {$::MovieMaker::userframe > 250 && $::MovieMaker::userframe <= 350} {
     material change opacity Material22 [expr 0.75-(0.0075*($::MovieMaker::userframe-250))]
 }
 if {$::MovieMaker::userframe > 300 && $::MovieMaker::userframe <= 400} {
     material change opacity Material23 [expr 1.0-(0.01*($::MovieMaker::userframe-300))]
     material change opacity Material24 [expr 0.0+(0.01*($::MovieMaker::userframe-300))]
 }

 if {$::MovieMaker::userframe > 100 && $::MovieMaker::userframe <= 200} {
     material change opacity Material22 [expr 1.0-(0.0025*($::MovieMaker::userframe-100))]
 }
 animate next
}

proc setupmovie {} {
 enablemoviecallback
}

proc finishmovie {} {
  disablemoviecallback
}


## Easy-to-use proc to enable the user-defined movie frame callback
proc enablemoviecallback { }  {
  trace add variable ::MovieMaker::userframe write moviecallback
}

## Easy-to-use proc to disable the user-defined movie frame callback
proc disablemoviecallback { } {
  trace remove variable ::MovieMaker::userframe write moviecallback
}
