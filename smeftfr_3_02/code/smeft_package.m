(* SmeftFR v3.0 package *)
(* loading of package and routines *)

(* load predefined variables *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_variables.m"}] ];

(* load input scheme initialization *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_input_scheme.m"}] ];

(* load model and parameter files initialization routines *)
Get[ FileNameJoin[{ SMEFT$Path, "code","smeft_initialization.m"}] ];

(* load I/O routines *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_io.m"}] ];
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_wcxf.m"}] ];

(* load functions for Lagrangian transformations *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_functions.m"}] ];

(* load functions for parameter manipulations *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_parameters.m"}] ];

(* load functions for Lagrangian rediagonalization *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_bilinears.m"}] ];

(* load functions for Feynman rules generation *)
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_gaugeint.m"}] ];
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_gaugefix.m"}] ];
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_4fermion.m"}] ];

(* load Latex generator of Feynman rules output *)
Get[ FileNameJoin[{ SMEFT$Path, "code","smeft_to_latex.m"}] ];


(* print package and setup info*)
If[ SMEFT$WelcomePrint, SMEFTPackageInfo[] ];

