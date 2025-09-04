(* SmeftFR v3.0 package *)
(* SmeftFR interfaces - run smeft_fr_init.m prior to this code! *)

Quiet[Remove["Global`*"]];
(* Off[General::stop]; *)

(* FeynRules and SmeftFR package installation paths - edit if necessary *)
$FeynRulesPath = FileNameJoin[{"/home","rosiek","FeynRules"}];
SMEFT$MajorVersion      = "3";
SMEFT$MinorVersion      = "01";
SMEFT$Path = FileNameJoin[{$FeynRulesPath, "Models", "SMEFT_" <> 
                SMEFT$MajorVersion <> "_" <> SMEFT$MinorVersion}];
If[ ! DirectoryQ[SMEFT$Path], 
  Print["Directory " <> SMEFT$Path <> "does not exist, please check package setup"];
  Abort[];
];

(* Load FeynRules and SMEFT packages *)
Get[ FileNameJoin[{$FeynRulesPath,"FeynRules.m"}] ];
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_package.m"}] ];

(* initialize time counter *)     
CPUTime = TimeUsed[];

(* initialize mass basis model files and Lagrangian
Expansion option can be "none", "smeft" or "user" *)
SMEFTInitializeMB[ Expansion->"user",
		   Include4Fermion->True,
		   IncludeBL4Fermion->False
		 ];
Print["\nMass basis Lagrangian loaded, time = ", TimeUsed[] - CPUTime,"\n"];

(* WCxf output *)
SMEFTToWCXF[ FileNameJoin[ { SMEFT$Path, "output", "smeft_par_MB.fr" } ],
             FileNameJoin[ { SMEFT$Path, "output", "smeft_wcxf_MB.json" } ] ];
Print["\nParameters stored in WCXF file, time = ", TimeUsed[] - CPUTime];

(* Latex output *)
SMEFTToLatex[ Expansion -> "smeft" (*, ScreenOutput->True*) ];
Print["\nLatex output generated, time = ", TimeUsed[] - CPUTime];  

(* UFO output *)
SMEFTToUFO[ SMEFT$MBLagrangian,
	    Output -> FileNameJoin[{SMEFT$Path, "output", "UFO"}],
	    CorrectIO -> True ];
Print["\nUFO output generated, time = ", TimeUsed[] - CPUTime,"\n"];

(* FeynArts output *)
WriteFeynArtsOutput[ SMEFT$MBLagrangian,
		     Output -> FileNameJoin[{SMEFT$Path, "output", "FeynArts", "FeynArts"}] ];
Print["\nVertices stored in FeynArts format in directory ", 
       Style[FileNameJoin[{SMEFT$Path, "output", "FeynArts"}] ,Bold] ];
Print["\nTotal CPU time used = ", TimeUsed[] - CPUTime];



