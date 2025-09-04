(* SmeftFR v3.0 package *)
(* main initialization sequence, generation of mass basis Lagrangian
and Feynman rules in Mathematica format *)

Quiet[Remove["Global`*"]];
(*Off[General::stop];*)

(* FeynRules and SmeftFR package installation paths - edit if necessary *)
$FeynRulesPath = FileNameJoin[{"/home","rosiek","FeynRules"}];
SMEFT$MajorVersion      = "3";
SMEFT$MinorVersion      = "02";
SMEFT$Path = FileNameJoin[{$FeynRulesPath, "Models", "SMEFT_" <> 
                SMEFT$MajorVersion <> "_" <> SMEFT$MinorVersion}];
If[ ! DirectoryQ[SMEFT$Path], 
  Print["Directory " <> SMEFT$Path <> "does not exist, please check package setup"];
  Abort[];
];
  
(* Load FeynRules and SMEFT packages *)
Get[ FileNameJoin[{$FeynRulesPath,"FeynRules.m"}] ];
Get[ FileNameJoin[{ SMEFT$Path, "code", "smeft_package.m"}] ];

(* Choose operators to include, full list at dimension 6 and allowed list at dim8 is:

OpList6 = { "G", "Gtilde", "W", "Wtilde", "phi", "phiBox", "phiD",
 "phiW", "phiB", "phiWB", "phiWtilde", "phiBtilde", "phiWtildeB",
 "phiGtilde", "phiG", "ephi", "dphi", "uphi", "eW", "eB", "uG", "uW",
 "uB", "dG", "dW", "dB", "phil1", "phil3", "phie", "phiq1", "phiq3",
 "phiu", "phid", "phiud", "ll", "qq1", "qq3", "lq1", "lq3", "ee",
 "uu", "dd", "eu", "ed", "ud1", "ud8", "le", "lu", "ld", "qe", "qu1",
 "qu8", "qd1", "qd8", "ledq", "quqd1", "quqd8", "lequ1", "lequ3",
 "vv", "duq", "qqu", "qqq", "duu" };

OpList8 = {"phi8", "phi6Box", "phi6D2", "G2phi4n1", "G2phi4n2",
"W2phi4n1", "W2phi4n2", "W2phi4n3", "W2phi4n4", "WBphi4n1",
"WBphi4n2", "B2phi4n1", "B2phi4n2", "G4n1", "G4n2", "G4n3", "G4n4",
"G4n5", "G4n6", "G4n7", "G4n8", "G4n9", "W4n1", "W4n2", "W4n3",
"W4n4", "W4n5", "W4n6", "B4n1", "B4n2", "B4n3", "G3Bn1", "G3Bn2",
"G3Bn3", "G3Bn4", "G2W2n1", "G2W2n2", "G2W2n3", "G2W2n4", "G2W2n5",
"G2W2n6", "G2W2n7", "G2B2n1", "G2B2n2", "G2B2n3", "G2B2n4", "G2B2n5",
"G2B2n6", "G2B2n7", "W2B2n1", "W2B2n2", "W2B2n3", "W2B2n4", "W2B2n5",
"W2B2n6", "W2B2n7", "phi4D4n1", "phi4D4n2", "phi4D4n3", "G3phi2n1",
"G3phi2n2", "W3phi2n1", "W3phi2n2", "W2Bphi2n1", "W2Bphi2n2",
"G2phi2D2n1", "G2phi2D2n2", "G2phi2D2n3", "W2phi2D2n1", "W2phi2D2n2",
"W2phi2D2n3", "W2phi2D2n4", "W2phi2D2n5", "W2phi2D2n6", "WBphi2D2n1",
"WBphi2D2n2", "WBphi2D2n3", "WBphi2D2n4", "WBphi2D2n5", "WBphi2D2n6",
"B2phi2D2n1", "B2phi2D2n2", "B2phi2D2n3", "Wphi4D2n1", "Wphi4D2n2",
"Wphi4D2n3", "Wphi4D2n4", "Bphi4D2n1", "Bphi4D2n2"};
   
*)

(* Example of subset of operators (using full list leads to lengthy
calculations (may never end for full dim-8 set)!) *)

OpList6 = {"phi", "phiBox", "phiD", "phiW", "phiWB", "eB", "uW",
   "dphi", "ll"};

OpList8 = {"phi8", "phi4D4n1", "phi4D4n3"};

OpList = DeleteDuplicates[ Join[OpList6, OpList8] ];

(* file with numerical values of SMEFT parameters in WCxf JSON format,
If equal to "wcxf_default.json" or does not exist, all WCs are by
initialized to 0 *)

WCXFInput = FileNameJoin[{SMEFT$Path, "definitions", "wcxf_example.json"}];

(* initialize time counter and calculate Feynman rules, consult manual for full list of options *)
CPUTime = TimeUsed[];
SMEFTInitializeModel[Operators -> OpList,
                     Gauge -> Rxi,
                     ExpansionOrder -> 2,
                     WCXFInitFile -> WCXFInput,
                     InputScheme -> "GF",
                     CKMInput -> "yes",
                     RealParameters -> True,
                     MaxParticles -> 4];

SMEFTLoadModel[];
Print[Style["Lagrangian generation completed, time = ", Bold ], TimeUsed[] - CPUTime];

SMEFTFindMassBasis[];
Print[Style["Transformations to physical basis calculated, time = ", Bold ], TimeUsed[] - CPUTime]; 

SMEFTFeynmanRules[];
Print[Style["Feynman rules evaluated, time = ", Bold], TimeUsed[] - CPUTime];

SMEFTOutput[];
Print[Style["Output file generated, time = ", Bold], TimeUsed[] - CPUTime];

(* Mass basis vertices in different parametrizations *)

(* By default, vertices in FeynRules format are given in "none"
scheme, e.g. for Higgs-photon-photon vertex one has: *)
Print["Higgs-photon-photon vertex in \"none\" scheme: ",
      SelectVertices[GaugeHiggsVertices, SelectParticles -> {H, A, A}]];

(* For expanded vertices, one needs to choose "input scheme"; note
"Exp" extension added at the end of the name of variable storing
vertices *)

(* "smeft" input scheme *)
SMEFTExpandVertices[Input -> "smeft", ExpOrder -> 2];
Print["Higgs-photon-photon vertex in \"smeft\" scheme: ",
      SelectVertices[GaugeHiggsVerticesExp, SelectParticles -> {H, A, A}]];

(* "user" input scheme, previously chosen to be "GF" scheme *)
SMEFTExpandVertices[Input -> "user", ExpOrder -> 2];
Print["Higgs-photon-photon vertex in \"user\" scheme: ",
      SelectVertices[GaugeHiggsVerticesExp, SelectParticles -> {H, A, A}]];

(* To produce other output formats (Latex, WCxf, FeynArts, UFO) it is
now necessary to quit this session and run SmeftFR-interfaces.nb in
new kernel *)