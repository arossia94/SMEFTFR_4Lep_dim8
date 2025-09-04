(* SmeftFR v3.0 package *)
(* routines for initialization of the model in weak and mass basis *)


SMEFTPackageInfo = Function[{},
(* print package information to screen *)

Print["*****************************************"];
Print[" - SmeftFR package - "];
Print["Version:    ", SMEFT$Version, " ("SMEFT$Date, ")"];
Print["Authors:    A. Dedes, J. Rosiek, M. Ryczkowski, K. Suxho, L. Trifyllis"];
Print["URL:        www.fuw.edu.pl/smeft"];
Print["Maintainer: janusz.rosiek@fuw.edu.pl"];
Print[" "];
Print["Please cite:"];
Print["    - JHEP 1706 (2017) 143 (1704.03888)"];
Print["    - Comput.Phys.Commun. 247 (2020) 106931, (1904.03204)"];
Print["    - arXiv:2302.01353"];
Print["*****************************************"];
Print[" "];

(* end of SMEFTPackageInfo *)
]



SMEFTInitializeModel[ OptionsPattern[{ Operators -> SMEFT$Dim6Operators,
                         Gauge -> "Unitary",
                         ExpansionOrder -> 1,
                         WCXFInitFile -> "",
			 RealParameters -> False,
			 InputScheme -> "GF",	       
			 CKMInput -> "yes",	       
                         MajoranaNeutrino -> False,
                         WBFirstLetter -> "c",
                         MBFirstLetter -> "C",
			 Correct4Fermion -> False,
			 MaxParticles -> 6 }] ] :=
Module[{},
(* routines generating and loading FeynRules parameter files and SMEFT Lagrangian in Weak Basis mode *)

If [ Complement[ {OptionValue[RealParameters], OptionValue[MajoranaNeutrino], OptionValue[Correct4Fermion]},
                 {True,False} ] != {},
   Print["SMEFTInitializeModel: options RealParameters, MajoranaNeutrino, Correct4Fermion can be only True or False!"];
   Abort[];
];

If [ Head[ OptionValue[InputScheme] ] =!= String, 
  Print["Option InputScheme must be of string type, check selected SMEFTInitializeModel options"];
  Abort[];
];

       
If [ ! MemberQ[ {"no","yes","force"}, OptionValue[CKMInput] ],
  Print["SMEFTInitializeModel: option CKMInput can be only \"no\", \"yes\" or \"force\"!"];
  Abort[];
];

(* naming convention - editable first letter of WC names in Warsaw and
mass basis. Defaults are c and C, respectively: *)
SMEFT$WB = OptionValue[WBFirstLetter];
If [ Head[ SMEFT$WB ] === String && StringLength[ SMEFT$WB ] > 0,
  SMEFT$WB = StringTake[SMEFT$WB,1];
,
  Print["Non-string value for WBFirstLetter, check selected SMEFTInitializeModel options"];
  Abort[];
];
SMEFT$MB = OptionValue[MBFirstLetter];
If [ Head[ SMEFT$MB ] === String && StringLength[ SMEFT$MB ] > 0,
  SMEFT$MB = StringTake[SMEFT$MB,1];
,
  Print["Non-string value for MBFirstLetter, check selected SMEFTInitializeModel options"];
  Abort[];
];

(* decides if to include corrections from WCs to CKM elements even if they exceed
SMEFT$CKMTreshold value *)
SMEFT$CKMInput = OptionValue[CKMInput];

(* maximal number of legs in vertices *)
SMEFT$MaxParticles = OptionValue[MaxParticles];

(* use complex or real only parameters - MadGraph5 requires real couplings only! *)
SMEFT$RealParameters = OptionValue[RealParameters];       
       
(* active operator list in string format and useful operator sublists *)
SMEFT$OperatorList = Intersection[DeleteDuplicates[ ToString /@ OptionValue[Operators] ], SMEFT$AllOperators];
GenerateOperatorLists[];
              
(* maximal expansion order is 1/Lambda^(2*ExpansionOrder) *)
SMEFT$ExpansionOrder = OptionValue[ExpansionOrder];
If [ ! MemberQ[ {0,1,2}, SMEFT$ExpansionOrder ],
  Print["ExpansionOrder can be only 0, 1 or 2, check selected SMEFTInitializeModel options"];
  Abort[];
];
    
(* gauge selection, rxi or unitary *)
SMEFT$RxiGaugeStatus = TrueQ[ ToLowerCase[ToString[OptionValue[Gauge]]] === "rxi" ];

(* status of neutrino fields, Weyl or Majorana *)
SMEFT$MajoranaNeutrino = MemberQ[SMEFT$OperatorList, "vv"] || OptionValue[MajoranaNeutrino];

SMEFT$Correct4FermionSign = OptionValue[Correct4Fermion];

SMEFT$WCXF = ToString[ OptionValue[WCXFInitFile] ];
If [ SMEFT$WCXF != "",
  If[ ! FileExistsQ[SMEFT$WCXF],
    Print["File " <> SMEFT$WCXF <> " does not exist, Wilson coefficient values will be initialized to 0"];
    SMEFT$WCXF = FileNameJoin[{SMEFT$Path, "definitions", "wcxf_default.json"}];
  ];
,
  SMEFT$WCXF = FileNameJoin[{SMEFT$Path, "definitions", "wcxf_default.json"}];
];

(* check if user defined input parameters name in the header of
smeft_input_scheme.m does not overlap with standard SmeftFR parameters
names *)
CheckClashingNames[];

(* print setup variables *)
Print["Operators included in Feynman rules: ", Style[ SMEFT$OperatorList, Bold ] ];
Print["Feynman rules will be calculated for ",
  Style[ If[SMEFT$RxiGaugeStatus,"R_xi","unitary"] <> " gauge", Bold ] ];
Print["Neutrinos are treated as ",
      Style[If[SMEFT$MajoranaNeutrino, "Majorana", "massless Weyl"], Bold ], " spinors\n" ];

(* generate parameter files in the Warsaw basis *)
SMEFTParFile[SMEFT$OperatorList, FileNameJoin[{SMEFT$Path, "output", "smeft_par_WB.fr"}] ];

(* read numerical input date for WCs from the file in the WCxf format *)
ReadWCXFInput[ SMEFT$WCXF, Operators -> OptionValue[Operators], Silent -> True ];
		
SMEFTInputScheme[ InputScheme -> OptionValue[InputScheme] ];

SMEFT$UserInputInitialization = True;
    
UpdateSMParameters[];
       
(* end of SMEFTInitializeModel *)
];






SMEFTInitializeMB[ OptionsPattern[{
    ModelFile -> FileNameJoin[{SMEFT$Path, "output", "smeft_par_MB.fr"}],
    InteractionFile -> FileNameJoin[{SMEFT$Path, "output", "smeft_feynman_rules.m"}],
    Expansion -> "smeft",
    Include4Fermion -> True,
    IncludeBL4Fermion -> False }] ] :=

(* Load mass basis parameter files and Lagrangian, update parameter
   values *)
Module[{ism, jsm, tmp, L2Fermion, LEW, LQCD, L4Fermion, LBL4Fermion,
LSMEFT, mlist},

mlist = {"mz", "mw" , "mh", "mqt", "mqc", "mqu", "mqb", "mqs", "mqd",
  "mlt", "mlm", "mle", "mvt", "mvm", "mve"};

If [ ! FileExistsQ[ OptionValue[InteractionFile] ],
  Print["File " <> OptionValue[InteractionFile] <> " does not exist, check its name and location"];
  Abort[];
];
       
If [ ! FileExistsQ[ OptionValue[ModelFile] ],
  Print["File " <> OptionValue[ModelFile] <> " does not exist, check its name and location"];
  Abort[];
];
       
(* Load SM fields and parameter definitions in mass basis*)
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_fields_MB.fr"}] ];
Get[ OptionValue[ModelFile] ];
LoadModel[];

(* load MB Lagrangian *)
Get[ OptionValue[InteractionFile] ];

(* useful active operator sub-lists *)
GenerateOperatorLists[];

(* Give warning if lequ3 operator is on the list *)
(*If [ MemberQ[ SMEFT$OperatorList, "lequ3" ],
  Print[ Style["\nWARNING: ", Bold, Red ], Style["4-fermion operator \"lequ3\" is included in SMEFT Lagrangian! It contains tensor coupling sigma_munu, not compatible with Madgraph 5.", Red] ];

  Print[Style["\nTo avoid Madgraph crashes, remove \"lequ3\" from operator list or keep it but edit manually UFO files generated by SmeftFR, replacing sigma_munu by explicit product of gamma matrices."  ,Red] ];
];*)

(* update masses of all particles *)
UpdateParameters[ ToExpression[ ToUpperCase[#] ] -> NumericalValue[ ToExpression[#] ] ] & /@ mlist;

(* update fermion mass parameters *)
UpdateParameters[ fml[1,1] -> NumericalValue[mle] ];
UpdateParameters[ fml[2,2] -> NumericalValue[mlm] ];
UpdateParameters[ fml[3,3] -> NumericalValue[mlt] ];
UpdateParameters[ fmu[1,1] -> NumericalValue[mqu] ];
UpdateParameters[ fmu[2,2] -> NumericalValue[mqc] ];
UpdateParameters[ fmu[3,3] -> NumericalValue[mqt] ];
UpdateParameters[ fmd[1,1] -> NumericalValue[mqd] ];
UpdateParameters[ fmd[2,2] -> NumericalValue[mqs] ];
UpdateParameters[ fmd[3,3] -> NumericalValue[mqb] ];

(* update ghost and Goldstone masses in Rxi gauge, unused in unitary gauge *)
UpdateParameters[ MgZ -> Sqrt[NumericalValue[xiZ]] NumericalValue[MZ],
                  MG0 -> Sqrt[NumericalValue[xiZ]] NumericalValue[MZ] ];
UpdateParameters[ MgW -> Sqrt[NumericalValue[xiW]] NumericalValue[MW],
                  MGP -> Sqrt[NumericalValue[xiW]] NumericalValue[MW] ];

UpdateParameters[ Lam -> 1 ];

NormalizationConstantsExpansion[ OptionValue[Expansion] ];

L2Fermion = OutputParametersExpansion[LeptonGaugeLagrangian] +
       OutputParametersExpansion[LeptonHiggsGaugeLagrangian] +
       OutputParametersExpansion[DeltaLTwoLagrangian] +
       OutputParametersExpansion[QuarkGaugeLagrangian] +
       OutputParametersExpansion[QuarkHiggsGaugeLagrangian]; 

LEW  = OutputParametersExpansion[GaugeSelfLagrangian] +
       OutputParametersExpansion[GaugeHiggsLagrangian];
LQCD = OutputParametersExpansion[QuarkGluonLagrangian] +
       OutputParametersExpansion[GluonSelfLagrangian] +
       OutputParametersExpansion[GluonHiggsLagrangian];

LSMEFT = L2Fermion + LEW + LQCD + If[ SMEFT$RxiGaugeStatus, OutputParametersExpansion[GhostLagrangian], 0];

(* 4-fermion vertices may not be correct due to relative vertex sign problems!  *)
If [ OptionValue[Include4Fermion] === True,
  L4Fermion = OutputParametersExpansion[FourLeptonLagrangian] +
              OutputParametersExpansion[TwoQuarkTwoLeptonLagrangian] +
              OutputParametersExpansion[FourQuarkLagrangian];
  LSMEFT = LSMEFT + L4Fermion; 
];
(* BL-violating 4-fermion vertices are even more problematic with FeynRules!  *)
If [ OptionValue[IncludeBL4Fermion] === True,
  LBL4Fermion = OutputParametersExpansion[BLViolatingLagrangian];
  LSMEFT = LSMEFT + LBL4Fermion; 
];

(* replace photon-Z mixing by non-matrix notation, irrelevant if not Expansion->"none" *)

LSMEFT = LSMEFT /. AZnorm[ Index[SU2W,1], Index[SU2W,1] ] -> AZnorm11 /.
                   AZnorm[ Index[SU2W,1], Index[SU2W,2] ] -> AZnorm12 /.
                   AZnorm[ Index[SU2W,2], Index[SU2W,1] ] -> AZnorm21 /.
                   AZnorm[ Index[SU2W,2], Index[SU2W,2] ] -> AZnorm22;

(* scale included in WC's, so remember that Lam = 1 in UFO parameter file! *)

SMEFT$MBLagrangian = LSMEFT;

]
(* end of SMEFTInitializeMB *)





SMEFTLoadModel = Function[{},
(* calculate SMEFT Lagrangian parts - can be time-consuming *)
Module[{},

SMEFTLoadModelFiles[ "smeft_par_WB.fr" ];
SMEFTLoadLagrangian[];

SMEFT$LGeff = LGauge + LHiggs;
SMEFT$LGeff6 = 0;
SMEFT$LGFmass = LFermions + LYukawa // FunctionExpand;
ttt = SMEFT$LGFmass;
If[ SMEFT$ExpansionOrder > 0,
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGeff = SMEFT$LGeff + Lam ToExpression["LQ"<>#] ] & /@ GaugeBilinearOperators;
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGeff6 = SMEFT$LGeff6 + Lam ToExpression["LQ"<>#] ] & /@ GaugeTripleOperators;
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGFmass = SMEFT$LGFmass + Lam ToExpression["LQ"<>#] ] & /@ TwoFermionMassOperators;
];
If[ SMEFT$ExpansionOrder > 1,
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGeff = SMEFT$LGeff + Lam^2 ToExpression["LQ"<>#] ] & /@ GaugeBilinearOperators8;
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGeff6 = SMEFT$LGeff6 + Lam^2 ToExpression["LQ"<>#] ] & /@ GaugeQuadrupleOperators;
];

SMEFT$LGferm = SMEFT$LGFmass;
If[ SMEFT$ExpansionOrder > 0,
  If[ MemberQ[ SMEFT$OperatorList, # ], SMEFT$LGferm = SMEFT$LGferm + Lam ToExpression["LQ"<>#] ] & /@ TwoFermionOperators;
];


];
(* end of SMEFTLoadModel *)
]

