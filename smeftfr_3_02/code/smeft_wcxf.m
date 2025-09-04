(* SmeftFR v3.0 package *)
(* routines for WCxf input/output operations *)


(* make complex WC out of Re and Im in JSON file, scientific notation *)
SMEFTToComplex[x_] := If[NumberQ[x], x, ("Re" + I "Im") /. x ];



SMEFTToWCXF[ SourceFile_, TargetFile_,
	     OptionsPattern[ { FirstLetter -> SMEFT$MB } ] ] :=
(* WCxf JSON export function - reads from SMEFT model file and
   generates Wilson coefficients file in JSON format *)

Module[{il, it, WC, WN, ig1, ig2, ig3, ig4, tmp, date, values, output,
oplist, stdl, first, temp, mspar, ActiveS, Active2T, Active4T},

stdl = {"G1", "GW", "GS", "Gf", "vev", "hlambda", "MZ", "MW" , "MH",
  "MQT", "MQC", "MQU", "MQB", "MQS", "MQD", "MLT", "MLM", "MLE",
  "MVT", "MVM", "MVE", "mz", "mw" , "mh", "mqt", "mqc", "mqu", "mqb",
  "mqs", "mqd", "mlt", "mlm", "mle", "mvt", "mvm", "mve", "K", "U",
  "G0norm", "g1norm", "Gnorm", "GPnorm", "gsnorm", "gwnorm", "Hnorm",
  "AZnorm", "Kq", "Lam", "MG0", "MGP", "MgW", "MgZ", "Ul", "Wnorm",
  "xiA", "xiG", "xiW", "xiZ" };

If [ DeleteDuplicates[ (Head[#] & /@ {SourceFile, TargetFile} ) ] != {String},
  Print["File names in SMEFTToWCXF are not of string type, aborting..."];
  Abort[];
];

first = ToString[ OptionValue[FirstLetter] ];       
(* first letter in WCs definition *)
If [ StringLength[first] != 1,
  Print["Option FirstLetter in SMEFTToWCXF is not a single character, aborting..."];
  Abort[];
];
	
If [ ! FileExistsQ[ SourceFile ],
  Print["File " <> SourceFile <> " does not exist, check its name and location"];
  Abort[];
];

Get[ SourceFile ];

(* find parameters with "Value" assigned *)
mspar = {};
temp = Value /. Values[ Rule @@@ M$Parameters ];
For[il=1, il < Length[temp] + 1, il++,
  If[ temp[[il]] =!= Value, AppendTo[ mspar, (Rule @@@ M$Parameters)[[il]] ] ];
];

(* find dim-5 and 6 WC names only, skip other parameters *)
temp = Complement[  ToString /@ Keys[mspar], stdl ];
oplist = Intersection[ SMEFT$AllOperators, StringDrop[ temp, 1 ] ];

ActiveS = Intersection[ToString /@ ScalarWC, oplist];
ScalarVal = If[ Length[ActiveS] > 0, Value /. (ToExpression[StringInsert[ActiveS, first, 1]] /. mspar), 0 ];

Active2T = Intersection[Tensor2WC[[All,1]], oplist];
Tensor2Val = If[ Length[Active2T] > 0, Association[ Flatten[ Value /. (ToExpression[StringInsert[Active2T, first, 1]] /. mspar) ] ], 0 ];

Active4T = Intersection[Tensor4WC[[All,1]], oplist];
Tensor4Val = If[ Length[Active4T] > 0, Association[ Flatten[ Value /. (ToExpression[StringInsert[Active4T, first, 1]] /. mspar) ] ], {} ];

Print["Wilson coefficients in WCXF format stored as ", Style[ TargetFile, Bold ] ];

output = {};

date = Date[];
AppendTo[ output, "eft" -> "SMEFT" ];
AppendTo[ output, "basis" -> "Warsaw mass" ];
AppendTo[ output, "metadata" -> { "generator" -> " SmeftFR v" <> SMEFT$Version <> ", SMEFTToWCXF routine",
			          "description" -> "Generated from the FeynRules model file " <>
	                           SourceFile <> " on " <> ToString[date[[3]]] <> "." <>
	                           ToString[date[[2]]] <> "." <> ToString[date[[1]]] <> " " <>
	                           ToString[date[[4]]] <> ":" <> ToString[date[[5]]] } ];

(* renormalization scale for WCs currently not defined in SMEFT model files

temp = If [ MemberQ[ Keys[mspar], Lam ], Sqrt[Value /. (Lam /. mspar)], Sqrt[Value] ];
AppendTo[ output, "scale" -> If[ NumberQ[temp], temp, 0 ] ];

*)
AppendTo[ output, "scale" -> 0 ];

values = {};

For[ il=1, il < Length[ActiveS] + 1, il++,
  If[ NumberQ[ ScalarVal[[il]]], AppendTo[ values, (ActiveS[[il]] /. ScalarJSON) -> ScalarVal[[il]] ] ];
];

Active2T = Complement[Active2T, {"vv"}];
For[ il=1, il < Length[Active2T] + 1, il++,
  tmp = Tensor2Ind[ Active2T[[il]] ];
  For[ it=1, it < Length[tmp] + 1, it++,
    ig1 = tmp[[it,1]];
    ig2 = tmp[[it,2]];
    WN = (Active2T[[il]] /. Tensor2JSON) <> "_" <> ToString[ig1] <> ToString[ig2];
    WC = Tensor2Val[ ToExpression[ first <> Active2T[[il]] ] @@ {ig1,ig2} ];
    If[ NumberQ[WC],
      If[ tmp[[it,3]] && NumberQ[WC],
        AppendTo[ values, WN ->  WC];
      ,	
        AppendTo[ values, WN -> {"Re" -> Re[WC], "Im" -> Im[WC]} ];
      ];
    ];
  ];
];

For[ il=1, il < Length[Active4T] + 1, il++,
  tmp = Tensor4Ind[ Active4T[[il]] ];
  For[ it=1, it < Length[tmp] + 1, it++,
    ig1 = tmp[[it,1]];
    ig2 = tmp[[it,2]];
    ig3 = tmp[[it,3]];
    ig4 = tmp[[it,4]];
    WN = (Active4T[[il]] /. Tensor4JSON) <> "_" <> ToString[ig1] <> ToString[ig2] <> ToString[ig3] <> ToString[ig4];
    WC = Tensor4Val[ ToExpression[ first <> Active4T[[il]] ] @@ {ig1,ig2,ig3,ig4} ];
    If[ NumberQ[WC],
      If[ tmp[[it,5]],
        AppendTo[ values, WN ->  WC];
      ,	
        AppendTo[ values, WN -> {"Re" -> Re[WC], "Im" -> Im[WC]} ];
      ];
    ];
  ];
];

(* special treatment for vv operator *)
If[ MemberQ[oplist, "vv"],
  tmp = Tensor2Ind[ "vv" ];
  For[ it=1, it < Length[tmp] + 1, it++,
    ig1 = tmp[[it,1]];
    ig2 = tmp[[it,2]];
    WN = ("vv" /. Tensor2JSON) <> "_" <> ToString[ig1] <> ToString[ig2];
    WC = Tensor2Val[ ToExpression[ first <> "vv"] @@ {ig1,ig2} ];
    AppendTo[ values, WN -> {"Re" -> Re[WC], "Im" -> Im[WC]} ];
  ];
];

AppendTo[output, "values" -> values ];
Export[TargetFile, output, "json"];

]
(* end of SMEFTToWCXF *)





ReadWCXFInput[ SourceFile_, OptionsPattern[ { Operators -> SMEFT$AllOperators,
                                              Silent -> False} ] ] :=

(* WCxf JSON import function - reads from JSON format and generates entries describing Wilson coefficients *)
Module[{JSONData,JSONWilson,JSONNum,il,ig1,ig2,ig3,ig4,tmp,WC,ActiveS,Active2T,Active4T,oplist,ind,ind1,wcname},

If[ Head[SourceFile] =!= String, Print["Source file name must be string!"]; Abort[]; ];

If [ ! FileExistsQ[ SourceFile ],
  Print["File " <> SourceFile <> " does not exist, check its name and location"];
  Abort[];
];

If [ FileFormat[ SourceFile ] != "JSON",
  Print["File " <> SourceFile <> " not in JSON format or does not have extension .json or .JSON, please correct!"];
  Abort[];
];

oplist = OptionValue[Operators];
If[ Head[oplist] =!= List, Print["Option 'Operators' must be a list, aborting!"];Abort[];];

If [ OptionValue[Silent] =!= True, Print["Wilson coefficient values initialized from WCxf file ", Style[ SourceFile, Bold ] ] ];
JSONData = Import[ SourceFile ];
JSONWilson = SMEFTToComplex /@ Association["values" /. JSONData];
JSONNum[tmp_] := If[ NumberQ[tmp], tmp, 0];

ActiveS = ToString /@ Intersection[ScalarWC, oplist];
ScalarVal = Table[ToExpression[SMEFT$MB <> ActiveS[[il]]] -> JSONNum[ JSONWilson[ ActiveS[[il]] /. ScalarJSON ] ], {il,1,Length[ActiveS]} ];

ind = Flatten[ Table[wcname[i,j], {i,1,3}, {j,1,3}] ];
ind1 = Table["_"<>ToString[i]<>ToString[j], {i,1,3}, {j,1,3}];
Active2T = ToString /@ Intersection[ Tensor2WC[[All,1]], oplist ];
Tensor2Val = {};
For[ il = 1, il < Length[Active2T] + 1, il++,
  WC = Table[ JSONNum[ JSONWilson[ (Active2T[[il]] /. Tensor2JSON) <> ind1[[ig1,ig2]] ] ], {ig1,1,3}, {ig2,1,3} ];
  WC = Flatten[ Tensor2WCSymmetrize[ WC, ToExpression[ Active2T[[il]] ] /. Tensor2Class ] ];
  tmp = Rule @@@ Partition[ Riffle[ ind /. wcname -> ToExpression[ SMEFT$MB <> Active2T[[il]]], WC], 2 ];
  AppendTo[Tensor2Val, Active2T[[il]] -> tmp];
];

ind = Flatten[ Table[wcname[i,j,k,l], {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}] ];
ind1 = Table["_"<>ToString[i]<>ToString[j]<>ToString[k]<>ToString[l], {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}];
Active4T = ToString /@ Intersection[Tensor4WC[[All,1]], oplist];
Tensor4Val = {};
For[ il = 1, il < Length[Active4T] + 1, il++,
  WC4 = Table[ JSONNum[ JSONWilson[ (Active4T[[il]] /. Tensor4JSON) <> ind1[[ig1,ig2,ig3,ig4]] ] ],
               {ig1,1,3},{ig2,1,3},{ig3,1,3},{ig4,1,3} ];
  WC4 = Flatten[ Tensor4WCSymmetrize[ WC4, ToExpression[ SMEFT$MB <>  Active4T[[il]] ] /. Tensor4Class ] ];
  tmp = Rule @@@ Partition[ Riffle[ ind /. wcname -> ToExpression[ SMEFT$MB <> Active4T[[il]]], WC4], 2 ];
  AppendTo[Tensor4Val, Active4T[[il]] -> tmp];
];

(* store all WC value substitutions in public variables *)
SMEFT$WCValues = Flatten[ Join[ ScalarVal, Values[Tensor2Val], Values[Tensor4Val], {Lam->1} ] ];
SMEFT$ScalarVal  = ScalarVal;
SMEFT$Tensor2Val = Tensor2Val;
SMEFT$Tensor4Val = Tensor4Val;
SMEFT$WCXFSource = SourceFile;

]
(* end of ReadWCXFInput *)





WCXFToSMEFT[ TargetFile_, OptionsPattern[ { OverwriteTarget -> False,
					    RealParameters -> False,
                                            Silent -> False} ] ] :=
(* generates SMEFT model file from variables generated by reading WCxf input file*)
Module[{date,cfile,TargetExistsQ,a},

If[ Head[TargetFile] =!= String, Print["WCXFToSMEFT: Target file name must be string!"]; Abort[]; ];

If [ FileExistsQ[ TargetFile ] && (! OptionValue[OverwriteTarget]),
    Print["File " <> TargetFile <> " exists!"];
    TargetExistsQ = AskFunction[ Ask["Overwrite (y/n)?"-> "String"] ][];
];
If[ TargetExistsQ != "y", Print[ "Change name of parameter file generated from WCXF format and rerun" ]; Abort[]; ];
If [ OptionValue[Silent] =!= True,
  Print["Parameter file in mass basis generated as ",  Style[TargetFile, Bold ] ];
];

If[ OptionValue[RealParameters],
  Print[Style["WARNING: Only real parts of SMEFT couplings are initialized, as accepted by MADGRAPH5; use option RealParameters->False in SMEFTInitializeModel to allow complex parameters in SMEFT model file",Bold]];
,
  Print[Style["WARNING: SMEFT couplings may have imaginary parts, not accepted by MADGRAPH5; use option RealParameters->True in SMEFTInitializeModel to force storing only parts of parameters in SMEFT model file",Bold]]
];

cfile = OpenWrite[TargetFile];

WriteLine[cfile, "(* Generated by WCXFToSMEFT routine from the file " <> SMEFT$WCXFSource <> " *)"];
date = Date[];
WriteLine[cfile, "(* " <> ToString[date[[3]]] <> "." <>
    ToString[date[[2]]] <> "." <> ToString[date[[1]]] <> " " <>
    ToString[date[[4]]] <> ":" <> ToString[date[[5]]] <> " *)\n\n"]
WriteLine[cfile, "(* Active operators included in Feynman Rules: " <> ToString[SMEFT$OperatorList] <> " *)\n"];
       
WriteString[cfile, ReadString[FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_head_MB.fr"}]] ];
WriteLine[cfile, "\n"];

(* add input observables in preferred scheme to fix SM Lagrangian
   parameters; see file smeft_input_parameters.m *)
WriteLine[cfile, "(* input observables *)\n"];
WriteInputParameter[cfile,#] & /@ SM$InputParameters;

(* standard SMEFT parameters and their numerical values *)
WriteLine[cfile, "(* standard SMEFT parametrization *)\n"];
WriteInputParameter[cfile,#] & /@ If [ OptionValue[RealParameters], SMEFT$SMParametersReal, SMEFT$SMParameters ];

WriteLine[cfile, "(* redefined (mass basis) WC coefficients *)\n"];
WriteLine[cfile, "(* flavor independent *)\n"];
WriteExternalScalarEntry[ cfile, # ] & /@ If[ OptionValue[RealParameters],
                                              Map[Re, SMEFT$ScalarVal, {2}]  /. Re[a_] -> a,
					      SMEFT$ScalarVal ];

WriteLine[cfile, "(* flavor dependent *)\n"];
WriteLine[cfile, "(* 2 fermion operators *)\n"];
WriteExternalTensor2Entry[ cfile, # ] & /@ If[ OptionValue[RealParameters],
					       Map[Re, SMEFT$Tensor2Val, {4}]  /. Re[a_] -> a,
					       SMEFT$Tensor2Val ] ;
WriteLine[cfile, "(* 4 fermion operators *)\n"];
WriteExternalTensor4Entry[ cfile, # ] & /@ If[ OptionValue[RealParameters],
                                               Map[Re, SMEFT$Tensor4Val, {4}]  /. Re[a_] -> a,
					       SMEFT$Tensor4Val ];

WriteLine[cfile, "(* Effective NP scale squared*)\n"];
WriteLine[cfile, "  Lam == {"];
WriteLine[cfile, "    ParameterType -> External,"];
(*
WriteLine[cfile, "    Value   -> " <> NumberToString[ ("scale" /. JSONData)^2 ] <> ","];
*)
WriteLine[cfile, "    Value   -> 1,"];
WriteLine[cfile, "    TeX           -> 1/\[CapitalLambda]^2,"];
WriteLine[cfile, "    Description   -> \"Effective NP scale squared\""];
WriteLine[cfile, "  }\n\n"];

WriteLine[cfile, "}\n"];

(* updates the correct maximal NP interaction order for model files *)
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_io.fr"}] ];
WriteLine[cfile, "M$InteractionOrderHierarchy = ", M$InteractionOrderHierarchy, "\n"];
WriteLine[cfile, "M$InteractionOrderLimit = ", M$InteractionOrderLimit /. {NP,a_} -> {NP,SMEFT$ExpansionOrder}, "\n"];


Close[TargetFile];

]
(* end of WCXFToSMEFT *)


